/*
 * ============================================================================
 *
 *       Filename:  main.cc
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  05/06/2016 09:56:26 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Prashant Pandey (ppandey@cs.stonybrook.edu)
 *                  Rob Patro (rob.patro@cs.stonybrook.edu)
 *                  Rob Johnson (rob@cs.stonybrook.edu)
 *   Organization:  Stony Brook University
 *
 * ============================================================================
 */

#include <iostream>
#include <algorithm>
#include <cstring>
#include <vector>
#include <set>
#include <unordered_set>
#include <bitset>
#include <cassert>
#include <fstream>

#include <boost/thread/thread.hpp>
#include <boost/lockfree/queue.hpp>
#include <boost/lockfree/spsc_queue.hpp>
#include <boost/atomic.hpp>

#include <time.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>

#include "clipp.h"
#include "threadsafe-gqf/gqf.h"
#include "hashutil.h"
#include "chunk.h"
#include "kmer.h"
#include "reader.h"

#define BITMASK(nbits) ((nbits) == 64 ? 0xffffffffffffffff : (1ULL << (nbits)) \
												- 1ULL)
#define QBITS_LOCAL_QF 16
#define SPARE_EMPTY_LOCAL_QFS 16

using namespace std;
using namespace kmercounting;

typedef struct {
	QF *local_qf;
	QF *main_qf;
	QF *main_qfDisk;
	uint32_t count {0};
	uint32_t ksize {28};

}flush_object;

volatile bool main_qf_lock=false;
volatile uint64_t main_qf_count=0;

bool exactCounting=false;

 #define Max_Main_QF_Load_Factor 0.9

struct file_pointer {
	std::unique_ptr<reader> freader{nullptr};
	char* part{nullptr};
	char* part_buffer{nullptr};
	int mode{0};
	uint64_t size{0};
	uint64_t part_filled{0};
};

/*create a multi-prod multi-cons queue for storing the chunk of fastq file.*/
boost::lockfree::queue<file_pointer*, boost::lockfree::fixed_sized<true> > ip_files(64);
boost::atomic<int> num_files {0};

/* Count distinct items in a sorted list */
uint64_t count_distinct_kmers(multiset<uint64_t> kmers)
{
	uint64_t cnt = 0;
	uint64_t curr_kmer = 0;

	for(uint64_t kmer: kmers) {
		if (kmer != curr_kmer) {
			curr_kmer = kmer;
			cnt++;
		}
	}
	return cnt;
}

/* Print elapsed time using the start and end timeval */
void print_time_elapsed(string desc, struct timeval* start, struct timeval* end)
{
	struct timeval elapsed;
	if (start->tv_usec > end->tv_usec) {
		end->tv_usec += 1000000;
		end->tv_sec--;
	}
	elapsed.tv_usec = end->tv_usec - start->tv_usec;
	elapsed.tv_sec = end->tv_sec - start->tv_sec;
	float time_elapsed = (elapsed.tv_sec * 1000000 + elapsed.tv_usec)/1000000.f;
	cout << desc << "\t" << to_string(time_elapsed) <<
		" seconds" << endl;
}

/* check if it's the end of the file. */
inline bool is_eof(reader &file_reader, int mode)
{
	if (mode == 0)
		return feof(file_reader.in) != 0;
	else if (mode == 1)
		return gzeof(file_reader.in_gzip) != 0;
	else if (mode == 2)
		return file_reader.bzerror == BZ_STREAM_END;

	return true;
}

/* move the pointer to the end of the next newline. */
bool skip_next_eol(char *part, int64_t &pos, int64_t max_pos)
{
	int64_t i;
	for(i = pos; i < max_pos-2; ++i)
		if((part[i] == '\n' || part[i] == '\r') && !(part[i+1] == '\n' ||
																								 part[i+1] == '\r'))
			break;

	if(i >= max_pos-2)
		return false;
	pos = i+1;

	return true;
}

/* read a part of the fastq file. */
static bool fastq_read_parts(int mode, file_pointer *fp)
{
	char *& _part = (fp->part);
	uint64_t& _size = fp->size;
	char*& part_buffer = (fp->part_buffer);
	uint64_t& part_filled = fp->part_filled;
	reader& file_reader = *(fp->freader.get());

	uint32_t OVERHEAD_SIZE = 65535;
	uint64_t part_size = 1ULL << 23;
	char *part = (char *)malloc((part_size + OVERHEAD_SIZE)*sizeof(char));
	memcpy(part, part_buffer, part_filled);

	if(is_eof(file_reader, mode))
		return false;

	uint64_t readed = 0;

	if (mode == 0)
		readed = fread(part+part_filled, 1, part_size, file_reader.in);
	else if (mode == 1)
		readed = gzread(file_reader.in_gzip, part+part_filled, (int) part_size);
	else if (mode == 2)
		readed = BZ2_bzRead(&file_reader.bzerror, file_reader.in_bzip2,
												part+part_filled, (int) part_size);
	else
		readed = 0;

	int64_t total_filled = part_filled + readed;
	int64_t i;
	if(part_filled >= OVERHEAD_SIZE)
	{
		cout << "Error: Wrong input file!\n";
		exit(EXIT_FAILURE);
	}
	if(is_eof(file_reader, mode))
	{
		_part = part;
		_size = total_filled;
		part = NULL;
		return true;
	}
	// Looking for a FASTQ record at the end of the area
	{
		int64_t line_start[9];
		int32_t j;
		i = total_filled - OVERHEAD_SIZE / 2;
		for(j = 0; j < 9; ++j)
		{
			if(!skip_next_eol(part, i, total_filled))
				break;
			line_start[j] = i;
		}
		_part = part;
		if(j < 9)
			_size = 0;
		else
		{
			int k;
			for(k = 0; k < 4; ++k)
			{
				if(part[line_start[k]+0] == '@' && part[line_start[k+2]+0] == '+')
				{
					if(part[line_start[k+2]+1] == '\n' || part[line_start[k+2]+1] == '\r')
						break;
					if(line_start[k+1]-line_start[k] == line_start[k+3]-line_start[k+2] &&
						 memcmp(part+line_start[k]+1, part+line_start[k+2]+1,
										line_start[k+3]-line_start[k+2]-1) == 0)
						break;
				}
			}
			if(k == 4)
				_size = 0;
			else
				_size = line_start[k];
		}
	}

	copy(_part+_size, _part+total_filled, part_buffer);
	part_filled = total_filled - _size;

	return true;
}
static inline bool qf_spin_lock(volatile int *lock, bool flag_spin)
{
	if (!flag_spin) {
		return !__sync_lock_test_and_set(lock, 1);
	} else {
		while (__sync_lock_test_and_set(lock, 1))
			while (*lock);
		return true;
	}

	return false;
}
uint64_t nMerges=0;
/* dump the contents of  the main QF  int Disk QF*/
static void dump_main_qf_to_disk(flush_object *obj)
{
	qf_spin_lock((int*)&main_qf_lock,true);
	main_qf_lock=true;
	nMerges++;
	QFi cfi;
	if (qf_iterator(obj->main_qf, &cfi, 0)) {

		do {
			uint64_t key = 0, value = 0, count = 0;
			qfi_get(&cfi, &key, &value, &count);
			qf_insert(obj->main_qfDisk, key, 0, count, true, true);
		} while (!qfi_next(&cfi));
		qf_reset(obj->main_qf);
		main_qf_count=0;
	}

	main_qf_lock=false;
}


/* dump the contents of a local QF into the main QF */
static void dump_local_qf_to_main(flush_object *obj)
{
	QFi local_cfi;

	if (qf_iterator(obj->local_qf, &local_cfi, 0)) {
		do {
			uint64_t key = 0, value = 0, count = 0;
			qfi_get(&local_cfi, &key, &value, &count);
			qf_spin_lock((int*)&main_qf_lock,true);
			main_qf_lock=false;
			qf_insert(obj->main_qf, key, 0, count, true, true);
			main_qf_count++;
			double loadFactor=(double)obj->main_qf->metadata->noccupied_slots/
																		(double)obj->main_qf->metadata->nslots;
			if(loadFactor>Max_Main_QF_Load_Factor){
					dump_main_qf_to_disk(obj);
			}

		} while (!qfi_next(&local_cfi));
		qf_reset(obj->local_qf);
	}
}


/* convert a chunk of the fastq file into kmers */
void reads_to_kmers(chunk &c, flush_object *obj)
{
	auto fs = c.get_reads();
	auto fe = c.get_reads();
	auto end = fs + c.get_size();
	while (fs && fs!=end) {
		fs = static_cast<char*>(memchr(fs, '\n', end-fs)); // ignore the first line
		fs++; // increment the pointer

		fe = static_cast<char*>(memchr(fs, '\n', end-fs)); // read the read
		string read(fs, fe-fs);
		/*cout << read << endl;*/

start_read:
		if (read.length() < obj->ksize) // start with the next read if length is smaller than K
			goto next_read;
		{
			uint64_t first = 0;
			uint64_t first_rev = 0;
			uint64_t item = 0;
			//cout << "K " << read.substr(0,K) << endl;
			for(int i=0; i<obj->ksize; i++) { //First kmer
				uint8_t curr = kmer::map_base(read[i]);
				if (curr > DNA_MAP::G) { // 'N' is encountered
					read = read.substr(i+1, read.length());
					goto start_read;
				}
				first = first | curr;
				first = first << 2;
			}
			first = first >> 2;
			first_rev = kmer::reverse_complement(first, obj->ksize);

			//cout << "kmer: "; cout << int_to_str(first);
			//cout << " reverse-comp: "; cout << int_to_str(first_rev) << endl;

			if (kmer::compare_kmers(first, first_rev))
				item = first;
			else
				item = first_rev;

			// hash the kmer using murmurhash/xxHash before adding to the list
			//item = HashUtil::MurmurHash64A(((void*)&item), sizeof(item),
				//														 obj->local_qf->metadata->seed);

			if(exactCounting)
			{
				item = HashUtil::hash_64(item, BITMASK(2*obj->ksize));
			}
			else{
				item = HashUtil::MurmurHash64A(((void*)&item), sizeof(item),
																				 obj->local_qf->metadata->seed);
			}

			/*
			 * first try and insert in the main QF.
			 * If lock can't be accuired in the first attempt then
			 * insert the item in the local QF.
			 */
			if (!main_qf_lock && !qf_insert(obj->main_qf, item%obj->main_qf->metadata->range, 0, 1,
										 true, false)) {
				qf_insert(obj->local_qf, item%obj->local_qf->metadata->range, 0, 1,
									false, false);
				obj->count++;
				// check of the load factor of the local QF is more than 50%
				if (obj->count > 1ULL<<(QBITS_LOCAL_QF-1)) {
					dump_local_qf_to_main(obj);
					obj->count = 0;
				}
			}
			else if(!main_qf_lock){
					// kmer is inserted to main qf.
					// Check if the main qf(memory) is full and dump it
					main_qf_count++;
					double loadFactor=(double)obj->main_qf->metadata->noccupied_slots/
																						(double)obj->main_qf->metadata->nslots;
					if(loadFactor>Max_Main_QF_Load_Factor){
						dump_main_qf_to_disk(obj);
					}
			}
			//cout<< "X " << bitset<64>(first)<<endl;

			uint64_t next = (first << 2) & BITMASK(2*obj->ksize);
			uint64_t next_rev = first_rev >> 2;

			for(uint32_t i=obj->ksize; i<read.length(); i++) { //next kmers
				//cout << "K: " << read.substr(i-K+1,K) << endl;
				uint8_t curr = kmer::map_base(read[i]);
				if (curr > DNA_MAP::G) { // 'N' is encountered
					read = read.substr(i+1, read.length());
					goto start_read;
				}
				next |= curr;
				uint64_t tmp = kmer::reverse_complement_base(curr);
				tmp <<= (obj->ksize*2-2);
				next_rev = next_rev | tmp;
				if (kmer::compare_kmers(next, next_rev))
					item = next;
				else
					item = next_rev;

			// hash the kmer using murmurhash/xxHash before adding to the list
				//item = HashUtil::MurmurHash64A(((void*)&item), sizeof(item),
				//															 obj->local_qf->metadata->seed);
				//item = XXH63 (((void*)&item), sizeof(item), seed);

				if(exactCounting)
				{
					item = HashUtil::hash_64(item, BITMASK(2*obj->ksize));
				}
				else{
					item = HashUtil::MurmurHash64A(((void*)&item), sizeof(item),
																				 obj->local_qf->metadata->seed);
				}

				/*
				 * first try and insert in the main QF.
				 * If lock can't be accuired in the first attempt then
				 * insert the item in the local QF.
				 */
				if (!main_qf_lock && !qf_insert(obj->main_qf, item%obj->main_qf->metadata->range, 0, 1, true,
											 false)) {
					qf_insert(obj->local_qf, item%obj->local_qf->metadata->range, 0, 1, false,
										false);
					obj->count++;
					// check of the load factor of the local QF is more than 50%
					if (obj->count > 1ULL<<(QBITS_LOCAL_QF-1)) {
						dump_local_qf_to_main(obj);
						obj->count = 0;
					}
				}
				else if(!main_qf_lock){
						// kmer is inserted to main qf.
						// Check if the main qf(memory) is full and dump it
						main_qf_count++;
						double loadFactor=(double)obj->main_qf->metadata->noccupied_slots/
																							(double)obj->main_qf->metadata->nslots;

						if(loadFactor>Max_Main_QF_Load_Factor){
							dump_main_qf_to_disk(obj);
						}
				}

				//cout<<bitset<64>(next)<<endl;
				//assert(next == str_to_int(read.substr(i-K+1,K)));

				next = (next << 2) & BITMASK(2*obj->ksize);
				next_rev = next_rev >> 2;
			}
		}

next_read:
		fs = ++fe;		// increment the pointer
		fs = static_cast<char*>(memchr(fs, '\n', end-fs)); // ignore one line
		fs++; // increment the pointer
		fs = static_cast<char*>(memchr(fs, '\n', end-fs)); // ignore one more line
		fs++; // increment the pointer
	}
	free(c.get_reads());
}

/* read a part of the fastq file, parse it, convert the reads to kmers, and
 * insert them in the CQF
 */
static bool fastq_to_uint64kmers_prod(flush_object* obj)
{
	file_pointer* fp;

	while (num_files) {
		while (ip_files.pop(fp)) {
			if (fastq_read_parts(fp->mode, fp)) {
				ip_files.push(fp);
				chunk c(fp->part, fp->size);
				reads_to_kmers(c, obj);
			} else {
				/* close the file */
				if (fp->mode == 0)
					fclose(fp->freader->in);
				else if (fp->mode == 1)
					gzclose(fp->freader->in_gzip);
				else if (fp->mode == 2)
					if (fp->freader->in) {
						BZ2_bzReadClose(&(fp->freader->bzerror), fp->freader->in_bzip2);
						fclose(fp->freader->in);
					}
				delete[] fp->part_buffer;
				delete fp;
				num_files--;
			}
		}
	}
	if (obj->count) {
		dump_local_qf_to_main(obj);
		obj->count = 0;
	}

	return true;
}

bool getFileReader(int mode, const char* fastq_file, reader* file_reader)
{
	uint64_t gzip_buffer_size = 1ULL << 26;
	uint64_t bzip2_buffer_size = 1ULL << 26;

	if (mode == 0) {
		if ((file_reader->in = fopen(fastq_file, "rb")) == NULL)
			return false;
	} else if (mode == 1) {
		if ((file_reader->in_gzip = gzopen(fastq_file, "rb")) == NULL)
			return false;
		gzbuffer(file_reader->in_gzip, gzip_buffer_size);
	} else if (mode == 2) {
		file_reader->in = fopen(fastq_file, "rb");
		if (!file_reader->in)
			return false;
		setvbuf(file_reader->in, NULL, _IOFBF, bzip2_buffer_size);
		if ((file_reader->in_bzip2 = BZ2_bzReadOpen(&file_reader->bzerror,
																								file_reader->in, 0, 0, NULL,
																								0)) == NULL) {
			fclose(file_reader->in);
			return false;
		}
	}
	return true;
}

/* main method */
int main(int argc, char *argv[])
{


  enum class file_type {fastq, gzip, bzip2};

  file_type in_type = file_type::fastq;
  int mode = 0;
  int ksize;
  int qbits;
	int qbitsM;
  int numthreads;
  std::string prefix = "./";
  std::vector<std::string> filenames;
  using namespace clipp;
  auto cli = (
              one_of(
                      required("-f").set(in_type, file_type::fastq) % "plain fastq",
                      required("-g").set(in_type, file_type::gzip) % "gzip compressed fastq",
                      required("-b").set(in_type, file_type::bzip2) % "bzip2 compressed fastq") % "format of the input",
              required("-k","--kmer") & value("k-size", ksize) % "length of k-mers to count",
							option("-e").set(exactCounting, true) % "exact counting",
              required("-s","--log-slots") & value("log-slots", qbits) % "log of number of slots in the CQF on disk",
							required("-m","--log-slots-memory") & value("log-slots-memory", qbitsM) % "log of number of slots in the CQF on the memory",
              required("-t","--threads") & value("num-threads", numthreads) % "number of threads to use to count",
              option("-o","--output-dir") & value("out-dir", prefix) % "directory where output should be written (default = \"./\")",
              values("files", filenames) % "list of files to be counted",
              option("-h", "--help")      % "show help"
              );

  auto res = parse(argc, argv, cli);

  if (!res) {
    std::cout << make_man_page(cli, argv[0]) << "\n";
    return 1;
  }

	if (exactCounting && ksize > 32) {
    cerr << "k-mer size greater than 32 is not supported for exact counting." << endl;
    return 1;
  }


  switch (in_type) {
  case file_type::fastq: mode = 0; break;
  case file_type::gzip: mode = 1; break;
  case file_type::bzip2: mode = 2; break;
  }
  /*
  if (argc == 2) {
		string arg_help(argv[1]);
		if (arg_help.compare("-h") != 0 || arg_help.compare("-help") != 0) {
			cout << "./squeakr-count [OPTIONS]" << endl
				   << "file format   : 0 - plain fastq, 1 - gzip compressed fastq, 2 - bzip2 compressed fastq" << endl
					 << "CQF size      : the log of the number of slots in the CQF" << endl
					 << "num of threads: number of threads to count" << endl
					 << "file(s)       : \"filename\" or \"dirname/\" for all the files in a directory" << endl;
			exit(0);
		}
	}
  */

  /*
	int mode = atoi(argv[1]);
	int ksize = atoi(argv[2]);
	int qbits = atoi(argv[3]);
	int numthreads = atoi(argv[4]);
  */

	string ser_ext(".ser");
	string log_ext(".log");
	string cluster_ext(".cluster");
	string freq_ext(".freq");
	struct timeval start1, start2, end1, end2;
	struct timezone tzp;
	uint32_t OVERHEAD_SIZE = 65535;

  for( auto& fn : filenames ) {
		auto* fr = new reader;
		if (getFileReader(mode, fn.c_str(), fr)) {
			file_pointer* fp = new file_pointer;
			fp->mode = mode;
			fp->freader.reset(fr);
			fp->part_buffer = new char[OVERHEAD_SIZE];
			ip_files.push(fp);
			num_files++;
		} else {
			delete fr;
		}
	}

  string filepath(filenames.front());
  auto const pos = filepath.find_last_of('/');
  string filename = filepath.substr(pos+1);
  if (prefix.back() != '/') {
    prefix += '/';
  }
	string ds_file =      prefix + filename + ser_ext;
	string log_file =     prefix + filename + log_ext;
	string cluster_file = prefix + filename + cluster_ext;
	string freq_file_name =    prefix + filename + freq_ext;

	QF cf;
	QF cfM;
	QFi cfi;
	QF local_qfs[50];
	uint32_t seed = 2038074761;
	int num_hash_bits;
	if(exactCounting){
			num_hash_bits = 2*ksize; // Each base 2 bits.
	}
	else{
			num_hash_bits=qbits+8;
	}
	//Initialize the main  QF
	qf_init(&cf, (1ULL<<qbits), num_hash_bits, 0, false, ds_file.c_str(), seed);
	qf_init(&cfM, (1ULL<<qbitsM), num_hash_bits, 0, true, "", seed);

	boost::thread_group prod_threads;
	flush_object* obj;
	for (int i = 0; i < numthreads; i++) {
		qf_init(&local_qfs[i], (1ULL << QBITS_LOCAL_QF), num_hash_bits, 0, true,
						"", seed);
		obj = (flush_object*)malloc(sizeof(flush_object));
		obj->local_qf = &local_qfs[i];
		obj->main_qf = &cfM;
		obj->main_qfDisk=&cf;
		obj->ksize = ksize;
		prod_threads.add_thread(new boost::thread(fastq_to_uint64kmers_prod,
																							obj));
	}

	cout << "Reading from the fastq file and inserting in the QF" << endl;
	gettimeofday(&start1, &tzp);
	prod_threads.join_all();
  dump_main_qf_to_disk(obj);
	//qf_serialize(&cf, ds_file.c_str()); no need because I am using mmap
	gettimeofday(&end1, &tzp);
	print_time_elapsed("InsertionTime", &start1, &end1);

	qf_destroy(&cfM, true);
	cout<<"nMerges\t"<<nMerges<<endl;

	cout << "Calc freq distribution: " << endl;
	ofstream freq_file;
	freq_file.open(freq_file_name.c_str());
	uint64_t max_cnt = 0;
	qf_iterator(&cf, &cfi, 0);
	gettimeofday(&start2, &tzp);
	do {
		uint64_t key = 0, value = 0, count = 0;
		qfi_get(&cfi, &key, &value, &count);
		if(exactCounting){
			freq_file << int_to_str(HashUtil::hash_64i(key,BITMASK(2*ksize)),ksize) << " " << count << endl;
		}
		else{
			freq_file << key << " " << count << endl;
		}
		if (max_cnt < count)
			max_cnt = count;
	} while (!qfi_next(&cfi));
	gettimeofday(&end2, &tzp);
	print_time_elapsed("QueryTime", &start2, &end2);

	cout << "Maximum freq: " << max_cnt << endl;
	freq_file.close();

	cout << "Num distinct elem: " << cf.metadata->ndistinct_elts << endl;
	cout << "Total num elems: " << cf.metadata->nelts << endl;

#ifdef LOG_WAIT_TIME
	ofstream wait_time_log;
	wait_time_log.open(log_file.c_str());
	wait_time_log << "Id\tTotalTimeSingle\tTotalTimeSpinning\tNumLocks\tNumSingleAttempt\tPercentageSpinningTimes"
		<< endl;
	for (uint32_t i=0; i<cf.num_locks; i++)
		wait_time_log << i << "\t" << cf.wait_times[i].total_time_single << "\t\t\t"
			<< cf.wait_times[i].total_time_spinning << "\t\t\t"
			<< cf.wait_times[i].locks_taken << "\t\t\t"
			<< cf.wait_times[i].locks_acquired_single_attempt << "\t\t\t"
			<< ((double)(cf.wait_times[i].locks_taken
									 -cf.wait_times[i].locks_acquired_single_attempt)
					/(double)cf.wait_times[i].locks_taken)*100
			<< endl;
	wait_time_log.close();
#endif

#ifdef LOG_CLUSTER_LENGTH
	ofstream cluster_len_log;
	cluster_len_log.open(cluster_file.c_str());
	cluster_len_log << "StartingIndex\tLength" << endl;
	for (uint32_t i = 0; i < cfi.num_clusters; i++)
		cluster_len_log << cfi.c_info[i].start_index << "\t\t" <<
			cfi.c_info[i].length << endl;
	cluster_len_log.close();

#endif

	//destroy the QF and reclaim the memory
	qf_destroy(&cf, false);

	return 0;
}
