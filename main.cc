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
#include <openssl/rand.h>

#define BITMASK(nbits) ((nbits) == 64 ? 0xffffffffffffffff : (1ULL << (nbits)) \
												- 1ULL)
#define QBITS_LOCAL_QF 16
#define SPARE_EMPTY_LOCAL_QFS 16

using namespace std;
using namespace kmercounting;

typedef struct {
	QF *local_qf;
	QF *main_qf;
	uint32_t count {0};
	uint32_t ksize {28};
}flush_object;

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
	cout << desc << "Total Time Elapsed: " << to_string(time_elapsed) <<
		"seconds" << endl;
}

/* dump the contents of a local QF into the main QF */
static void dump_local_qf_to_main(flush_object *obj)
{
	QFi local_cfi;

	if (qf_iterator(obj->local_qf, &local_cfi, 0)) {
		do {
			uint64_t key = 0, value = 0, count = 0;
			qfi_get(&local_cfi, &key, &value, &count);
			qf_insert(obj->main_qf, key, 0, count, true, true);
		} while (!qfi_next(&local_cfi));
		qf_reset(obj->local_qf);
	}
}
/* main method */
int main(int argc, char *argv[])
{
	QF cf;
	QFi cfi;
	printf("Testing cqf\n");

	uint32_t seed = 2038074761;
		//Initialize the main  QF
	uint64_t qbits=16;
	uint64_t num_hash_bits=qbits+8;
	uint64_t nSlots=(1ULL<<qbits);
	uint64_t nValues=(uint64_t)((double)nSlots*0.20);
	uint64_t *vals;
	bool mem=false;


	qf_init(&cf, nSlots, num_hash_bits, 0, true, "", seed);
	printf("cqf initialized nhashbots=%d\n",num_hash_bits);
	vals = (uint64_t*)malloc(nValues*sizeof(vals[0]));
	RAND_pseudo_bytes((unsigned char *)vals, sizeof(*vals) * nValues);
	for (uint64_t i = 0; i < nValues; i++) {
			vals[i] = (1 * vals[i]);
		}
	for (uint64_t i = 0; i < nValues; i++) {
			  vals[i]=vals[i]%cf.metadata->range;
				qf_insert(&cf, vals[i], 0, 50,false,true);
	}
	printf("%d are inserted\n",nValues);

	QF cfdisk;
	uint64_t qbitsDisk=19;
	uint64_t num_hash_bitsDisk=qbitsDisk+5;
	uint64_t nSlotsDisk=(1ULL<<qbitsDisk);


	qf_init(&cfdisk, nSlotsDisk, num_hash_bitsDisk, 0, mem, "/home/mostafa/Documents/squeakr/cqf.mmap", seed);
	printf("cqf initialized nhashbots=%d\n",num_hash_bitsDisk);


	/* Initialize an iterator */
	qf_iterator(&cf, &cfi, 0);
	do {
		uint64_t key, value, count;
		qfi_get(&cfi, &key, &value, &count);
		qf_insert(&cfdisk, key, 0, count,false,true);
	} while(!qfi_next(&cfi));


	for (uint64_t i = 0; i < nValues; i++) {
		uint64_t count = qf_count_key_value(&cfdisk, vals[i], 0);
		if (count < 50) {
			fprintf(stderr, "failed lookup after insertion for %lx %ld.\n", vals[i],
							count);
			abort();
		}
	}

	/* Initialize an iterator */
	qf_iterator(&cfdisk, &cfi, 0);
	do {
		uint64_t key, value, count;
		qfi_get(&cfi, &key, &value, &count);
		if (qf_count_key_value(&cfdisk, key, 0) < 50) {
			fprintf(stderr, "Failed lookup from A for: %ld. Returned count: %ld\n",
							key, qf_count_key_value(&cf, key, 0));
			abort();
		}
	} while(!qfi_next(&cfi));

	fprintf(stdout, "Validated the CQF.\n");


	//destroy the QF and reclaim the memory
	qf_destroy(&cf, true);
	qf_destroy(&cfdisk, false);
	return 0;
}
