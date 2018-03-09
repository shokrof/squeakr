#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include<iostream>
using namespace std;

TEST_CASE( "get/set fixed counters" ,"[!hide]") {
  QF qf;
  for(int counter_size=1;counter_size<=5;counter_size++){
    uint64_t qbits=15;
    uint64_t num_hash_bits=qbits+8;
    uint64_t maximum_count=(1ULL<<counter_size)-1;
    INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
    qf_init(&qf, (1ULL<<qbits), num_hash_bits, 0,counter_size, true, "", 2038074761);
    uint64_t c;
    // test many set and get
    uint64_t last=0;
    for(int i=1;i<=maximum_count;i++){
      REQUIRE( get_fixed_counter(&qf,99) == 0 );
      CHECK( get_fixed_counter(&qf,100) == last );
      REQUIRE( get_fixed_counter(&qf,101) == 0 );
      set_fixed_counter(&qf,100,i);
      REQUIRE( get_fixed_counter(&qf,99) == 0 );
      CHECK( get_fixed_counter(&qf,100) == i );
      REQUIRE( get_fixed_counter(&qf,101) == 0 );
      last=i;
    }

    //test on block boundaries

    REQUIRE( get_fixed_counter(&qf,63) == 0 );
    REQUIRE( get_fixed_counter(&qf,64) == 0 );
    REQUIRE( get_fixed_counter(&qf,65) == 0 );

    c=1;
    set_fixed_counter(&qf,63,c);
    c=(c+1)%maximum_count;
    set_fixed_counter(&qf,64,c);
    c=(c+1)%maximum_count;
    set_fixed_counter(&qf,65,c);

    c=1;
    REQUIRE( get_fixed_counter(&qf,63) == c );
    c=(c+1)%maximum_count;
    REQUIRE( get_fixed_counter(&qf,64) == c );
    c=(c+1)%maximum_count;
    REQUIRE( get_fixed_counter(&qf,65) == c );

    // test special slots
    c=1;
    uint64_t special_slots[]={
      0,
      (1ULL<qbits),
      qf.metadata->xnslots-1
    };
    for(int i=0;i<3;i++){
      INFO("Testing Special Slot "<<special_slots[i]);
      REQUIRE( get_fixed_counter(&qf,special_slots[i]) == 0 );
      set_fixed_counter(&qf,special_slots[i],c);
      REQUIRE( get_fixed_counter(&qf,special_slots[i]) == c );
    }

    qf_destroy(&qf,true);
  }
}

TEST_CASE( "shift fixed counters","[!hide]" ) {
  QF qf;
  for(int counter_size=1;counter_size<=5;counter_size++){
    uint64_t qbits=15;
    uint64_t num_hash_bits=qbits+8;
    uint64_t maximum_count=(1ULL<<counter_size)-1;
    INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
    qf_init(&qf, (1ULL<<qbits), num_hash_bits, 0,counter_size, true, "", 2038074761);
    uint64_t c;
    // test  shift one item
    uint64_t last=0;
    set_fixed_counter(&qf,100,1);
    CHECK( get_fixed_counter(&qf,100) == 1 );
    shift_fixed_counters(&qf,100,100,1);
    CHECK( get_fixed_counter(&qf,100) == 0 );
    CHECK( get_fixed_counter(&qf,101) == 1 );

    //test shift many items
    set_fixed_counter(&qf,100,maximum_count);
    set_fixed_counter(&qf,99,maximum_count);
    shift_fixed_counters(&qf,99,101,5);

    CHECK( get_fixed_counter(&qf,99) == 0 );
    CHECK( get_fixed_counter(&qf,100) == 0 );
    CHECK( get_fixed_counter(&qf,101) == 0 );

    CHECK( get_fixed_counter(&qf,104) == maximum_count );
    CHECK( get_fixed_counter(&qf,105) == maximum_count );
    CHECK( get_fixed_counter(&qf,106) == 1 );

    // test shift on block boundary
    set_fixed_counter(&qf,63,maximum_count);
    set_fixed_counter(&qf,64,maximum_count);
    set_fixed_counter(&qf,65,maximum_count);

    shift_fixed_counters(&qf,63,65,5);

    CHECK( get_fixed_counter(&qf,63) == 0 );
    CHECK( get_fixed_counter(&qf,64) == 0 );
    CHECK( get_fixed_counter(&qf,65) == 0 );

    CHECK( get_fixed_counter(&qf,68) == maximum_count );
    CHECK( get_fixed_counter(&qf,69) == maximum_count );
    CHECK( get_fixed_counter(&qf,70) == maximum_count );

    qf_destroy(&qf,true);
  }
}


TEST_CASE( "Add fixed counters to items","[!hide]" ) {
  QF qf;
  for(int counter_size=1;counter_size<=5;counter_size++){
    uint64_t qbits=15;
    uint64_t num_hash_bits=qbits+8;
    uint64_t maximum_count=(1ULL<<counter_size)-1;
    INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
    qf_init(&qf, (1ULL<<qbits), num_hash_bits, 0,counter_size, true, "", 2038074761);

    qf_insert(&qf,150,0,50,false,false);
    CHECK( qf_count_key_value(&qf,150,0)==50);
    qf_set_fixed_counter(&qf,150,maximum_count);
    CHECK( qf_get_fixed_counter(&qf,150)==maximum_count);
    CHECK( qf_count_key_value(&qf,150,0)==50);

    for(int i=120;i<=149;i++){
      qf_insert(&qf,i,0,50,false,false);
      CHECK( qf_count_key_value(&qf,i,0)==50);
      qf_set_fixed_counter(&qf,i,1);
      CHECK( qf_get_fixed_counter(&qf,i)==1);
      CHECK( qf_count_key_value(&qf,i,0)==50);
    }
    CHECK( qf_get_fixed_counter(&qf,150)==maximum_count);
    CHECK( qf_count_key_value(&qf,150,0)==50);

    qf_insert(&qf,1500,0,50,false,false);
    qf_set_fixed_counter(&qf,1500,maximum_count);
    CHECK( qf_get_fixed_counter(&qf,1500)==maximum_count);

    qf_insert(&qf,3000,0,1,false,false);
    qf_set_fixed_counter(&qf,3000,maximum_count);
    CHECK( qf_get_fixed_counter(&qf,3000)==maximum_count);

    qf_insert(&qf,1500000,0,1,false,false);
    qf_set_fixed_counter(&qf,1500000,maximum_count);
    CHECK( qf_get_fixed_counter(&qf,1500000)==maximum_count);



    qf_destroy(&qf,true);
  }
}

TEST_CASE( "simple counting test" ,"[!hide]" ) {
  //except first item is inserted 5 times to full test _insert1
  QF qf;
  int counter_size=2;
  uint64_t qbits=8;
  uint64_t num_hash_bits=qbits+8;
  uint64_t maximum_count=(1ULL<<counter_size)-1;
  uint64_t count,fixed_counter;
  INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
  qf_init(&qf, (1ULL<<qbits), num_hash_bits, 0,counter_size, true, "", 2038074761);

  for(int i=0;i<=10;i++){
    qf_insert(&qf,100,0,1,false,false);
    count = qf_count_key_value(&qf, 100, 0);
    fixed_counter=qf_get_fixed_counter(&qf,100);
    INFO("Counter = "<<count<<" fixed counter = "<<fixed_counter)
    CHECK(count == (1+i));
  }

  qf_insert(&qf,1500,0,50,false,false);
  count = qf_count_key_value(&qf, 1500, 0);
  fixed_counter=qf_get_fixed_counter(&qf,1500);
  INFO("Counter = "<<count<<" fixed counter = "<<fixed_counter)
  CHECK(count == (50));

  qf_insert(&qf,1600,0,50,false,false);
  count = qf_count_key_value(&qf, 1500, 0);
  fixed_counter=qf_get_fixed_counter(&qf,1500);
  INFO("Counter = "<<count<<" fixed counter = "<<fixed_counter)
  CHECK(count == (50));

}
TEST_CASE( "Inserting items( repeated 1 time) in cqf(90% load factor )" ,"[!hide]") {
  //except first item is inserted 5 times to full test _insert1
  QF qf;
  int counter_size=2;
  uint64_t qbits=16;
  uint64_t num_hash_bits=qbits+8;
  uint64_t maximum_count=(1ULL<<counter_size)-1;
  INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
  qf_init(&qf, (1ULL<<qbits), num_hash_bits, 0,counter_size, true, "", 2038074761);

  uint64_t nvals = (1ULL<<qbits)*2;
  uint64_t *vals;
  vals = (uint64_t*)malloc(nvals*sizeof(vals[0]));
  for(int i=0;i<nvals;i++)
  {
    vals[i]=rand();
    vals[i]=(vals[i]<<32)|rand();
    vals[i]=vals[i]%(qf.metadata->range);
  }

  double loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
  uint64_t insertedItems=0;

    qf_insert(&qf,vals[0],0,1,false,false);
    qf_insert(&qf,vals[0],0,1,false,false);
    qf_insert(&qf,vals[0],0,1,false,false);
    qf_insert(&qf,vals[0],0,1,false,false);
  // for(int i=0;i<32;i++)
  // {
  //   cout<<get_fixed_counter(&qf,i)<<"-";
  // }
  //cout<<endl;
  while(loadFactor<0.9){

    qf_insert(&qf,vals[insertedItems],0,1,false,false);
    // for(int i=0;i<32;i++)
    // {
    //   cout<<get_fixed_counter(&qf,i)<<"-";
    // }
    // cout<<endl;
    insertedItems++;
    loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
  }
  INFO("Inserted Items = "<<insertedItems);

  uint64_t count;
  INFO("Fixed counter = "<<qf_get_fixed_counter(&qf,vals[0]));
  count = qf_count_key_value(&qf, vals[0], 0);
  CHECK(count >= 5);

  for(int i=1;i<insertedItems;i++)
  {
    count = qf_count_key_value(&qf, vals[i], 0);
    CHECK(count >= 1);
  }
  QFi qfi;
  qf_iterator(&qf, &qfi, 0);
  do {
    uint64_t key, value, count;
    qfi_get(&qfi, &key, &value, &count);
    count=qf_count_key_value(&qf, key, 0);
    if(key==vals[0]){
      CHECK(count >= 5);
    }
    else{
      CHECK(count >= 1);
    }

  } while(!qfi_next(&qfi));

  qf_destroy(&qf,true);

}


TEST_CASE( "Inserting items( repeated 50 times) in cqf(90% load factor )" ,"[!hide]") {
  QF qf;
  int counter_size=2;
  uint64_t qbits=16;
  uint64_t num_hash_bits=qbits+8;
  uint64_t maximum_count=(1ULL<<counter_size)-1;
  INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
  qf_init(&qf, (1ULL<<qbits), num_hash_bits, 0,counter_size, true, "", 2038074761);

  uint64_t nvals = (1ULL<<qbits);
  uint64_t *vals;
  vals = (uint64_t*)malloc(nvals*sizeof(vals[0]));
  for(int i=0;i<nvals;i++)
  {
    vals[i]=rand();
    vals[i]=(vals[i]<<32)|rand();
    vals[i]=vals[i]%(qf.metadata->range);
  }
  double loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
  uint64_t insertedItems=0;
  while(loadFactor<0.9){

    qf_insert(&qf,vals[insertedItems],0,50,false,false);
    insertedItems++;
    loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
  }
  uint64_t count;
  for(int i=0;i<insertedItems;i++)
  {
    count = qf_count_key_value(&qf, vals[i], 0);
    CHECK(count >= 50);
  }
  QFi qfi;
  qf_iterator(&qf, &qfi, 0);
  do {
    uint64_t key, value, count;
    qfi_get(&qfi, &key, &value, &count);
    count=qf_count_key_value(&qf, key, 0);
    CHECK(count >= 50);
  } while(!qfi_next(&qfi));

  qf_destroy(&qf,true);

}


TEST_CASE( "Inserting items( repeated 1-1000 times) in cqf(90% load factor )","[!mayfail]" ) {
  QF qf;
  int counter_size=2;
  uint64_t qbits=11;
  uint64_t num_hash_bits=qbits+8;
  uint64_t maximum_count=(1ULL<<counter_size)-1;
  INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
  qf_init(&qf, (1ULL<<qbits), num_hash_bits, 0,counter_size, true, "", 2038074761);

  uint64_t nvals = (1ULL<<qbits);
  //uint64_t nvals = 3;
  uint64_t *vals;
  uint64_t *nRepetitions;
  vals = (uint64_t*)malloc(nvals*sizeof(vals[0]));
  nRepetitions= (uint64_t*)malloc(nvals*sizeof(nRepetitions[0]));
  uint64_t count;

  for(int i=0;i<nvals;i++)
  {
    vals[i]=rand();
    vals[i]=(vals[i]<<32)|rand();
    vals[i]=vals[i]%(qf.metadata->range);

    nRepetitions[i]=rand()%258;
  }
  double loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
  uint64_t insertedItems=0;
  while(insertedItems<nvals && loadFactor<0.9){
    INFO("Inserting "<< vals[insertedItems] << " Repeated "<<nRepetitions[insertedItems]);
    qf_insert(&qf,vals[insertedItems],0,nRepetitions[insertedItems],false,false);
    qf_dump(&qf);
    INFO("Load factor = "<<loadFactor <<" inserted items = "<<insertedItems);
    count = qf_count_key_value(&qf, vals[insertedItems], 0);
    //CHECK(count >= nRepetitions[insertedItems]);
    insertedItems++;
    loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;

  }
  INFO("Load factor = "<<loadFactor <<" inserted items = "<<insertedItems);

  for(int i=0;i<insertedItems;i++)
  {
    count = qf_count_key_value(&qf, vals[i], 0);
    INFO("value = "<<vals[i]<<" Repeated " <<nRepetitions[i]);
    CHECK(count >= nRepetitions[i]);
  }

  qf_destroy(&qf,true);

}

TEST_CASE( "Writing and Reading to/from Disk","[!hide]") {
  QF qf;
  int counter_size=2;
  uint64_t qbits=16;
  uint64_t num_hash_bits=qbits+8;
  uint64_t maximum_count=(1ULL<<counter_size)-1;
  INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
  qf_init(&qf, (1ULL<<qbits), num_hash_bits, 0,counter_size, true, "", 2038074761);

  uint64_t nvals = (1ULL<<qbits);
  //uint64_t nvals = 3;
  uint64_t *vals;
  uint64_t *nRepetitions;
  vals = (uint64_t*)malloc(nvals*sizeof(vals[0]));
  nRepetitions= (uint64_t*)malloc(nvals*sizeof(nRepetitions[0]));
  uint64_t count;

  for(int i=0;i<nvals;i++)
  {
    vals[i]=rand();
    vals[i]=(vals[i]<<32)|rand();
    vals[i]=vals[i]%(qf.metadata->range);

    nRepetitions[i]=rand()%258;
  }
  double loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
  uint64_t insertedItems=0;
  while(insertedItems<nvals && loadFactor<0.9){
    qf_insert(&qf,vals[insertedItems],0,nRepetitions[insertedItems],false,false);
    count = qf_count_key_value(&qf, vals[insertedItems], 0);
    CHECK(count >= nRepetitions[insertedItems]);
    insertedItems++;
    loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
  }
  INFO("Load factor = "<<loadFactor <<" inserted items = "<<insertedItems);
  INFO("nslots ="<<qf.metadata->nslots);
  qf_serialize(&qf,"tmp.ser");
  qf_destroy(&qf,true);

  SECTION("Reading using qf_read(mmap)"){
    QF qf2;
    qf_read(&qf2,"tmp.ser");
    INFO("nslots ="<<qf2.metadata->nslots);
    for(int i=0;i<insertedItems;i++)
    {
      count = qf_count_key_value(&qf2, vals[i], 0);
      INFO("value = "<<vals[i]<<" Repeated " <<nRepetitions[i]);
      CHECK(count >= nRepetitions[i]);
    }

    qf_destroy(&qf2,false);
  }

  SECTION("Reading using deserialize "){
    qf_deserialize(&qf,"tmp.ser");

    for(int i=0;i<insertedItems;i++)
    {
      count = qf_count_key_value(&qf, vals[i], 0);
      INFO("value = "<<vals[i]<<" Repeated " <<nRepetitions[i]);
      CHECK(count >= nRepetitions[i]);
    }

    qf_destroy(&qf,true);
  }



}

TEST_CASE( "MMap test","[!hide]") {
  QF qf;
  int counter_size=2;
  uint64_t qbits=16;
  uint64_t num_hash_bits=qbits+8;
  uint64_t maximum_count=(1ULL<<counter_size)-1;
  INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
  qf_init(&qf, (1ULL<<qbits), num_hash_bits, 0,counter_size, false, "tmp.ser", 2038074761);

  uint64_t nvals = (1ULL<<qbits);
  //uint64_t nvals = 3;
  uint64_t *vals;
  uint64_t *nRepetitions;
  vals = (uint64_t*)malloc(nvals*sizeof(vals[0]));
  nRepetitions= (uint64_t*)malloc(nvals*sizeof(nRepetitions[0]));
  uint64_t count;

  for(int i=0;i<nvals;i++)
  {
    vals[i]=rand();
    vals[i]=(vals[i]<<32)|rand();
    vals[i]=vals[i]%(qf.metadata->range);

    nRepetitions[i]=rand()%258;
  }
  double loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
  uint64_t insertedItems=0;
  while(insertedItems<nvals && loadFactor<0.9){
    qf_insert(&qf,vals[insertedItems],0,nRepetitions[insertedItems],false,false);
    count = qf_count_key_value(&qf, vals[insertedItems], 0);
    CHECK(count >= nRepetitions[insertedItems]);
    insertedItems++;
    loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
  }
  INFO("Load factor = "<<loadFactor <<" inserted items = "<<insertedItems<<" "<<vals[0]);

  for(int i=0;i<insertedItems;i++)
  {
    INFO("Check = "<<vals[i]);
    count = qf_count_key_value(&qf, vals[i], 0);
    INFO("value = "<<vals[i]<<" Repeated " <<nRepetitions[i]);
    CHECK(count >= nRepetitions[i]);
  }

  qf_destroy(&qf,false);

}


TEST_CASE( "Removing items from cqf(90% load factor )","[!hide]") {
  QF qf;
  int counter_size=2;
  uint64_t qbits=16;
  uint64_t num_hash_bits=qbits+8;
  uint64_t maximum_count=(1ULL<<counter_size)-1;
  qf_init(&qf, (1ULL<<qbits), num_hash_bits, 0,counter_size, true, "", 2038074761);

  uint64_t nvals = (1ULL<<qbits);
  uint64_t *vals;
  vals = (uint64_t*)malloc(nvals*sizeof(vals[0]));
  for(int i=0;i<nvals;i++)
  {
    uint64_t newvalue=0;
    while(newvalue==0){
      newvalue=rand();
      newvalue=(newvalue<<32)|rand();
      newvalue=newvalue%(qf.metadata->range);
      for(int j=0;j<i;j++)
      {
        if(vals[j]==newvalue)
        {
          newvalue=0;
          break;
        }
      }
    }
    vals[i]=newvalue;
  }
  double loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
  uint64_t insertedItems=0;
  while(loadFactor<0.9){

    qf_insert(&qf,vals[insertedItems],0,50,false,false);
    insertedItems++;
    loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
  }
  uint64_t count;
  for(int i=0;i<insertedItems;i++)
  {
    if(i%2==0){
      count = qf_count_key_value(&qf, vals[i], 0);
      if(count==100){
        printf("coubn ==100\n" );
      }
    _remove(&qf,vals[i],50);
    count = qf_count_key_value(&qf, vals[i], 0);
    CHECK(count ==0);
    }
  }
  for(int i=0;i<insertedItems;i++)
  {
    count = qf_count_key_value(&qf, vals[i], 0);
    if(i%2==1){
    CHECK(count >= 50);
    }
    else{
      if(count!=0){
        INFO("ERROR "<<vals[i]<<" Not deleted index= "<<i)
        //printf("%lu not delete at index %lu\n", vals[i],i);
      }
      CHECK(count ==0);

    }
  }
  QFi qfi;
  qf_iterator(&qf, &qfi, 0);
  do {
    uint64_t key, value, count;
    qfi_get(&qfi, &key, &value, &count);
    count=qf_count_key_value(&qf, key, 0);
    CHECK(count >= 50);
  } while(!qfi_next(&qfi));

  qf_destroy(&qf,true);

}
TEST_CASE( "Merging Cqf","[!hide]") {
  QF cf,cf1,cf2;
 QFi cfi;
 uint64_t qbits = 18;
 uint64_t small_qbits=qbits;
 uint64_t nhashbits = qbits + 8;
 uint64_t small_nhashbits=small_qbits+8;
 uint64_t nslots = (1ULL << qbits);
 uint64_t small_nslots=(1ULL << small_qbits);
 uint64_t nvals = 250*nslots/1000;
 uint64_t *vals;
 uint64_t counter_size=3;
 /* Initialise the CQF */

 INFO("Initialize first cqf size ="<<nslots<<", hashbits="<<nhashbits);
 qf_init(&cf, nslots, nhashbits, 0,counter_size, true, "", 2038074761);
 INFO("Initialize second cqf size ="<<small_nslots<<", hashbits="<<small_nhashbits);

 qf_init(&cf1, small_nslots, small_nhashbits, 0,counter_size, true, "", 2038074761);
 INFO("Initialize third cqf size ="<<small_nslots<<", hashbits="<<small_nhashbits);
 qf_init(&cf2, small_nslots, small_nhashbits, 0,counter_size, true, "", 2038074761);
 /* Generate random values */
 vals = (uint64_t*)malloc(nvals*sizeof(vals[0]));

 for(int i=0;i<nvals;i++)
 {
   vals[i]=rand();
   vals[i]=(vals[i]<<32)|rand();
 }


 /* Insert vals in the CQF */
 for (uint64_t i = 0; i < (nvals*2)/3; i++) {
   vals[i]=vals[i]%cf1.metadata->range;
   if(i%2==1){
     qf_insert(&cf2, vals[i], 0, 50,false,false);
   }
   else{
     qf_insert(&cf1, vals[i], 0, 50,false,false);
   }

 }
 qf_merge(&cf1,&cf2,&cf);

 for (uint64_t i = (nvals*2)/3; i <nvals; i++) {
   vals[i]=vals[i]%cf.metadata->range;
   qf_insert(&cf, vals[i], 0, 50,false,false);
   }

 for (uint64_t i = 0; i < nvals; i++) {

   uint64_t count = qf_count_key_value(&cf, vals[i]%cf.metadata->range, 0);
   CHECK(count>=50);
 }

 /* Initialize an iterator */
 qf_iterator(&cf, &cfi, 0);
 do {
   uint64_t key, value, count;
   qfi_get(&cfi, &key, &value, &count);
   CHECK(count>=50);
 } while(!qfi_next(&cfi));

}



TEST_CASE( "Inserting items( repeated 50 times)  and set fixed size counters in cqf(90% load factor )","[!hide]") {
  QF qf;
  int counter_size=3;
  uint64_t qbits=7;
  uint64_t num_hash_bits=qbits+8;
  uint64_t maximum_count=(1ULL<<counter_size)-1;
  INFO("Counter size = "<<counter_size<<" max count= "<<maximum_count);
  qf_init(&qf, (1ULL<<qbits), num_hash_bits, 0,counter_size, true, "", 2038074761);

  uint64_t nvals = (1ULL<<qbits);
  uint64_t *vals;
  vals = (uint64_t*)malloc(nvals*sizeof(vals[0]));
  for(int i=0;i<nvals;i++)
  {
    vals[i]=rand();
    vals[i]=(vals[i]<<32)|rand();
    vals[i]=vals[i]%(qf.metadata->range);
  }
  double loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
  uint64_t insertedItems=0;
  uint64_t count;
  while(loadFactor<0.9){

    qf_insert(&qf,vals[insertedItems],0,50,false,false);
    qf_set_fixed_counter(&qf,vals[insertedItems],vals[insertedItems]%(maximum_count+1));

    count = qf_get_fixed_counter(&qf,vals[insertedItems]);
    CHECK(count == vals[insertedItems]%(maximum_count+1));
    insertedItems++;
    loadFactor=(double)qf.metadata->noccupied_slots/(double)qf.metadata->nslots;
    for(int i=0;i<insertedItems;i++)
    {
      count = qf_count_key_value(&qf, vals[i], 0);
      REQUIRE(count >= 50);
      count = qf_get_fixed_counter(&qf,vals[i]);
      INFO("bug in  "<<i<<" of "<<insertedItems<<" loadFactor "<<loadFactor);
      REQUIRE(count == vals[i]%(maximum_count+1));
    }

  }

  for(int i=0;i<insertedItems;i++)
  {
    count = qf_count_key_value(&qf, vals[i], 0);
    CHECK(count >= 50);
    count = qf_get_fixed_counter(&qf,vals[i]);
    CHECK(count == vals[i]%(maximum_count+1));
  }
  QFi qfi;
  qf_iterator(&qf, &qfi, 0);
  do {
    uint64_t key, value, count;
    qfi_get(&qfi, &key, &value, &count);
    count=qf_count_key_value(&qf, key, 0);
    CHECK(count >= 50);
    count = qf_get_fixed_counter(&qf,key);
    CHECK(count == key%(maximum_count+1));
  } while(!qfi_next(&qfi));

  qf_destroy(&qf,true);

}
