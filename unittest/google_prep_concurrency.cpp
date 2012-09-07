/** "Concurrent read/write into data structure S"

  */
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>
#include <iostream>
#include <vector>
#include <cstdio>
using namespace std;

boost::mutex  GlobalMutex;

struct Worker {

  Worker():key(0),L(NULL)
  {
  }
  void operator()() { // do stuff
      boost::mutex::scoped_lock  lock(GlobalMutex, boost::try_to_lock);
      if( lock ) {
        int prev_size = (int) L->size();
        for( int crunch = 0; crunch < 1024; crunch++ ) {
          for( auto i = L->begin(); i != L->end(); ++i ) {
            *i = key; // writing to the shared memory
          }

          L->push_back( L->back() + 1 ); // writing to the shared memory
          if ( !L->empty() && (rand()%2 == 0) ) {
            L->pop_front();
          }
        }
        fprintf( stdout, "key=%d, |L|=%d, prev|L|=%d \n",
                             key,(int)L->size(),prev_size);
      } else {
        fprintf( stdout, "key=%d, Failed to get lock\n",key );
      }
  }
  int key;          // number of calls to this worker
  std::list<int>*  L; // a pointer to shared memory object
};

/**
 * 
 * Sample Output, |L| is the list size, key is the # of thread

key=1, |L|=506, prev|L|=1
key=3, Failed to get lock
key=2, |L|=1023, prev|L|=506
key=5, Failed to get lock
key=4, |L|=1491, prev|L|=1023
key=7, Failed to get lock
key=8, Failed to get lock
key=6, |L|=2015, prev|L|=1491
key=9, |L|=2550, prev|L|=2015
 * */

int main( int argc, char* argv[] )
{
  std::list<int> shared;
  shared.push_back(0);

  int k    = 0;
  int kmax = 10;
  std::vector<boost::shared_ptr<boost::thread> > threads;

  while( ++k < kmax )  {
    Worker w;
    w.L   = &shared;
    w.key = k;
    boost::shared_ptr<boost::thread>  T(new boost::thread(w) );
    threads.push_back( T );
    usleep(10000);
  }

  for(int i = 0; i < (int)threads.size(); i++ ) {
    threads[i]->join();
  }
  cout << endl;
  return 0;
}


