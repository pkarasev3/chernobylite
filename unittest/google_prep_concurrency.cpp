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


