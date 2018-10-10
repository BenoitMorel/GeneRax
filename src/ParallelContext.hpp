#ifndef _JOINTSEARCH_PARALLEL_CONTEXT_HPP_
#define _JOINTSEARCH_PARALLEL_CONTEXT_HPP_

#include <fstream>
#include <mpi.h>

using namespace std;

class ParallelContext {
public:
  static void init();
  static void finalize();
  static int getRank();
  static int getSize();
  static void setComm(MPI_Comm newComm);
private:
  static ofstream sink;
  static MPI_Comm comm;

};


#endif
