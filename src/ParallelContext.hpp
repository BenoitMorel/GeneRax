#ifndef _JOINTSEARCH_PARALLEL_CONTEXT_HPP_
#define _JOINTSEARCH_PARALLEL_CONTEXT_HPP_

#include <fstream>
#include <vector>
#include <mpi.h>

using namespace std;

class ParallelContext {
public:
  static void init(void *commPtr);
  static void finalize();
  static int getRank();
  static int getSize();
  static void setComm(MPI_Comm newComm);
  static void allGatherDouble(double localValue, vector<double> &allValues);
  static void broadcoastInt(int fromRank, int &value);
  static void broadcoastDouble(int fromRank, double &value);
  static int getBestLL(double &bestLL, int &bestRank);
  

  static int getBegin(int elems);
  static int getEnd(int elems);

private:
  static ofstream sink;
  static MPI_Comm comm;
  static bool ownMPIContext;
  
};


#endif
