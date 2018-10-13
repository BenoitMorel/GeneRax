#include "ParallelContext.hpp"
#include <algorithm>

ofstream ParallelContext::sink("/dev/null");
MPI_Comm ParallelContext::comm(MPI_COMM_WORLD);
bool ParallelContext::ownMPIContext(true);

void ParallelContext::init(void *commPtr)
{
  if (commPtr) {
    comm = *((MPI_Comm*)commPtr);
    ownMPIContext = false; 
  }

  if (ownMPIContext) {
    MPI_Init(0, 0);
  }
}


void ParallelContext::finalize()
{
  if (ownMPIContext) {
    MPI_Finalize();
  }
}

int ParallelContext::getRank() 
{
  int rank = 0;
  MPI_Comm_rank(comm, &rank);
  return rank;
}

int ParallelContext::getSize() 
{
  int size = 0;
  MPI_Comm_size(comm, &size);
  return size;
}

void ParallelContext::setComm(MPI_Comm newComm)
{
  comm = newComm;
}



int ParallelContext::getBegin(int elems)
{
  int perRankElements = elems / getSize();
  return getRank() * perRankElements;
}
 
int ParallelContext::getEnd(int elems)
{
  int perRankElements = elems / getSize();
  return min(elems, (getRank() + 1) * perRankElements);

}

void ParallelContext::allGatherDouble(double localValue, vector<double> &allValues) {
  allValues.resize(getSize());
  MPI_Allgather(
    &localValue,
    1,
    MPI_DOUBLE,
    &(allValues[0]),
    1,
    MPI_DOUBLE,
    comm);
}
  
void ParallelContext::broadcoastInt(int fromRank, int &value)
{
  MPI_Bcast(
    &value,
    1,
    MPI_INT,
    fromRank,
    comm);
}

void ParallelContext::broadcoastDouble(int fromRank, double &value)
{
  MPI_Bcast(
    &value,
    1,
    MPI_DOUBLE,
    fromRank,
    comm);
}

int ParallelContext::getBestLL(double &bestLL, int &bestRank)
{
  vector<double> allValues;
  allGatherDouble(bestLL, allValues);
  bestRank = 0;
  for (int i = 0; i < allValues.size(); ++i) {
    if (allValues[i] > allValues[bestRank]) {
      bestRank = i;
    }
  }
  bestLL = allValues[bestRank];
  return bestRank;
}

