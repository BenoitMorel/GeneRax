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
  return (getRank() * elems) / getSize();
}
 
int ParallelContext::getEnd(int elems)
{
  return min(elems, ((getRank() + 1) * elems) / getSize());
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
  
void ParallelContext::broadcastInt(int fromRank, int &value)
{
  MPI_Bcast(
    &value,
    1,
    MPI_INT,
    fromRank,
    comm);
}

void ParallelContext::broadcastDouble(int fromRank, double &value)
{
  MPI_Bcast(
    &value,
    1,
    MPI_DOUBLE,
    fromRank,
    comm);
}

int ParallelContext::getMax(double &value, int &bestRank)
{
  vector<double> allValues;
  allGatherDouble(value, allValues);
  bestRank = 0;
  for (unsigned int i = 0; i < allValues.size(); ++i) {
    if (allValues[i] > allValues[bestRank]) {
      bestRank = i;
    }
  }
  value = allValues[bestRank];
  return bestRank;
}
  
void ParallelContext::abort(int errorCode)
{ 
  if (ownMPIContext) {
    MPI_Finalize();
    exit(errorCode);
  } else {
    throw ParallelException(errorCode); 
  }
}

