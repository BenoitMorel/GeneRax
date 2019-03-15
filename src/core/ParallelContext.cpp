#include "ParallelContext.hpp"
#include <algorithm>

ofstream ParallelContext::sink("/dev/null");
bool ParallelContext::ownMPIContext(true);
stack<MPI_Comm> ParallelContext::_commStack;
stack<bool> ParallelContext::_ownsMPIContextStack;


void ParallelContext::init(void *commPtr)
{
  if (commPtr) {
    setComm(*((MPI_Comm*)commPtr));
    setOwnMPIContext(false);
  } else {
    MPI_Init(0, 0);
    setComm(MPI_COMM_WORLD);
    setOwnMPIContext(true);
  }

}


void ParallelContext::finalize()
{
  if (_ownsMPIContextStack.top()) {
    MPI_Finalize();
  }
  _commStack.pop();
  _ownsMPIContextStack.pop();
}

int ParallelContext::getRank() 
{
  int rank = 0;
  MPI_Comm_rank(getComm(), &rank);
  return rank;
}

int ParallelContext::getSize() 
{
  int size = 0;
  MPI_Comm_size(getComm(), &size);
  return size;
}
  
void ParallelContext::setOwnMPIContext(bool own)
{
  _ownsMPIContextStack.push(own);
  ownMPIContext = own;
  
}

void ParallelContext::setComm(MPI_Comm newComm)
{
  _commStack.push(newComm);
}



int ParallelContext::getBegin(int elems)
{
  return (getRank() * elems) / getSize();
}
 
int ParallelContext::getEnd(int elems)
{
  return min(elems, ((getRank() + 1) * elems) / getSize());
}
  
void ParallelContext::sumDouble(double &value)
{
  double sum = 0;
  barrier();
  MPI_Allreduce(&value, &sum, 1, MPI_DOUBLE, MPI_SUM, getComm());
  barrier();
  value = sum;
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
    getComm());
}

void ParallelContext::concatenateIntVectors(const vector<int> &localVector, vector<int> &globalVector)
{
  globalVector.resize(getSize() * localVector.size(), 0);
  MPI_Allgather(
    &localVector[0],
    localVector.size(),
    MPI_INT,
    &(globalVector[0]),
    localVector.size(),
    MPI_INT,
    getComm());
}
  
void ParallelContext::broadcastInt(int fromRank, int &value)
{
  MPI_Bcast(
    &value,
    1,
    MPI_INT,
    fromRank,
    getComm());
}

void ParallelContext::broadcastDouble(int fromRank, double &value)
{
  MPI_Bcast(
    &value,
    1,
    MPI_DOUBLE,
    fromRank,
    getComm());
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

void ParallelContext::barrier()
{
  MPI_Barrier(getComm());
}

void ParallelContext::abort(int errorCode)
{ 
  if (_ownsMPIContextStack.top()) {
    MPI_Finalize();
    exit(errorCode);
  } else {
    throw ParallelException(errorCode); 
  }
}


