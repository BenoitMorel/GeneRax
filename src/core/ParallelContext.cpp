#include "ParallelContext.hpp"
#include <algorithm>

std::ofstream ParallelContext::sink("/dev/null");
bool ParallelContext::ownMPIContext(true);
std::stack<MPI_Comm> ParallelContext::_commStack;
std::stack<bool> ParallelContext::_ownsMPIContextStack;
bool ParallelContext::_mpiEnabled = true;

void ParallelContext::init(void *commPtr)
{
  if (commPtr) {
    if (*static_cast<int *>(commPtr) == -1) {
      _mpiEnabled = false;
      return;
    }
    _mpiEnabled = true;
    setComm(*(static_cast<MPI_Comm*>(commPtr)));
    setOwnMPIContext(false);
  } else {
    _mpiEnabled = true;
    MPI_Init(0, 0);
    setComm(MPI_COMM_WORLD);
    setOwnMPIContext(true);
  }

}


void ParallelContext::finalize()
{
  if (!_mpiEnabled) {
    return;
  }
  if (_ownsMPIContextStack.top()) {
    MPI_Finalize();
  }
  _commStack.pop();
  _ownsMPIContextStack.pop();
}

unsigned int ParallelContext::getRank() 
{
  if (!_mpiEnabled) {
    return 0;
  }
  int rank = 0;
  MPI_Comm_rank(getComm(), &rank);
  return static_cast<unsigned int>(rank);
}

unsigned int ParallelContext::getSize() 
{
  if (!_mpiEnabled) {
    return 1;
  }
  int size = 0;
  MPI_Comm_size(getComm(), &size);
  return static_cast<unsigned int>(size);
}
  
void ParallelContext::setOwnMPIContext(bool own)
{
  if (!_mpiEnabled) {
    return;
  }
  _ownsMPIContextStack.push(own);
  ownMPIContext = own;
  
}

void ParallelContext::setComm(MPI_Comm newComm)
{
  if (!_mpiEnabled) {
    return;
  }
  _commStack.push(newComm);
}



unsigned int ParallelContext::getBegin(unsigned int elems)
{
  return (getRank() * elems) / getSize();
}
 
unsigned int ParallelContext::getEnd(unsigned int elems)
{
  return std::min(elems, ((getRank() + 1) * elems) / getSize());
}
  
void ParallelContext::sumDouble(double &value)
{
  if (!_mpiEnabled) {
    return;
  }
  double sum = 0;
  barrier();
  barrier();
  MPI_Allreduce(&value, &sum, 1, MPI_DOUBLE, MPI_SUM, getComm());
  value = sum;
}


void ParallelContext::allGatherDouble(double localValue, std::vector<double> &allValues) {
  if (!_mpiEnabled) {
    allValues.clear();
    allValues.push_back(localValue);
    return;
  }
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

void ParallelContext::concatenateIntVectors(const std::vector<int> &localVector, std::vector<int> &globalVector)
{
  if (!_mpiEnabled) {
    globalVector = localVector;
    return;
  }

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
  
void ParallelContext::concatenateUIntVectors(const std::vector<unsigned int> &localVector, 
  std::vector<unsigned int> &globalVector)
{
  if (!_mpiEnabled) {
    globalVector = localVector;
    return;
  }

  globalVector.resize(getSize() * localVector.size(), 0);
  MPI_Allgather(
    &localVector[0],
    localVector.size(),
    MPI_UNSIGNED,
    &(globalVector[0]),
    localVector.size(),
    MPI_UNSIGNED,
    getComm());
}
  
void ParallelContext::broadcastInt(unsigned int fromRank, int &value)
{
  if (!_mpiEnabled) {
    return;
  }
  MPI_Bcast(
    &value,
    1,
    MPI_INT,
    fromRank,
    getComm());
}

void ParallelContext::broadcastUInt(unsigned int fromRank, unsigned int &value)
{
  if (!_mpiEnabled) {
    return;
  }
  MPI_Bcast(
    &value,
    1,
    MPI_UNSIGNED,
    fromRank,
    getComm());
}

void ParallelContext::broadcastDouble(unsigned int fromRank, double &value)
{
  if (!_mpiEnabled) {
    return;
  }
  MPI_Bcast(
    &value,
    1,
    MPI_DOUBLE,
    fromRank,
    getComm());
}

int ParallelContext::getMax(double &value, unsigned int &bestRank)
{
  if (!_mpiEnabled) {
    bestRank = 0;
    return bestRank;
  }
  std::vector<double> allValues;
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
  if (!_mpiEnabled) {
    return;
  }
  MPI_Barrier(getComm());
}

void ParallelContext::abort(int errorCode)
{ 
  if (!_mpiEnabled) {
    exit(errorCode);
  }
  if (_ownsMPIContextStack.top()) {
    MPI_Finalize();
    exit(errorCode);
  } else {
    throw ParallelException(errorCode); 
  }
}


