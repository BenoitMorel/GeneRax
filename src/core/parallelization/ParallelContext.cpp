#include "parallelization/ParallelContext.hpp"
#include <algorithm>
#include <cassert>

std::ofstream ParallelContext::sink("/dev/null");
bool ParallelContext::ownMPIContext(true);
std::stack<MPI_Comm> ParallelContext::_commStack;
std::stack<bool> ParallelContext::_ownsMPIContextStack;
bool ParallelContext::_mpiEnabled = false;

void ParallelContext::init(void *commPtr)
{
  if (commPtr && *static_cast<int *>(commPtr) == -1) {
    _mpiEnabled = false;
    return;
  }
#ifdef WITH_MPI
  if (commPtr) {
    _mpiEnabled = true;
    setComm(*(static_cast<MPI_Comm*>(commPtr)));
    setOwnMPIContext(false);
  } else {
    _mpiEnabled = true;
    MPI_Init(0, 0);
    setComm(MPI_COMM_WORLD);
    setOwnMPIContext(true);
    if (getSize() == 1) {
      finalize();
      _mpiEnabled = false;
      setOwnMPIContext(false);
    }   
  }
#else
  assert(false);
#endif  


}
  
ParallelException::~ParallelException() 
{
}


void ParallelContext::finalize()
{
  if (!_mpiEnabled) {
    return;
  }
#ifdef WITH_MPI
  if (_ownsMPIContextStack.top()) {
    MPI_Finalize();
  }
  _commStack.pop();
  _ownsMPIContextStack.pop();
#else
  assert(false);
#endif
}

unsigned int ParallelContext::getRank() 
{
#ifdef WITH_MPI
  if (!_mpiEnabled) {
    return 0;
  }
  int rank = 0;
  MPI_Comm_rank(getComm(), &rank);
  return static_cast<unsigned int>(rank);
#else
  return 0;
#endif
}

unsigned int ParallelContext::getSize() 
{
#ifdef WITH_MPI
  if (!_mpiEnabled) {
    return 1;
  }
  int size = 0;
  MPI_Comm_size(getComm(), &size);
  return static_cast<unsigned int>(size);
#else
  return 1;
#endif
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
#ifdef WITH_MPI
  if (!_mpiEnabled) {
    return;
  }
  double sum = 0;
  barrier();
  MPI_Allreduce(&value, &sum, 1, MPI_DOUBLE, MPI_SUM, getComm());
  value = sum;
#endif
}

void ParallelContext::sumUInt(unsigned int &value)
{
#ifdef WITH_MPI
  if (!_mpiEnabled) {
    return;
  }
  unsigned int sum = 0;
  barrier();
  MPI_Allreduce(&value, &sum, 1, MPI_UNSIGNED, MPI_SUM, getComm());
  value = sum;
#endif
}

void ParallelContext::sumVectorDouble(std::vector<double> &value)
{
#ifdef WITH_MPI
  if (!_mpiEnabled) {
    return;
  }
  std::vector<double> sum(value.size());
  barrier();
  MPI_Allreduce(&(value[0]), &(sum[0]), static_cast<int>(value.size()), MPI_DOUBLE, MPI_SUM, getComm());
  value = sum;
#endif
}

void ParallelContext::parallelAnd(bool &value)
{
#ifdef WITH_MPI
  if (!_mpiEnabled) {
    return;
  }
  int v = value ? 1 : 0;
  int res = 0;
  MPI_Allreduce(&v, &res, 1, MPI_INT, MPI_MIN, getComm());
  value = (res == 1);
#endif
}


void ParallelContext::allGatherDouble(double localValue, std::vector<double> &allValues) {
  if (!_mpiEnabled) {
    allValues.clear();
    allValues.push_back(localValue);
    return;
  }
#ifdef WITH_MPI
  allValues.resize(getSize());
  MPI_Allgather(
    &localValue,
    1,
    MPI_DOUBLE,
    &(allValues[0]),
    1,
    MPI_DOUBLE,
    getComm());
#else
  assert(false);
#endif
}

void ParallelContext::allGatherInt(int localValue, std::vector<int> &allValues)
{
  if (!_mpiEnabled) {
    allValues.clear();
    allValues.push_back(localValue);
    return;
  }
#ifdef WITH_MPI
  allValues.resize(getSize());
  MPI_Allgather(
    &localValue,
    1,
    MPI_INT,
    &(allValues[0]),
    1,
    MPI_INT,
    getComm());
#else
  assert(false);
#endif

}

void ParallelContext::concatenateIntVectors(const std::vector<int> &localVector, std::vector<int> &globalVector)
{
  if (!_mpiEnabled) {
    globalVector = localVector;
    return;
  }

#ifdef WITH_MPI
  globalVector.resize(getSize() * localVector.size(), 0);
  MPI_Allgather(
    &localVector[0],
    static_cast<int>(localVector.size()),
    MPI_INT,
    &(globalVector[0]),
    static_cast<int>(localVector.size()),
    MPI_INT,
    getComm());
#else
  assert(false);
#endif
}
  
void ParallelContext::concatenateUIntVectors(const std::vector<unsigned int> &localVector, 
  std::vector<unsigned int> &globalVector)
{
  if (!_mpiEnabled) {
    globalVector = localVector;
    return;
  }

#ifdef WITH_MPI
  globalVector.resize(getSize() * localVector.size(), 0);
  MPI_Allgather(
    &localVector[0],
    static_cast<int>(localVector.size()),
    MPI_UNSIGNED,
    &(globalVector[0]),
    static_cast<int>(localVector.size()),
    MPI_UNSIGNED,
    getComm());
#else
  assert(false);
#endif
}
  
void ParallelContext::broadcastInt(unsigned int fromRank, int &value)
{
  if (!_mpiEnabled) {
    return;
  }
#ifdef WITH_MPI
  MPI_Bcast(
    &value,
    1,
    MPI_INT,
    static_cast<int>(fromRank),
    getComm());
#endif
}

void ParallelContext::broadcastUInt(unsigned int fromRank, unsigned int &value)
{
  if (!_mpiEnabled) {
    return;
  }
#ifdef WITH_MPI
  MPI_Bcast(
    &value,
    1,
    MPI_UNSIGNED,
    static_cast<int>(fromRank),
    getComm());
#endif
}

void ParallelContext::broadcastDouble(unsigned int fromRank, double &value)
{
  if (!_mpiEnabled) {
    return;
  }
#ifdef WITH_MPI
  MPI_Bcast(
    &value,
    1,
    MPI_DOUBLE,
    static_cast<int>(fromRank),
    getComm());
#endif
}

void ParallelContext::maxUInt(unsigned int &value)
{
#ifdef WITH_MPI
  if (!_mpiEnabled) {
    return;
  }
  unsigned int sum = 0;
  barrier();
  MPI_Allreduce(&value, &sum, 1, MPI_UNSIGNED, MPI_MAX, getComm());
  value = sum;
#endif
}

unsigned int ParallelContext::getMax(double &value, unsigned int &bestRank)
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
#ifdef WITH_MPI
  MPI_Barrier(getComm());
#endif
}

void ParallelContext::abort(int errorCode)
{ 
  if (!_mpiEnabled) {
    exit(errorCode);
  }
#ifdef WITH_MPI
  if (_ownsMPIContextStack.top()) {
    MPI_Finalize();
    exit(errorCode);
  } else {
    throw ParallelException(errorCode); 
  }
#else
  assert(false);
#endif
}

bool ParallelContext::allowSchedulerSplitImplementation()
{
  return getSize() > 4;
}

bool ParallelContext::isRandConsistent()
{
#ifdef WITH_MPI
  std::vector<int> rands(getSize());
  int localValue = rand();
  allGatherInt(localValue, rands);
  for (auto value: rands) {
    if (value != rands[0]) {
      return false;
    }
  }
#endif
  return true;
}
