#include "ParallelContext.hpp"


ofstream ParallelContext::sink("/dev/null");
MPI_Comm ParallelContext::comm(MPI_COMM_WORLD);

void ParallelContext::init()
{
  MPI_Init(0, 0);
  if (getRank() != 0) {
    std::cout.rdbuf(sink.rdbuf());
    std::cerr.rdbuf(sink.rdbuf());      
  }
}


void ParallelContext::finalize()
{
  MPI_Finalize();
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


