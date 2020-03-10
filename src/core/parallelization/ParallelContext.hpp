#pragma once

#include <fstream>
#include <vector>
#include <stack>
#ifdef WITH_MPI
  #include <mpi.h>
#else
  typedef int MPI_Comm;
#endif


/**
 * Singleton class that handles parallelization routines
 */
class ParallelContext {
public:
  ParallelContext() = delete;

  /**
   *  Initialize the parallel context. Must always be called at the start of the programm
   *  @param commPtr a pointer to the MPI communicator used by this programm
   */
  static void init(void *commPtr);

  /**
   *  Terminates the parallel context. Must always be called at the end of the programm
   */
  static void finalize();

  /**
   *  @return the MPI rank
   */
  static unsigned int getRank();

  /**
   *  @return the number of MPI ranks
   */
  static unsigned int getSize();

  /**
   * Gather the values of each rank into a global std::vector
   *  @param localValue input value for this rank
   *  @param allValues output values for all rank (indexed with the rank index)
   */
  static void allGatherDouble(double localValue, std::vector<double> &allValues);
  static void allGatherInt(int localValue, std::vector<int> &allValues);

  static bool isRandConsistent();
  static bool isIntEqual(int value);
  static bool isDoubleEqual(double value);

  static void concatenateIntVectors(const std::vector<int> &localVector, std::vector<int> &globalVector);
  static void concatenateUIntVectors(const std::vector<unsigned int> &localVector, 
    std::vector<unsigned int> &globalVector);

  static void sumDouble(double &value);
  static void sumUInt(unsigned int &value);
  static void sumVectorDouble(std::vector<double> &value);
  static void maxUInt(unsigned int &value);
  
  static void parallelAnd(bool &value);

  /**
   *  broadcast a value from a given rank
   *  @param fromRank rank from which we want the value
   *  @param value input value for this rank, output value for the other ranks
   */
  static void broadcastInt(unsigned int fromRank, int &value);
  static void broadcastUInt(unsigned int fromRank, unsigned int &value);
  static void broadcastDouble(unsigned int fromRank, double &value);

  /**
   *  Get the highest value from all ranks
   *  @param value as input, the value for the current rank. As output, 
   *    the highest value among the ranks
   *  @param bestRank the rank that has the highest value
   */
  static unsigned int getMax(double &value, unsigned int &bestRank);
  
  

  /**
   *  When having a given number of independant tasks, subdivide
   *  them into chunks, and assign each chunk to one rank
   *  @param elems the number of elements
   *  @return the begin/end of the interval of the chunk for this rank
   */
  static unsigned int getBegin(unsigned int elems);
  static unsigned int getEnd(unsigned int elems);

  static void barrier();
  static void abort(int errorCode);

  /**
   *  Do we have enough cores to allow the split implementation
   *  in the scheduler?
   */
  static bool allowSchedulerSplitImplementation();
  static MPI_Comm &getComm() {return _commStack.top();}

  static void pushSequentialContext();
  static void popContext();

private:
  static void setComm(MPI_Comm newComm);
  static void setOwnMPIContext(bool own);
  static std::ofstream sink;
  static bool ownMPIContext;
  static std::stack<MPI_Comm> _commStack;
  static std::stack<bool> _ownsMPIContextStack;
  static bool _mpiEnabled;
  class ParallelException: public std::exception
  {
  public:
    ParallelException(int errorCode) 
    {
      msg_ = "Program failed with error " + std::to_string(errorCode);
    }
    
    virtual ~ParallelException() {}

    virtual const char* what() const noexcept
    {
      return msg_.c_str();
    }
  private:
    std::string msg_;
  };
};

