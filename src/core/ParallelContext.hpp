#pragma once

#include <fstream>
#include <vector>
#include <mpi.h>
#include <stack>

using namespace std;
class ParallelException: public exception
{
public:
  ParallelException(int errorCode):errorCode_(errorCode) 
  {
    msg_ = "Program failed with error " + to_string(errorCode);
  }

  virtual const char* what() const throw()
  {
    return msg_.c_str();
  }
private:
  int errorCode_;
  string msg_;
};

class ParallelContext {
public:
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
  static int getRank();

  /**
   *  @return the number of MPI ranks
   */
  static int getSize();

  /**
   * Gather the values of each rank into a global vector
   *  @param localValue: input value for this rank
   *  @param allValues: output values for all rank (indexed with the rank index)
   */
  static void allGatherDouble(double localValue, vector<double> &allValues);

  static void concatenateIntVectors(const vector<int> &localVector, vector<int> &globalVector);

  static void sumDouble(double &value);

  /**
   *  broadcast a value from a given rank
   *  @param fromRank: rank from which we want the value
   *  @param value: input value for this rank, output value for the other ranks
   */
  static void broadcastInt(int fromRank, int &value);
  static void broadcastDouble(int fromRank, double &value);

  /**
   *  Get the highest value from all ranks
   *  @param value: as input, the value for the current rank. As output, 
   *    the highest value among the ranks
   *  @param bestRank: the rank that has the highest value
   */
  static int getMax(double &value, int &bestRank);
  

  /**
   *  When having a given number of independant tasks, subdivide
   *  them into chunks, and assign each chunk to one rank
   *  @param elems: the number of elements
   *  @return the begin/end of the interval of the chunk for this rank
   */
  static int getBegin(int elems);
  static int getEnd(int elems);

  static void barrier();
  static void abort(int errorCode);

  static MPI_Comm getComm() {return _commStack.top();}
private:
  static void setComm(MPI_Comm newComm);
  static void setOwnMPIContext(bool own);
  static ofstream sink;
  static bool ownMPIContext;
  static stack<MPI_Comm> _commStack;
  static stack<bool> _ownsMPIContextStack;
  static bool _mpiEnabled;
};

