#include "IO/Logger.hpp"


Logger Logger::info;
Logger Logger::error;
Logger Logger::timed;
TimePoint Logger::start;
ofstream *Logger::logFile = 0;

Logger::Logger(): _os(&cout) {
  setType(lt_info);  
}


void Logger::init() {
  info.setType(lt_info);
  error.setStream(cout);
  error.setType(lt_error);
  error.setStream(cout);
  timed.setType(lt_timed);
  timed.setStream(cout);
  start = chrono::high_resolution_clock::now(); 
}
  
void Logger::initFileOutput(const string &output)
{
  if (ParallelContext::getRank()) {
    return;
  } 
  logFile = new ofstream(output + ".log");
}

