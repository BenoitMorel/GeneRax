#include "IO/Logger.hpp"


Logger Logger::info;
Logger Logger::error;
Logger Logger::timed;
TimePoint Logger::start;
std::ofstream *Logger::logFile = 0;
std::ofstream *Logger::saveLogFile = 0;

Logger::Logger(): _os(&std::cout) {
  setType(lt_info);  
}

void Logger::init() {
  info.setType(lt_info);
  error.setStream(std::cout);
  error.setType(lt_error);
  error.setStream(std::cout);
  timed.setType(lt_timed);
  timed.setStream(std::cout);
  start = std::chrono::high_resolution_clock::now(); 
}
  
void Logger::initFileOutput(const std::string &output)
{
  if (ParallelContext::getRank()) {
    return;
  } 
  std::string log = output + ".log";
  Logger::info << "Logs will also be printed into " << log << std::endl;
  logFile = new std::ofstream(log);
  saveLogFile = logFile;
}

void Logger::close() {
  delete(logFile);
  logFile = saveLogFile = 0;
}

