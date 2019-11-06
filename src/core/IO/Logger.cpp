#include "IO/Logger.hpp"


Logger Logger::info;
Logger Logger::error;
Logger Logger::timed;
Logger Logger::perrank;
TimePoint Logger::start;
std::string Logger::outputdir;
std::ofstream *Logger::logFile = nullptr;
std::ofstream *Logger::rankLogFile = nullptr;
std::ofstream *Logger::saveLogFile = nullptr;
bool Logger::inited = false;

Logger::Logger(): _os(&std::cout), _silent(false) {
  setType(lt_info);  
}

void Logger::init() {
  if (inited) {
    return;
  }
  inited = true;
  info.setType(lt_info);
  error.setStream(std::cout);
  error.setType(lt_error);
  error.setStream(std::cout);
  timed.setType(lt_timed);
  timed.setStream(std::cout);
  perrank.setType(lt_perrank);
  start = std::chrono::high_resolution_clock::now(); 
  ParallelContext::barrier();
}

void Logger::initRankFileOutput()
{
  if (!rankLogFile) {
    std::string rankLog = outputdir + "rank_" + std::to_string(ParallelContext::getRank())
      + ".log";
    rankLogFile = new std::ofstream(rankLog);
  }
}


void Logger::initFileOutput(const std::string &output)
{
  Logger::outputdir = output;
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
  logFile = saveLogFile = nullptr;
}

