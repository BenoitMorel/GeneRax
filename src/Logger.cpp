#include "Logger.hpp"
#include "ParallelContext.hpp"

Logger Logger::info;
Logger Logger::error;
Logger Logger::timed;
TimePoint Logger::start;

Logger::Logger(): _os(&cout) {
  setType(lt_info);  
}


void Logger::init() {
  if (ParallelContext::getRank()) {
    info.setType(lt_silent);
    error.setType(lt_silent);
    timed.setType(lt_silent);
    return;
  }
  info.setType(lt_info);
  error.setStream(cout);
  error.setType(lt_error);
  error.setStream(cerr);
  timed.setType(lt_timed);
  timed.setStream(cout);
  start = chrono::high_resolution_clock::now(); 
}

