#pragma once

#include <chrono>
#include <iostream>
#include <fstream>
#include "ParallelContext.hpp"

using namespace std;

using TimePoint = chrono::high_resolution_clock::time_point;

/**
 *  All logs are printed to cout.
 *  In parallel runs, only the master can log. 
 *  init() must be called before using the Logger
 *  if initFileOutput is called, logs are also printed
 *  to the given file
 *
 *  Logger::info is used for normal logs
 *  Logger::error is used for errors and prefixes messages with [Error]
 *  Logger::time also pints the elapsed time since the start of the programm
 *
 */

class Logger: public ofstream
{
private:
  enum LoggerType {
    lt_info, lt_error, lt_timed
  };
  LoggerType _type;
  ostream *_os; 
  Logger();

  void setStream(ostream &os) {_os = &os;}
  void setType(LoggerType type) {_type = type;}


public:
  static void init();
  
  static void initFileOutput(const string &output);

  bool isSilent() {
    return (_type == lt_timed || _type == lt_info) && ParallelContext::getRank(); 
  }

  template <typename T>
    Logger& operator<<(T&& t)
    {
      if (isSilent()) {
        return *this;
      }

      if (_type == lt_timed) {
        auto finish = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = finish - start;
        int seconds = chrono::duration_cast<chrono::seconds>(elapsed).count();
        int hours  = seconds / 3600;
        int minutes = (seconds % 3600) / 60;
        seconds = seconds % 60;
        char s[25];
        sprintf(s, "%02d:%02d:%02d", hours, minutes, seconds);
        *_os << "[" << s << "] " << t;
        if (logFile) 
          *logFile << "[" << s << "] " << t;
        return Logger::info;
      } else if (_type == lt_error) {
        *_os << "[Error] " << t;
        if (logFile)
          *logFile << "[Error] " << t;
        return Logger::info;
      } else {
        *_os << t;
        _os->flush();
        if (logFile) 
          *logFile << t;
        return *this;
      }
    }
    
  Logger& operator<<(std::ostream & (*manip)(std::ostream &)) 
  {
    if (isSilent()) {
      return *this;
    }
    manip(*_os);
    if (logFile) 
      manip(*logFile);
    return *this;
  }

  static Logger info;
  static Logger error;
  static Logger timed;
  static TimePoint start;
  static ofstream *logFile;
};

