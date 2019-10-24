#pragma once

#include <chrono>
#include <iostream>
#include <fstream>
#include "parallelization/ParallelContext.hpp"



using TimePoint = std::chrono::high_resolution_clock::time_point;

/**
 *  All logs are printed to std::cout.
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

class Logger: public std::ofstream
{
private:
  enum LoggerType {
    lt_info, lt_error, lt_timed, lt_perrank
  };
  LoggerType _type;
  std::ostream *_os; // I not own this one
  bool _silent;
  Logger();
  void setStream(std::ostream &os) {_os = &os;}
  void setType(LoggerType type) {_type = type;}


public:
  static void init();
  static void close();

  static void initFileOutput(const std::string &output);
  
  static void mute() {info._silent = timed._silent = true;}
  static void unmute() {info._silent = timed._silent = false;}

  bool isSilent() {
    if (_silent) {
      return true;
    }
    if (_type == lt_perrank) {
      return false;
    }
    return (_type == lt_timed || _type == lt_info) && ParallelContext::getRank(); 
  }

  static void enableLogFile(bool enable) {
    logFile = enable ? saveLogFile : 0;
  }

  static long getElapsedSec() {
      auto finish = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> elapsed = finish - start;
      return std::chrono::duration_cast<std::chrono::seconds>(elapsed).count();
  }

  template <typename T>
    Logger& operator<<(T&& t)
    {
      if (isSilent()) {
        return *this;
      }
      if (_type == lt_timed) {
        auto seconds = getElapsedSec();
        auto hours  = seconds / 3600;
        auto minutes = (seconds % 3600) / 60;
        seconds = seconds % 60;
        char s[25];
        sprintf(s, "%02ld:%02ld:%02ld", hours, minutes, seconds);
        *_os << "[" << s << "] " << t;
        if (logFile) 
          *logFile << "[" << s << "] " << t;
        return Logger::info;
      } else if (_type == lt_error) {
        *_os << "[Error] " << t;
        if (logFile)
          *logFile << "[Error] " << t;
        return Logger::info;
      } else if (_type == lt_perrank) {
        initRankFileOutput();
        *rankLogFile << t;
        return *this;
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
    if (_type != lt_perrank) {
      manip(*_os);
    }
    if (logFile && _type != lt_perrank) { 
      manip(*logFile);
    }
    if (_type == lt_perrank) {
      manip(*rankLogFile);
    }
    return *this;
  }

  static void initRankFileOutput();

  static Logger info;
  static Logger error;
  static Logger timed;
  static Logger perrank;
  static TimePoint start;
  static std::string outputdir;
  static std::ofstream *logFile;
  static std::ofstream *rankLogFile;
  static std::ofstream *saveLogFile;
  static bool inited;
};

