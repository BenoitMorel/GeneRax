#ifndef _JOINTSEQRCH_LOGGER_HPP_
#define _JOINTSEQRCH_LOGGER_HPP_

#include <chrono>
#include <iostream>
#include <fstream>
using namespace std;

using TimePoint = chrono::high_resolution_clock::time_point;

class Logger: public ofstream
{
private:
  enum LoggerType {
    lt_info, lt_error, lt_timed, lt_silent
  };
  LoggerType _type;
  ostream *_os;
  
  Logger();

  void setStream(ostream &os) {_os = &os;}
  void setType(LoggerType type) {_type = type;}


public:
  static void init();

  template <typename T>
    Logger& operator<<(T&& t)
    {
      if (_type == lt_silent) {
        return *this;
      } else if (_type == lt_timed) {
        auto finish = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = finish - start;
        auto seconds = chrono::duration_cast<chrono::seconds>(elapsed).count();
        auto hours  = seconds / 3600;
        auto minutes = (seconds % 3600) / 60;
        seconds = seconds % 60;
        char s[25];
        sprintf(s, "%02d:%02d:%02d", hours, minutes, seconds);
        *_os << "[" << s << "] ";
        *_os << t;
        return Logger::info;
      } else {
        *_os << t;
        return *this;
      }
    }
    
  Logger& operator<<(std::ostream & (*manip)(std::ostream &)) 
  {
    if (_type != lt_silent) {
      manip(*_os);
    }
    return *this;
  }

  static Logger info;
  static Logger error;
  static Logger timed;
  static TimePoint start;
};

#endif

