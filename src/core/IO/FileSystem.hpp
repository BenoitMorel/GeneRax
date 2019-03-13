#pragma once

#include <string> 
#include <sys/types.h> 
#include <sys/stat.h> 

#include <ParallelContext.hpp>
#include <IO/Logger.hpp>

class FileSystem {
public:

  static void mkdir(const string &path, bool masterRankOnly) {
    if (masterRankOnly && ParallelContext::getRank() != 0) {
      return;
    }
    int nError = 1;
#if defined(_WIN32)
    nError = _mkdir(path.c_str()); // can be used on Windows
#else 
    mode_t nMode = 0733; // UNIX style permissions
    nError = ::mkdir(path.c_str(),nMode); // can be used on non-Windows
#endif
    if (nError != 0) {
      Logger::info << "Failed to create directory " << path << endl; 
      ParallelContext::abort(1);
    }
  }

  static string joinPaths(const string &p1, const string &p2) {
    string sep = "/";
#ifdef _WIN32
    string sep = "\\";
#endif
    return p1 + sep + p2;
  }
};

