#pragma once

#include <string> 
#include <sys/types.h> 
#include <sys/stat.h> 

#include <ParallelContext.hpp>
#include <IO/Logger.hpp>

class FileSystem {
public:

  static void mkdir(const std::string &path, bool masterRankOnly) {
    if (masterRankOnly && ParallelContext::getRank() != 0) {
      return;
    }
#if defined(_WIN32)
    _mkdir(path.c_str()); // can be used on Windows
#else 
    mode_t nMode = 0733; // UNIX style permissions
    ::mkdir(path.c_str(),nMode); // can be used on non-Windows
#endif

  }

  static std::string joinPaths(const std::string &p1, const std::string &p2) {
    std::string sep = "/";
#ifdef _WIN32
    std::string sep = "\\";
#endif
    return p1 + sep + p2;
  }

  static void getFileContent(const std::string &filePath, std::string &content) 
  {
    std::ifstream ifs;
    content.assign((std::istreambuf_iterator<char>(ifs)),
      (std::istreambuf_iterator<char>()) );
  }

};

