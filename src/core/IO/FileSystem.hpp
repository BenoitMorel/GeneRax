#pragma once

#include <fstream>
#include <string> 
#include <sys/types.h> 
#include <sys/stat.h> 

#include <parallelization/ParallelContext.hpp>
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


  static bool exists (const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
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
    std::ifstream ifs(filePath);
    content.assign((std::istreambuf_iterator<char>(ifs)),
      (std::istreambuf_iterator<char>()) );
  }

  static void replaceWithContentIfFile(std::string &str)
  {
    std::ifstream ifs(str);
    if (ifs.good()) {
      str.assign((std::istreambuf_iterator<char>(ifs)),
        (std::istreambuf_iterator<char>()) );
    }
  }
  
  static void copy(const std::string &f1, const std::string &f2, bool masterRankOnly) 
  {
    if (masterRankOnly && ParallelContext::getRank() != 0) {
      return;
    }
    std::ifstream  src(f1, std::ios::binary);
    std::ofstream  dst(f2,   std::ios::binary);
    dst << src.rdbuf();
  }
};

