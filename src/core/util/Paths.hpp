#pragma once


class Paths {
public: 

  static std::string getFile(const std::string &outputDir,
      const std::string &file) {
    return joinPaths(outputDir, file);
  }

  static std::string getTempFile(const std::string &outputDir, 
      unsigned int index) {
    return joinPaths(outputDir, "tmp_" + std::to_string(index));

  }
  
  static std::string getSpeciesTreesDir(const std::string &outputDir) {
    return joinPaths(outputDir, "species_trees");
  }

  static std::string getSpeciesTreeFile(const std::string &outputDir,
      const std::string &file) {
    return joinPaths(getSpeciesTreesDir(outputDir), file);
  }


  static std::vector<std::string> getDirectoriesToCreate(const std::string &outputDir) {
    std::vector<std::string> dirs;
    dirs.push_back(getSpeciesTreesDir(outputDir));
    return dirs;
  }
private:
  static std::string joinPaths(const std::string &p1, const std::string &p2) {
    std::string sep = "/";
#ifdef _WIN32
    std::string sep = "\\";
#endif
    return p1 + sep + p2;
  }
  
};
