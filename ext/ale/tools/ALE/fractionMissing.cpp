// Modified by N.Comte for Treerecs – Copyright © by INRIA – All rights reserved – 2017
#include "fractionMissing.h"
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Text/TextTools.h>

inline bool fexists(const std::string& filename)
{
  std::ifstream ifile(filename);
  return ifile.good();
}

std::map<std::string, scalar_type> readFractionMissingFile(std::string fractionMissingFile) {
  std::map<std::string, scalar_type> toReturn;
  if (fractionMissingFile=="" )
  {
    std::cout << "No file providing the fraction of missing genes per species, we assume that all species have 100% of their genes."<<std::endl;
  }
  else
  {
    if (!fexists(fractionMissingFile)) {
      std::cout << "Error, file "<< fractionMissingFile << " does not seem accessible." <<  std::endl;
      exit(1);
    }
    std::ifstream inCoverage (fractionMissingFile.c_str());
    std::vector <std::string> listCoverages;
    std::string line;

    while(getline(inCoverage,line))
    {
      listCoverages.push_back(line);
    }
    for(std::vector<std::string >::iterator it = listCoverages.begin(); it != listCoverages.end(); it++)
    {
      bpp::StringTokenizer st1(*it, ":", true);
      toReturn[st1.getToken(0)] = bpp::TextTools::toDouble ( st1.getToken( 1 ) ) ;
    }
  }
  return (toReturn);

}
