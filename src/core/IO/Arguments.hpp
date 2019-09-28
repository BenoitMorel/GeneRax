#pragma once

#include <string>
#include <util/enums.hpp>
#include <IO/Logger.hpp>



class Arguments {
public:
  static std::string strategyToStr(Strategy s)
  {
    switch(s) {
    case SPR:
      return "SPR";
    case EVAL:
      return "EVAL";
    }
    exit(41);
  }

  static Strategy strToStrategy(const std::string &str) 
  {
    if (str == "SPR") {
      return SPR;
    }  else if (str == "EVAL") {
      return EVAL;
    } else {
      Logger::info << "Invalid strategy " << str << std::endl;
      exit(41);
    }
  }
 
  static std::string speciesRaxStrategyToStr(SpeciesRaxStrategy s)
  {
    switch(s) {
    case SIMPLE_SEARCH:
      return "SimpleSearch";
    }
    exit(41);
  }

  static SpeciesRaxStrategy strToSpeciesRaxStrategy(const std::string &str) 
  {
    if (str == "SimpleSearch") {
      return SIMPLE_SEARCH;
    } else {
      Logger::info << "Invalid strategy " << str << std::endl;
      exit(41);
    }
  }

  static std::string recModelToStr(RecModel recModel)
  {
    switch(recModel) {
    case UndatedDL:
      return "UndatedDL";
    case UndatedDTL:
      return "UndatedDTL";
    case UndatedDTLAdvanced:
      return "UndatedDTLAdvanced";
    };
    exit(41);
  }

  static bool isValidRecModel(const std::string &str)
  {
    return (str == "UndatedDL") || (str == "UndatedDTL") || (str == "UndatedDTLAdvanced");
  }

  static RecModel strToRecModel(const std::string &str)
  {
    if (str == "UndatedDL") {
      return UndatedDL;
    } else if (str == "UndatedDTL") {
      return UndatedDTL;
    } else if (str == "UndatedDTLAdvanced") {
      return UndatedDTLAdvanced;
    } else {
      Logger::info << "Invalid reconciliation model " << str << std::endl;
      exit(41);
    }
  }

  static std::string recOptToStr(RecOpt recOpt) {
    switch(recOpt) {
    case Grid:
      return "grid";
    case Simplex:
      return "simplex";
    case Gradient:
      return "gradient";
    }
    exit(41);
  }

  static RecOpt strToRecOpt(const std::string &str) {
    if (str == "grid") {
      return Grid;
    } else if (str == "simplex") {
      return Simplex;
    } else if (str == "gradient") {
      return Gradient;
    } else {
      Logger::info << "Invalid reconciliation optimization method " << str << std::endl;
      exit(41);
    }
  }

};
