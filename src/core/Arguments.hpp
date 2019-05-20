#pragma once

#include <string>
#include <enums.hpp>
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

  static std::string recModelToStr(RecModel recModel)
  {
    switch(recModel) {
    case UndatedDL:
      return "UndatedDL";
    case UndatedDTL:
      return "UndatedDTL";
    };
    exit(41);
  }

  static RecModel strToRecModel(const std::string &str)
  {
    if (str == "UndatedDL") {
      return UndatedDL;
    } else if (str == "UndatedDTL") {
      return UndatedDTL;
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
    }
    exit(41);
  }

  static RecOpt strToRecOpt(const std::string &str) {
    if (str == "grid") {
      return Grid;
    } else if (str == "simplex") {
      return Simplex;
    } else {
      Logger::info << "Invalid reconciliation optimization method " << str << std::endl;
      exit(41);
    }
  }

};
