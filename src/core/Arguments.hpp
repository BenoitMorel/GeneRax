#pragma once

#include <string>
#include <enums.hpp>
#include <IO/Logger.hpp>

using namespace std;

class Arguments {
public:
  static string strategyToStr(Strategy s)
  {
    switch(s) {
    case SPR:
      return "SPR";
    case EVAL:
      return "EVAL";
    }
    exit(41);
  }

  static Strategy strToStrategy(const string &str) 
  {
    if (str == "SPR") {
      return SPR;
    }  else if (str == "EVAL") {
      return EVAL;
    } else {
      Logger::info << "Invalid strategy " << str << endl;
      exit(41);
    }
  }

  static string recModelToStr(RecModel recModel)
  {
    switch(recModel) {
    case UndatedDL:
      return "UndatedDL";
    case UndatedDTL:
      return "UndatedDTL";
    case DatedDL:
      return "DatedDL";
    };
    exit(41);
  }

  static RecModel strToRecModel(const string &str)
  {
    if (str == "UndatedDL") {
      return UndatedDL;
    } else if (str == "UndatedDTL") {
      return UndatedDTL;
    } else if (str == "DatedDL") {
      return DatedDL;
    } else {
      Logger::info << "Invalid reconciliation model " << str << endl;
      exit(41);
    }
  }

  static string recOptToStr(RecOpt recOpt) {
    switch(recOpt) {
    case Grid:
      return "grid";
    case Simplex:
      return "simplex";
    }
    exit(41);
  }

  static RecOpt strToRecOpt(const string &str) {
    if (str == "grid") {
      return Grid;
    } else if (str == "simplex") {
      return Simplex;
    } else {
      Logger::info << "Invalid reconciliation optimization method " << str << endl;
      exit(41);
    }
  }

};
