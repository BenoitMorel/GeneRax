#pragma once

#include <string>
#include <util/enums.hpp>
#include <IO/Logger.hpp>



/**
 * Static methods for helping argument parsing
 */
class ArgumentsHelper {
public:
  ArgumentsHelper() = delete;

  static std::string strategyToStr(GeneSearchStrategy s)
  {
    switch(s) {
    case GeneSearchStrategy::SPR:
      return "SPR";
    case GeneSearchStrategy::EVAL:
      return "EVAL";
    }
    exit(41);
  }

  static GeneSearchStrategy strToStrategy(const std::string &str) 
  {
    if (str == "SPR") {
      return GeneSearchStrategy::SPR;
    }  else if (str == "EVAL") {
      return GeneSearchStrategy::EVAL;
    } else {
      Logger::info << "Invalid strategy " << str << std::endl;
      exit(41);
    }
  }
  
  static std::string speciesStrategyToStr(SpeciesSearchStrategy s)
  {
    switch(s) {
    case SpeciesSearchStrategy::SPR:
      return "SPR";
    case SpeciesSearchStrategy::TRANSFERS:
      return "EVAL";
    case SpeciesSearchStrategy::HYBRID:
      return "HYBRID";
    case SpeciesSearchStrategy::REROOT:
      return "REROOT";
    case SpeciesSearchStrategy::EVAL:
      return "EVAL";
    }
    exit(41);
  }

  static SpeciesSearchStrategy strToSpeciesSearchStrategy(const std::string &str) 
  {
    if (str == "SPR") {
      return SpeciesSearchStrategy::SPR;
    }  else if (str == "TRANSFERS") {
      return SpeciesSearchStrategy::TRANSFERS;
    }  else if (str == "HYBRID") {
      return SpeciesSearchStrategy::HYBRID;
    }  else if (str == "REROOT") {
      return SpeciesSearchStrategy::REROOT;
    }  else if (str == "EVAL") {
      return SpeciesSearchStrategy::EVAL;
    } else {
      Logger::info << "Invalid species strategy " << str << std::endl;
      exit(41);
    }
  }
 
  static std::string recModelToStr(RecModel recModel)
  {
    switch(recModel) {
    case RecModel::ParsimonyD:
      return "ParsimonyD";
    case RecModel::UndatedDL:
      return "UndatedDL";
    case RecModel::UndatedDTL:
      return "UndatedDTL";
    case RecModel::SimpleDS:
      return "SimpleDS";
    };
    exit(41);
  }

  static bool isValidRecModel(const std::string &str)
  {
    return (str == "UndatedDL") || 
      (str == "UndatedDTL") || 
      (str == "SimpleDS") || 
      (str == "ParsimonyD")
      ;
  }

  static RecModel strToRecModel(const std::string &str)
  {
    if (str == "UndatedDL") {
      return RecModel::UndatedDL;
    } else if (str == "ParsimonyD") {
      return RecModel::ParsimonyD;
    } else if (str == "UndatedDTL") {
      return RecModel::UndatedDTL;
    } else if (str == "SimpleDS") {
      return RecModel::SimpleDS;
    } else {
      Logger::info << "Invalid reconciliation model " << str << std::endl;
      exit(41);
    }
  }

  static std::string recOptToStr(RecOpt recOpt) {
    switch(recOpt) {
    case RecOpt::Grid:
      return "grid";
    case RecOpt::Simplex:
      return "simplex";
    case RecOpt::Gradient:
      return "gradient";
    }
    exit(41);
  }

  static RecOpt strToRecOpt(const std::string &str) {
    if (str == "grid") {
      return RecOpt::Grid;
    } else if (str == "simplex") {
      return RecOpt::Simplex;
    } else if (str == "gradient") {
      return RecOpt::Gradient;
    } else {
      Logger::info << "Invalid reconciliation optimization method " << str << std::endl;
      exit(41);
    }
  }

};
