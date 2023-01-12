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
    case GeneSearchStrategy::SKIP:
      return "SKIP";
    }
    exit(41);
  }

  static GeneSearchStrategy strToStrategy(const std::string &str) 
  {
    if (str == "SPR") {
      return GeneSearchStrategy::SPR;
    }  else if (str == "EVAL") {
      return GeneSearchStrategy::EVAL;
    }  else if (str == "SKIP") {
      return GeneSearchStrategy::SKIP;
    }  else if (str == "RECONCILE") {
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
    case SpeciesSearchStrategy::SKIP:
      return "SKIP";
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
    }  else if (str == "SKIP") {
      return SpeciesSearchStrategy::SKIP;
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

  static std::string transferConstraintToStr(TransferConstaint tc) {
    switch(tc) {
    case TransferConstaint::NONE:
      return "NONE";
    case TransferConstaint::PARENTS:
      return "PARENTS";
    case TransferConstaint::SOFTDATED:
      return "SOFTDATED";
    }
    exit(41);
  }

  static TransferConstaint strToTransferConstraint(const std::string &str) {
    if (str == "NONE") {
      return TransferConstaint::NONE;
    } else if (str == "PARENTS") {
      return TransferConstaint::PARENTS;
    } else if(str == "SOFTDATED") {
      return TransferConstaint::SOFTDATED;
    } else {
      Logger::info << "Invalid transfer constraint " << str << std::endl;
      exit(41);
      return TransferConstaint::NONE;
    }
  }

  static std::string recOptToStr(RecOpt recOpt) {
    switch(recOpt) {
    case RecOpt::Grid:
      return "GRID";
    case RecOpt::Simplex:
      return "SIMPLEX";
    case RecOpt::Gradient:
      return "GRADIENT";
    case RecOpt::None:
      return "NONE";
    }
    exit(41);
  }

  static RecOpt strToRecOpt(const std::string &str) {
    if (str == "GRID") {
      return RecOpt::Grid;
    } else if (str == "SIMPLEX") {
      return RecOpt::Simplex;
    } else if (str == "GRADIENT") {
      return RecOpt::Gradient;
    } else if (str == "NONE") {
      return RecOpt::None;
    } else {
      Logger::info << "Invalid reconciliation optimization method " << str << std::endl;
      exit(41);
    }
  }

  static CCPRooting strToCCPRooting(const std::string &str) {
    if (str == "UNIFORM") {
      return CCPRooting::UNIFORM;
    } else if (str == "ROOTED") {
      return CCPRooting::ROOTED;
    } else if (str == "MAD") {
      return CCPRooting::MAD;
    } else {
      Logger::info << "Invalid rooting mode for conditional clade probabilities: " << str << std::endl;
      exit(41);
    }
  }
};
