// Treerecs – Copyright © by INRIA – All rights reserved – 2017
// Created by Nicolas Comte on 28/12/16.
//

#ifndef PHYLASOLVER_COST_H
#define PHYLASOLVER_COST_H

#include <ostream>
#include <vector>
#include <Constants.h>

/// Types of events in genes history.
enum Event { duplication /// Duplication.
              , loss /// Loss.
              , speciation /// Speciation.
              , speciationLoss /// Speciation resulting in one loss son.
              , extant /// Extant.
              , bifurcationOut /// Bifurcation out.
              , none /// Nothing/ Unknown.
              };

/// Get event in std::string.
inline std::string event_to_str(const Event& event) {
  switch (event) {
    case duplication:
      return DUPLICATION_STR_FLAG;
    case loss:
      return LOSS_STR_FLAG;
    case speciation:
      return SPECIATION_STR_FLAG;
    case extant:
      return EXTANT_STR_FLAG;
    case bifurcationOut:
      return BIFURCATION_OUT_STR_FLAG;
    case speciationLoss:
      return SPECIATION_LOSS_STR_FLAG;
    default:
      return "none";
  }
}

inline std::ostream& operator<< (std::ostream& os, const Event& cost_type) {
  os << event_to_str(cost_type);
  return os;
}

typedef struct CostStruct {
  /*!
   * \struct CostStruct
   * \brief This structure contains cost informations for a Cost Table: the value of the cost and what (enum Event :
   * loss or dup, spec, ..., none defaults).
   */
  double value = 0.0;
  std::vector<Event> path = {none};
  std::vector<std::size_t> occurence = {1};
  friend std::ostream& operator<<(std::ostream& os, const struct CostStruct& cost);
} Cost;

inline std::ostream& operator<< (std::ostream& os, const Cost& cost) {
  /*!
   * Streams all pairs of an unordered_map (std).
   */
  os << "{";
  os << cost.value;
  os << " : ";
  if(cost.path.empty())
    os << " ";
  else
    for(std::size_t i = 0; i < cost.path.size(); i++){
      auto& type = cost.path[i];
      auto& n = cost.occurence[i];
      os << type << ":" + std::to_string(n);
      if(i != cost.path.size() - 1)
        os << ", ";
    }
  os << "}";
  return os;
}

#endif //PHYLASOLVER_COST_H
