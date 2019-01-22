#pragma once

#include <vector>
#include <IO/Logger.hpp>
#include <IO/ParallelOfstream.hpp>

class Scenario {
public:  
  enum EventType {
    D = 0 , T, S, SL, None, Invalid
  };

  static const char *eventNames[];

  Scenario(const string &fileName): _eventsCount(int(Invalid), 0), _os(fileName) {}
  
  void addEvent(EventType type, int geneNode, int speciesNode);

  void printEventsCounts() {
    for (int i = 0; i < int(Invalid); ++i) {
      Logger::info << _eventsCount[i] << " " << eventNames[i] << " events." << endl;
    }
  }

private:
  vector<int> _eventsCount;
  ParallelOfstream _os;
};



