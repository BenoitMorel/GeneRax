#pragma once

#include <vector>
#include <IO/Logger.hpp>
#include <IO/ParallelOfstream.hpp>

class Scenario {
public:  
  enum EventType {
    S = 0 , SL, D, T, None, Invalid
  };

  struct Event {
    EventType type;
    int geneNode;
    int speciesNode;
  };


  Scenario(const string &fileName): _eventsCount(int(Invalid), 0), _os(fileName) {}
  
  void addEvent(EventType type, int geneNode, int speciesNode);

  void saveEventsCounts() {
    for (int i = 0; i < int(Invalid); ++i) {
      _os << eventNames[i] << ":" << _eventsCount[i] << endl;
      Logger::info << _eventsCount[i] << " " << eventNames[i] << " events." << endl;
    }
  }

private:
  static const char *eventNames[];
  vector<Event> _events;
  vector<int> _eventsCount;
  ParallelOfstream _os;
};



