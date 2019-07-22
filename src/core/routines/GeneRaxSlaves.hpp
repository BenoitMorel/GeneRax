#pragma once

extern "C" int static_scheduled_main(int argc, char** argv, void* comm);
class GeneRaxSlaves {
public:
  static bool is_slave(int argc, char** argv);
};

