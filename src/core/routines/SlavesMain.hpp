#pragma once

extern "C" int static_scheduled_main(int argc, char** argv, void* comm);

class SlavesMain {
public:
  SlavesMain() = delete;
  static bool isSlave(int argc, char** argv);
};
