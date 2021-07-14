set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")


# Add _USE_MATH_DEFINES for cmath.
add_definitions(-D_USE_MATH_DEFINES)
set(MPISCHEDULER_BUILD_AS_LIBRARY True)
set(MPISCHEDULER_STATIC_SCHEDULED_MAIN True)
# Dependencies
add_subdirectory(MPIScheduler)
# the two following lines disable corax tests build
set(CORAX_TESTS_ENABLED OFF)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
add_subdirectory(coraxlib)

set(PLLMOD_EXPORTED_INCLUDE_PATH ${PLLMOD_EXPORTED_INCLUDE_PATH} ${PLL_EXPORTED_INCLUDE_PATH} PARENT_SCOPE)
