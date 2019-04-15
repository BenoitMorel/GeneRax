set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")


# Add _USE_MATH_DEFINES for cmath.
add_definitions(-D_USE_MATH_DEFINES)
set(MPISCHEDULER_BUILD_AS_LIBRARY True)
set(MPISCHEDULER_STATIC_SCHEDULED_MAIN True)
# Add Bio++ and Boost sources and create a libraries
add_subdirectory(pll-modules)
add_subdirectory(MPIScheduler)

set(PLLMOD_EXPORTED_INCLUDE_PATH ${PLLMOD_EXPORTED_INCLUDE_PATH} ${PLL_EXPORTED_INCLUDE_PATH} PARENT_SCOPE)
