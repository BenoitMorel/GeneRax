# Treerecs – Copyright © by INRIA – All rights reserved – 2017
# Comte Nicolas, 2017

# Add _USE_MATH_DEFINES for cmath.
add_definitions(-D_USE_MATH_DEFINES)

# Add Bio++ and Boost sources and create a libraries
add_subdirectory(pll-modules)

set(PLLMOD_EXPORTED_INCLUDE_PATH ${PLLMOD_EXPORTED_INCLUDE_PATH} ${PLL_EXPORTED_INCLUDE_PATH} PARENT_SCOPE)
