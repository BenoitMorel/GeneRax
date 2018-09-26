# CMake script for Bio++ and Treerecs Library
# Authors:
#   Nicolas COMTE
# Created: 27/07/2017

cmake_minimum_required(VERSION 3.0.1)

project(bpp)

# Add Bpp sources
add_subdirectory(bpp-core)
add_subdirectory(bpp-seq)
add_subdirectory(bpp-phyl)
