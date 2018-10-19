
set(ale_SOURCE_FILES
  containers/AbstractGrid.cpp
  containers/Cost.cpp
  containers/GeneMap.cpp
  containers/NodeProperty.cpp
  tools/SpeciesGeneMapper.cpp
  tools/Utils.cpp
  tools/PhyloTreeToolBox.cpp
  tools/IO/IO.cpp
  tools/IO/Newick.cpp
  tools/IO/Nhx.cpp
  tools/ALE/ALE.cpp
  tools/ALE/exODT_DL.cpp
  )

add_library(ale
  ${ale_SOURCE_FILES}
  )
target_include_directories(ale
  PUBLIC ${JOINTSEARCH_INCLUDE_DIRS}
  PUBLIC ${BOOST_INCLUDE_DIRS}
  PUBLIC ${BPP-PHYL_INCLUDE_DIRS}
  PUBLIC ${PLLMOD_EXPORTED_INCLUDE_PATH}
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
  )

set(PLL_LIBRARIES 
  pllmodalgorithm_static
  pllmodbinary_static
  pllmodmsa_static
  pllmodoptimize_static
  pllmodtree_static
  pllmodutil_static
)

target_link_libraries(ale
  PUBLIC bpp-phyl
  PUBLIC ${PLL_LIBRARIES})
