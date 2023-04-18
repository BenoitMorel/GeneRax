
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wextra")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -Wextra -pedantic -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 -Wundef -Wno-unused -Wno-variadic-macros -Wno-parentheses -fdiagnostics-show-option")


  



macro(add_program program_name program_sources)
  if (JOINTSEARCH_BUILD_AS_LIBRARY)
    add_library(${program_name} SHARED ${program_sources})
  else()
    add_executable(${program_name} ${program_sources})
  endif()
  target_sources(${program_name} PRIVATE ${program_sources})
  target_link_libraries(${program_name} 
    PRIVATE
    generaxcore
    ${CORAX_LIB_TARGET}
    ${MPI_LIBRARIES}
    )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  if (JOINTSEARCH_BUILD_AS_LIBRARY AND NOT APPLE)
    target_link_libraries(${program_name} PRIVATE gcc_s gcc)
  endif()

  set_target_properties(${program_name}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )
  set_target_properties(${program_name} 
    PROPERTIES
    COMPILE_FLAGS "-g ${MPI_COMPILE_FLAGS}"
    )
  set_target_properties(${program_name} 
    PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}"
    )

  target_include_directories(${program_name}
    PUBLIC ${JOINTSEARCH_INCLUDE_DIRS}
    )
endmacro()

macro(add_program_corax program_name program_sources)
  if (JOINTSEARCH_BUILD_AS_LIBRARY)
    add_library(${program_name} SHARED ${program_sources})
  else()
    add_executable(${program_name} ${program_sources})
  endif()
  target_link_libraries(${program_name} 
    generaxcore
    ${CORAX_LIB_TARGET}
    ${MPI_LIBRARIES}
    )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
  if (JOINTSEARCH_BUILD_AS_LIBRARY AND NOT APPLE)
    target_link_libraries(${program_name} gcc_s gcc)
  endif()

  set_target_properties(${program_name}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )
  set_target_properties(${program_name} 
    PROPERTIES
    COMPILE_FLAGS "-g ${MPI_COMPILE_FLAGS}"
    )
  set_target_properties(${program_name} 
    PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}"
    )

  target_include_directories(${program_name}
    PUBLIC ${JOINTSEARCH_INCLUDE_DIRS}
    )
endmacro()



add_subdirectory(GeneRax)
add_subdirectory(ale)
#add_subdirectory(Asteroid)
add_subdirectory(NJRax)
add_subdirectory(LightSearch)
add_subdirectory(tools)
add_subdirectory(tests)

