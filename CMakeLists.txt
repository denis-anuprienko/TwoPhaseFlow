cmake_minimum_required(VERSION 2.8)

project(TwoPhaseFlow)

add_executable(twophase main.cpp)
add_executable(gridgen gridgen.cpp)

find_package(INMOST REQUIRED)
if(NOT INMOST_FOUND)
    message("INMOST not found!")
endif()

include_directories(${INMOST_INCLUDE_DIRS})
target_link_libraries(twophase ${INMOST_LIBRARIES})
target_link_libraries(gridgen ${INMOST_LIBRARIES})

if(USE_MPI)
  message("Dealing with MPI")
  find_package(MPI REQUIRED)
  if(NOT MPI_FOUND)
      message("MPI not found!")
  endif()

  include_directories(${MPI_INCLUDE_PATH})
  target_link_libraries(twophase ${MPI_CXX_LIBRARIES})
  target_link_libraries(gridgen ${MPI_CXX_LIBRARIES})

  if(MPI_LINK_FLAGS)
      set_target_properties(twophase PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
      set_target_properties(gridgen PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif()
endif(USE_MPI)

