cmake_minimum_required(VERSION 2.8)

project(TwoPhaseFlow)

add_executable(twophase main.cpp)

target_include_directories(twophase PUBLIC ../INMOST/Source/Headers)
target_link_libraries(twophase ../../INMOST/build/inmost)

if(USE_SOLVER)
  if(USE_SOLVER_PETSC)
    message("linking twophase with PETSc")
    target_link_libraries(twophase ${PETSC_LIBRARIES})
  endif()
  if(USE_SOLVER_TRILINOS)
    message("linking twophase with Trilinos")
    target_link_libraries(twophase ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
  endif()
endif()


if(USE_PARTITIONER)
  if(USE_PARTITIONER_ZOLTAN)
    message("linking twophase with Zoltan")
    target_link_libraries(twophase ${ZOLTAN_LIBRARIES})
  endif()
  if(USE_PARTITIONER_PARMETIS)
    message("linking twophase with ParMETIS")
    target_link_libraries(twophase ${PARMETIS_LIBRARIES})
  endif()
endif()

if(USE_MPI)
  message("linking twophase with MPI")
  target_link_libraries(twophase ${MPI_LIBRARIES})
  if(MPI_LINK_FLAGS)
    set_target_properties(twophase PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif()
endif(USE_MPI)

#install(TARGETS twophase EXPORT inmost-targets RUNTIME DESTINATION bin)
#install(TARGETS MFDELAST EXPORT inmost-targets RUNTIME DESTINATION bin)

