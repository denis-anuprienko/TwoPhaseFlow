cmake_minimum_required(VERSION 2.8)

project(TwwoPhase)

add_executable(MAIN main.cpp)

target_include_directories(MAIN PUBLIC ../INMOST/Source/Headers)
target_link_libraries(MAIN ../../INMOST/build/inmost)

if(USE_SOLVER)
  if(USE_SOLVER_PETSC)
    message("linking MAIN with PETSc")
    target_link_libraries(MAIN ${PETSC_LIBRARIES})
  endif()
  if(USE_SOLVER_TRILINOS)
    message("linking MAIN with Trilinos")
    target_link_libraries(MAIN ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
  endif()
endif()


if(USE_PARTITIONER)
  if(USE_PARTITIONER_ZOLTAN)
    message("linking MAIN with Zoltan")
    target_link_libraries(MAIN ${ZOLTAN_LIBRARIES})
  endif()
  if(USE_PARTITIONER_PARMETIS)
    message("linking MAIN with ParMETIS")
    target_link_libraries(MAIN ${PARMETIS_LIBRARIES})
  endif()
endif()

if(USE_MPI)
  message("linking MAIN with MPI")
  target_link_libraries(MAIN ${MPI_LIBRARIES})
  if(MPI_LINK_FLAGS)
    set_target_properties(MAIN PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif()
endif(USE_MPI)

#install(TARGETS MAIN EXPORT inmost-targets RUNTIME DESTINATION bin)
#install(TARGETS MFDELAST EXPORT inmost-targets RUNTIME DESTINATION bin)

