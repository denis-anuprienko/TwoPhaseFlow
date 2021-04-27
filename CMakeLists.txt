cmake_minimum_required(VERSION 2.8)

project(TwoPhaseFlow)

option(USE_MPI "Compile with MPI support" ON)
option(USE_TRILINOS "Use Trilinos solver" OFF)

add_executable(twophase main.cpp)
add_executable(gridrefine gridrefine.cpp)
add_executable(gmshgrid gmshgrid.cpp)

find_package(INMOST REQUIRED)
if(NOT INMOST_FOUND)
    message("INMOST not found!")
endif()

include_directories(${INMOST_INCLUDE_DIRS})
message(${INMOST_INCLUDE_DIRS})
message(${INMOST_LIBRARIES})
target_link_libraries(twophase ${INMOST_LIBRARIES})
target_link_libraries(gridrefine ${INMOST_LIBRARIES})
target_link_libraries(gmshgrid ${INMOST_LIBRARIES})
#target_link_libraries(twophase ../../INMOST/build-debug/inmost)#${INMOST_LIBRARIES})
#target_link_libraries(gridrefine ../../INMOST/build-debug/inmost)#${INMOST_LIBRARIES})

if(USE_MPI)
    message("Dealing with MPI")
    find_package(MPI REQUIRED)
    if(NOT MPI_FOUND)
        message("MPI not found!")
    endif()

    message((${MPI_INCLUDE_PATH}))
    include_directories(${MPI_INCLUDE_PATH})
    target_link_libraries(twophase ${MPI_CXX_LIBRARIES})

    if(MPI_LINK_FLAGS)
        set_target_properties(twophase PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
    endif()
endif()

if(USE_TRILINOS)
    set(CMAKE_PREFIX_PATH ${TRILINOS_PATH} ${CMAKE_PREFIX_PATH})
    find_package(Trilinos PATHS ${TRILINOS_PATH}/lib/cmake/Trilinos ${TRILINOS_PATH} REQUIRED)
    if(NOT Trilinos_FOUND)
        message("Trilinos NOT FOUND")
    else()

        MESSAGE("\nFound Trilinos!  Here are the details: ")
        MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
        MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
        MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
        MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
        MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
        MESSAGE("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
        MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
        MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
        MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
        MESSAGE("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
        MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
        MESSAGE("End of Trilinos details\n")

        include_directories(${Trilinos_INCLUDE_DIRS})
        include_directories(${Trilinos_TPL_INCLUDE_DIRS})
        #link_directories(${Trilinos_LIBRARY_DIRS})
        #link_directories(${Trilinos_TPL_LIBRARY_DIRS})
        target_link_libraries(twophase ${Trilinos_LIBRARIES})
        target_link_libraries(twophase ${Trilinos_TPL_LIBRARIES})


        message("Trilinos FOUND")
    endif()
endif()

