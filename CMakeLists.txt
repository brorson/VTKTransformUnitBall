set(VTK_DIR /usr/local/)

cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

#-----------------------------------------------------
project(TransformUnitBall)

find_package(VTK CONFIG)
find_package(VTK REQUIRED)

if (NOT VTK_FOUND)
  message("Skipping UnitBall: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()

message (STATUS "VTK_VERSION: ${VTK_VERSION}")
message (STATUS "VTK_LIBRARIES: ${VTK_LIBRARIES}")

include_directories(
		"/usr/local/include/vtk-9.0/"
)

link_directories(
		"/usr/local/lib/"
)

# set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES "CMakeLists.txt~")

add_executable(TransformUnitBall TransformUnitBall.cxx )
target_link_libraries(TransformUnitBall ${VTK_LIBRARIES})

vtk_module_autoinit(
    TARGETS TransformUnitBall
    MODULES ${VTK_LIBRARIES} 
) 

set_property(
        TARGET TransformUnitBall
        APPEND
        PROPERTY ADDITIONAL_CLEAN_FILES *~
)

