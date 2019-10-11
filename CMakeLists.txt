################################################################################
cmake_minimum_required(VERSION 3.15)

include("cmake/HunterGate.cmake")
HunterGate(
	URL "https://github.com/cpp-pm/hunter/archive/v0.23.219.tar.gz"
	SHA1 "9cf34f61c8d34cfc2b1396a63aaaedaa5e09441e"
)

################################################################################

project(stripes)

if(POLICY CMP0091)
	cmake_policy(SET CMP0091 NEW)
endif()
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

get_directory_property(HAS_PARENT PARENT_DIRECTORY)
if(HAS_PARENT)
	set(STRIPES_TOPLEVEL_PROJECT OFF)
else()
	set(STRIPES_TOPLEVEL_PROJECT ON)
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(UseColors)
include(Warnings)

################################################################################

# Stripe library
add_library(stripes
	src/stripes/Complex.cpp
	src/stripes/Complex.h
	src/stripes/DenseMatrix.cpp
	src/stripes/DenseMatrix.h
	src/stripes/DenseMatrix.inl
	src/stripes/Edge.cpp
	src/stripes/Edge.h
	src/stripes/Face.cpp
	src/stripes/Face.h
	src/stripes/HalfEdge.cpp
	src/stripes/HalfEdge.h
	src/stripes/Image.cpp
	src/stripes/Image.h
	src/stripes/LinearContext.cpp
	src/stripes/LinearContext.h
	src/stripes/LinearEquation.cpp
	src/stripes/LinearEquation.h
	src/stripes/LinearPolynomial.cpp
	src/stripes/LinearPolynomial.h
	src/stripes/LinearSystem.cpp
	src/stripes/LinearSystem.h
	src/stripes/Mesh.cpp
	src/stripes/Mesh.h
	src/stripes/MeshIO.cpp
	src/stripes/MeshIO.h
	src/stripes/Quaternion.cpp
	src/stripes/Quaternion.h
	src/stripes/Real.cpp
	src/stripes/Real.h
	src/stripes/SparseMatrix.cpp
	src/stripes/SparseMatrix.h
	src/stripes/SparseMatrix.inl
	src/stripes/Types.h
	src/stripes/Utility.h
	src/stripes/Variable.cpp
	src/stripes/Variable.h
	src/stripes/Vector.cpp
	src/stripes/Vector.h
	src/stripes/Vector.inl
	src/stripes/Vertex.cpp
	src/stripes/Vertex.h
)

add_library(stripes::stripes ALIAS stripes)
target_include_directories(stripes PUBLIC src/)
target_compile_features(stripes PUBLIC cxx_std_17)
target_compile_definitions(stripes PUBLIC SP_COMPILE_COMMAND_LINE)
target_link_libraries(stripes PRIVATE warnings::all)

if(WIN32)
	target_compile_definitions(stripes PUBLIC _USE_MATH_DEFINES)
endif()

# Local folders
set(SHADER_FOLDER "${CMAKE_CURRENT_SOURCE_DIR}/shaders/")
target_compile_definitions(stripes PUBLIC "-DSHADER_FOLDER=\"${SHADER_FOLDER}\"")
set(DATA_FOLDER "${CMAKE_CURRENT_SOURCE_DIR}/data/")
target_compile_definitions(stripes PUBLIC "-DDATA_FOLDER=\"${DATA_FOLDER}\"")

# Dependencies
hunter_add_package(SuiteSparse)
find_package(SuiteSparse CONFIG REQUIRED)
target_link_libraries(stripes PUBLIC SuiteSparse::cholmod SuiteSparse::umfpack SuiteSparse::spqr)

################################################################################

find_package(Eigen3 QUIET NO_MODULE)
if(TARGET Eigen3::Eigen)
	# Libigl-style wrapper using Eigen types
	add_library(stripes_wrapper
		src/stripes/wrapper/Stripes.cpp
		src/stripes/wrapper/Stripes.h
	)
	add_library(stripes::wrapper ALIAS stripes_wrapper)

	# Dependencies
	target_link_libraries(stripes_wrapper PUBLIC stripes::stripes Eigen3::Eigen)
else()
	message(WARNING "Eigen3::Eigen target not available. stripes::wrapper will not be built.")
endif()

################################################################################

if(STRIPES_TOPLEVEL_PROJECT)
	find_package(OpenGL)
	find_package(GLUT)
	if(NOT (OPENGL_GLU_FOUND AND GLUT_FOUND))
		message(WARNING "OpenGL or GLUT could not be found. The viewer will not be built.")
		return()
	endif()

	# Demo application
	add_executable(StripePatterns
		src/viewer/main.cpp
		src/viewer/Camera.cpp
		src/viewer/Camera.h
		src/viewer/Shader.cpp
		src/viewer/Shader.h
		src/viewer/Viewer.cpp
		src/viewer/Viewer.h
	)

	# Dependencies
	add_subdirectory(ext/glad)
	target_link_libraries(StripePatterns
		PUBLIC
			stripes::stripes
			GLUT::GLUT
			glad::glad
			OpenGL::GLU
	)
endif()
