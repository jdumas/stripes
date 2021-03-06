################################################################################
cmake_minimum_required(VERSION 3.11)
project(stripes)

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
add_library(stripes STATIC
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

set_target_properties(stripes PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(WIN32)
	target_compile_definitions(stripes PUBLIC _USE_MATH_DEFINES)
endif()

# Local folders
set(SHADER_FOLDER "${CMAKE_CURRENT_SOURCE_DIR}/shaders/")
target_compile_definitions(stripes PUBLIC "-DSHADER_FOLDER=\"${SHADER_FOLDER}\"")
set(DATA_FOLDER "${CMAKE_CURRENT_SOURCE_DIR}/data/")
target_compile_definitions(stripes PUBLIC "-DDATA_FOLDER=\"${DATA_FOLDER}\"")

# Dependencies
if(TARGET SuiteSparse::cholmod)
	target_link_libraries(stripes PUBLIC SuiteSparse::cholmod)
else()
	if(DEFINED ENV{CONDA_PREFIX})
		set(ENV{SUITESPARSE_ROOT} $ENV{CONDA_PREFIX})
	endif()
	find_package(SuiteSparse REQUIRED)
	find_package(Umfpack REQUIRED)
	target_include_directories(stripes PUBLIC ${CHOLMOD_INCLUDE_DIR})
	target_link_libraries(stripes PUBLIC ${CHOLMOD_LIBRARY})
endif()

if(STRIPES_TOPLEVEL_PROJECT)
	# Dummy executable for debugging
	add_executable(dummy src/main.cpp)
	target_link_libraries(dummy PUBLIC stripes::stripes)
endif()

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
	if(APPLE)
		# Workaround bug in CMake FindGLUT?
		# https://gitlab.kitware.com/cmake/cmake/-/issues/20184
		list(INSERT CMAKE_FRAMEWORK_PATH 0 /System/Library/Frameworks)
	endif()

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
