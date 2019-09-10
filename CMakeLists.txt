cmake_minimum_required(VERSION 3.11)
project(stripes)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(UseColors)

# Dependencies
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(SuiteSparse REQUIRED)
find_package(Umfpack REQUIRED)
add_subdirectory(ext/glad)

# Stripe library
add_library(stripes
	include/Camera.h
	include/Complex.h
	include/DenseMatrix.h
	include/Edge.h
	include/Face.h
	include/HalfEdge.h
	include/Image.h
	include/LinearContext.h
	include/LinearEquation.h
	include/LinearPolynomial.h
	include/LinearSystem.h
	include/Mesh.h
	include/MeshIO.h
	include/Quaternion.h
	include/Real.h
	include/Shader.h
	include/SparseMatrix.h
	include/Types.h
	include/Utility.h
	include/Variable.h
	include/Vector.h
	include/Vertex.h
	include/Viewer.h
	src/Camera.cpp
	src/Complex.cpp
	src/DenseMatrix.cpp
	src/DenseMatrix.inl
	src/Edge.cpp
	src/Face.cpp
	src/HalfEdge.cpp
	src/Image.cpp
	src/LinearContext.cpp
	src/LinearEquation.cpp
	src/LinearPolynomial.cpp
	src/LinearSystem.cpp
	src/Mesh.cpp
	src/MeshIO.cpp
	src/Quaternion.cpp
	src/Real.cpp
	src/Shader.cpp
	src/SparseMatrix.cpp
	src/SparseMatrix.inl
	src/Variable.cpp
	src/Vector.cpp
	src/Vertex.cpp
	src/Viewer.cpp
)

target_include_directories(stripes PUBLIC include/ src/)
target_compile_features(stripes PUBLIC cxx_std_17)

# Shader folder
set(SHADER_FOLDER "${CMAKE_CURRENT_SOURCE_DIR}/shaders/")
target_compile_definitions(stripes PUBLIC "-DSHADER_FOLDER=\"${SHADER_FOLDER}\"")

# Set dependencies
target_link_libraries(stripes PUBLIC GLUT::GLUT glad::glad OpenGL::GLU)
target_include_directories(stripes PUBLIC ${CHOLMOD_INCLUDE_DIR})
target_link_libraries(stripes PUBLIC ${CHOLMOD_LIBRARY})
target_include_directories(stripes PUBLIC ${SUITESPARSEQR_INCLUDE_DIR})
target_link_libraries(stripes PUBLIC ${SUITESPARSEQR_LIBRARY})
target_include_directories(stripes PUBLIC ${UMFPACK_INCLUDES})
target_link_libraries(stripes PUBLIC ${UMFPACK_LIBRARIES})

# Demo application
add_executable(StripePatterns src/main.cpp)
target_link_libraries(StripePatterns PUBLIC stripes)
