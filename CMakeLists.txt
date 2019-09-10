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
	src/stripes/Vertex.cpp
	src/stripes/Vertex.h
)

target_include_directories(stripes PUBLIC src/)
target_compile_features(stripes PUBLIC cxx_std_17)
target_compile_definitions(stripes PUBLIC SP_COMPILE_COMMAND_LINE)

# Local folders
set(SHADER_FOLDER "${CMAKE_CURRENT_SOURCE_DIR}/shaders/")
target_compile_definitions(stripes PUBLIC "-DSHADER_FOLDER=\"${SHADER_FOLDER}\"")
set(DATA_FOLDER "${CMAKE_CURRENT_SOURCE_DIR}/data/")
target_compile_definitions(stripes PUBLIC "-DDATA_FOLDER=\"${DATA_FOLDER}\"")

# Set dependencies
target_include_directories(stripes PUBLIC ${CHOLMOD_INCLUDE_DIR})
target_link_libraries(stripes PUBLIC ${CHOLMOD_LIBRARY})
target_include_directories(stripes PUBLIC ${SUITESPARSEQR_INCLUDE_DIR})
target_link_libraries(stripes PUBLIC ${SUITESPARSEQR_LIBRARY})
target_include_directories(stripes PUBLIC ${UMFPACK_INCLUDES})
target_link_libraries(stripes PUBLIC ${UMFPACK_LIBRARIES})

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

target_link_libraries(StripePatterns PUBLIC stripes GLUT::GLUT glad::glad OpenGL::GLU)
