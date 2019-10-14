hunter_config(SuiteSparse
    URL "https://github.com/jdumas/suitesparse-metis-for-windows.git"
    SHA1 "bf5e0b8565ffcbdc887f7f7933f6a17f47523064"
    CMAKE_ARGS
        USE_MKL=ON
        MKL_ROOT=${CMAKE_CURRENT_LIST_DIR}/ext/mkl
)
