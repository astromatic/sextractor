# CMake file for levmar's MEX-file; see http://www.cmake.org
# Requires FindMatlab.cmake included with cmake

PROJECT(LEVMARMEX)
#CMAKE_MINIMUM_REQUIRED(VERSION 1.4)

INCLUDE("C:/Program Files/CMake 2.4/share/cmake-2.4/Modules/FindMatlab.cmake")

# f2c is sometimes equivalent to libF77 & libI77; in that case, set HAVE_F2C to 0
SET(HAVE_F2C 1 CACHE BOOL "Do we have f2c or F77/I77?" )

# the directory where the lapack/blas/f2c libraries reside
SET(LAPACKBLAS_DIR /usr/lib CACHE PATH "Path to lapack/blas libraries")

# the directory where levmar.h resides
SET(LM_H_DIR .. CACHE PATH "Path to levmar.h")
# the directory where the levmar library resides
SET(LEVMAR_DIR .. CACHE PATH "Path to levmar library")

# actual names for the lapack/blas/f2c libraries
SET(LAPACK_LIB lapack CACHE STRING "The name of the lapack library")
SET(BLAS_LIB blas CACHE STRING "The name of the blas library")
IF(HAVE_F2C)
  SET(F2C_LIB f2c CACHE STRING "The name of the f2c library")
ELSE(HAVE_F2C)
  SET(F77_LIB libF77 CACHE STRING "The name of the F77 library")
  SET(I77_LIB libI77 CACHE STRING "The name of the I77 library")
ENDIF(HAVE_F2C)

########################## NO CHANGES BEYOND THIS POINT ##########################

INCLUDE_DIRECTORIES(${LM_H_DIR})
LINK_DIRECTORIES(${LAPACKBLAS_DIR} ${LEVMAR_DIR})

SET(SRC levmar.c)

# naming conventions for the generated file's suffix
IF(WIN32)
  SET(SUFFIX ".mexw32")
ELSE(WIN32)
  SET(SUFFIX ".mexglx")
ENDIF(WIN32)

SET(OUTNAME "levmar${SUFFIX}")

ADD_LIBRARY(${OUTNAME} MODULE ${SRC})

IF(HAVE_F2C)
	ADD_CUSTOM_COMMAND(OUTPUT ${OUTNAME}
                   DEPENDS ${SRC}
                   COMMAND mex
                   ARGS -O -I${LM_H_DIR} ${SRC} -I${MATLAB_INCLUDE_DIR} -L${LAPACKBLAS_DIR} -L${LEVMAR_DIR} -L${MATLAB_MEX_LIBRARY} -llevmar -l${LAPACK_LIB} -l${BLAS_LIB} -l${F2C_LIB} -output ${MATLAB_LIBRARIES} ${OUTNAME})
ELSE(HAVE_F2C)
	ADD_CUSTOM_COMMAND(OUTPUT ${OUTNAME}
                   DEPENDS ${SRC}
                   COMMAND mex
                   ARGS -O -I${LM_H_DIR} ${SRC} -I${MATLAB_INCLUDE_DIR} -L${LAPACKBLAS_DIR} -L${LEVMAR_DIR} -L${MATLAB_MEX_LIBRARY} -llevmar -l${LAPACK_LIB} -l${BLAS_LIB} -l${F77_LIB} -l${I77_LIB} ${MATLAB_LIBRARIES} -output ${OUTNAME})
ENDIF(HAVE_F2C)
