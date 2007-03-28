# LEVMAR cmake file (http://www.cmake.org)

#CMAKE_MINIMUM_REQUIRED(VERSION 1.4) 

PROJECT(LEVMAR)

# ROOTDIR should be set to the directory of your compiler lapack libs
IF (WIN32)
SET(ROOTDIR C:/local/lib)
ELSE (WIN32)
SET(ROOTDIR /usr/local/lib)
ENDIF (WIN32)

ADD_EXECUTABLE(lmdemo lmdemo_sources)

SOURCE_FILES(lmdemo_sources
	lm.c Axb.c lmdemo.c lm.h
)

#INCLUDE_DIRECTORIES(${ROOTDIR}/../include)

LINK_DIRECTORIES(${ROOTDIR})

TARGET_LINK_LIBRARIES(lmdemo clapack.lib blas.lib)
