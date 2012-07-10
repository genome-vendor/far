# Add FindTBB directory to CMake's module path
# set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${FAR_SOURCE_DIR}/FindTBB/CMake/FindTBB/")


INCLUDE_DIRECTORIES(${FAR_SOURCE_DIR}/include)
LINK_DIRECTORIES(${FAR_SOURCE_DIR}/lib)

# Execute the FindTBB macro
# If you want FindTBB to fail if TBB is not found, use this form:
#find_package(TBB REQUIRED)

#if(NOT TBB_FOUND)
 #       MESSAGE (STATUS "TBB not found. Example project will not be built.")
#else(NOT TBB_FOUND)

ADD_EXECUTABLE(far seqan_trimmer.cpp)
TARGET_LINK_LIBRARIES(far tbb)

ADD_EXECUTABLE(PairedreadFinder PairedreadFinder.cpp)
TARGET_LINK_LIBRARIES(PairedreadFinder tbb)
