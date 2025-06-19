cmake_minimum_required(VERSION 3.10)
project(kNN_Index)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add compiler flags from README
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Compile predata.cpp to pre
add_executable(pre predata.cpp)

# Compile index_bng.cpp to bng
# Note: The README shows an execution of ./bng, but no compile command.
# Assuming it's built from index_bng.cpp, similar to the others.
add_executable(bng index_bng.cpp)

# Compile query_up.cpp to qu
add_executable(qu query_up.cpp)

