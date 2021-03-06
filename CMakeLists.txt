project(PARTICLE_FILTER)

cmake_minimum_required (VERSION 3.5)

add_definitions(-std=c++11)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")

set(sources src/main.cpp src/particle.cpp src/particle.h 
src/landmark.cpp src/landmark.h
src/particle_filter.cpp src/particle_filter.h src/helper_functions.h ) 

# gtest 
#set(test_sources test/gtest_main.cpp src/ukf.cpp src/ukf.h src/tools.cpp src/sensor_data.cpp #src/sensor_data.h test/ukf_test.cpp test/util.cpp test/util.h)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 

include_directories(/usr/local/include)
include_directories(/usr/local/opt/openssl/include)
include_directories(src)
link_directories(/usr/local/lib)
link_directories(/usr/local/opt/openssl/lib)
link_directories(/usr/local/Cellar/libuv/1.19.2/lib)

endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 


add_executable(particle_filter ${sources})

# gtest 
#add_executable(Test_UnscentedKF ${test_sources})
#set_target_properties(Test_UnscentedKF PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(particle_filter z ssl uv uWS)

# gtest 
#target_link_libraries(Test_UnscentedKF gtest gmock z ssl uv uWS)