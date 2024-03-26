cmake_minimum_required(VERSION 3.0.2)
project(rosneuro_pwelch)


find_package(catkin REQUIRED COMPONENTS
  roscpp
  rosneuro_windows
  rosneuro_filters
  rosneuro_windows_flattop
  rosneuro_windows_blackman
  rosneuro_windows_hann
  rosneuro_windows_hamming
  rosneuro_buffers_ringbuffer
  Eigen3
)

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES 
  ${PROJECT_NAME}
	
  CATKIN_DEPENDS
  	roscpp
	  roslib
	  std_msgs
	  rosneuro_filters
    rosneuro_windows_flattop
    rosneuro_windows_blackman
    rosneuro_windows_hann
    rosneuro_windows_hamming
    rosneuro_windows
    rosneuro_buffers_ringbuffer
  DEPENDS
)
include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
)

catkin_add_gtest(TestPwelch test/TestPwelch.cpp)
target_link_libraries(TestPwelch ${GTEST_BOTH_LIBRARIES} pthread gmock fftw3)
target_link_libraries(TestPwelch ${PROJECT_NAME} ${catkin_LIBRARIES})

catkin_add_gtest(TestPwelchNode test/TestPwelchNode.cpp)
target_link_libraries(TestPwelchNode ${GTEST_BOTH_LIBRARIES} pthread gmock fftw3)
target_link_libraries(TestPwelchNode ${PROJECT_NAME} ${catkin_LIBRARIES})

add_library(${PROJECT_NAME} src/Pwelch.cpp)
target_link_libraries(${PROJECT_NAME} ${Eigen3_LIBRARIES} 
									  ${catkin_LIBRARIES}
)


add_executable(test_pwelch test/test_pwelch.cpp)
target_link_libraries(test_pwelch ${catkin_LIBRARIES} fftw3)

add_executable(pwelch_node src/pwelch_node.cpp)
target_link_libraries(pwelch_node ${catkin_LIBRARIES} fftw3)
