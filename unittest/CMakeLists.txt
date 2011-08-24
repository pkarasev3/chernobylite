# Define executable: test simple tracking on webcam
add_executable(test_opencv_planar_H test_opencv_planar_H.cpp)
  target_link_libraries(test_opencv_planar_H  chernobylite ${Boost_LIBRARIES}  ${OpenCV_LIBS})

