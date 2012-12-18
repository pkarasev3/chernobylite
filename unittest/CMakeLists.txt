
# Define executable: test simple tracking on webcam
add_executable(test_opencv_planar_H test_opencv_planar_H.cpp)
target_link_libraries(test_opencv_planar_H  chernobylite
                      ${Boost_LIBRARIES}  ${OpenCV_LIBS})


# test app for "concurrent access to data structure"
add_executable(google_prep_concurrency  google_prep_concurrency.cpp)
target_link_libraries(google_prep_concurrency
                      ${Boost_LIBRARIES}
                      ${Boost_SYSTEM_LIBRARY}
                      ${Boost_THREAD_LIBRARY} )

