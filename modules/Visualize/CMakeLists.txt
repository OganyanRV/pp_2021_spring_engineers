project("Filter visualsaion")
find_package(OpenCV REQUIRED)
add_executable(Visual visualize.cpp)
target_link_libraries(Visual ${OpenCV_LIBS})
add_library(mark_comp ../task_1/oganyan_r_mark_components/mark_components.h ../task_1/oganyan_r_mark_components/mark_components.cpp)
target_link_libraries(Visual mark_comp)

message( STATUS "-- " Visualisation )