project(lab-8)
add_executable(main-8 main.cpp)
add_library(lib-8 lib.cpp)
target_link_libraries(main-8 lib-8)