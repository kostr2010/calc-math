project(lab-10)

add_executable(main-10 main.cpp)
add_library(lib-10 lib.cpp)
target_link_libraries(main-10 lib-10)

