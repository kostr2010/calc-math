project(lab-9)

add_executable(main-9 main.cpp)
add_library(lib-9 lib.cpp)
target_link_libraries(main-9 lib-9)

