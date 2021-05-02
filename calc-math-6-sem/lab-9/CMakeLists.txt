project(lab-9)

add_executable(main main.cpp)
add_library(lib lib.cpp)
target_link_libraries(main lib)