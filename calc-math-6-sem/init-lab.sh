is_number='^[0-9]+$'

dir=lab-"$1"

readme_md="""# DESCRIPTION
how to build & run:

\`./lab-$1/main-$1\`

\`res\` folder contains an example of the result of this lab

\`task\` folder contains the task of this lab

# REQUIREMENTS
c/c++: cmake
"""

cmake_lists_txt="""project(${dir})

add_executable(main-$1 main.cpp)
add_library(lib-$1 lib.cpp)
target_link_libraries(main-$1 lib-$1)
"""

lib_h="""#ifndef LIB_H_INCLUDED
#define LIB_H_INCLUDED

#include <cstddef>

#endif
"""

lib_inl_h="""#ifndef LIB_INL_H_INCLUDED
#define LIB_INL_H_INCLUDED

#include \"lib.h\"

#endif
"""

lib_cpp="""#include \"lib-inl.h\"
"""

main_cpp="""#include \"lib-inl.h\"

int main(int argc, char* argv[]) {

    return 0;
}
"""

if ! [[ $1 =~ $is_number ]] ; then
    echo "> error, expected lab's number" >&2; exit 1
elif mkdir $dir ; then
    mkdir $dir/res
    mkdir $dir/task
    touch $dir/main.cpp
    touch $dir/lib.h
    touch $dir/lib-inl.h
    touch $dir/lib.cpp
    touch $dir/README.md
    touch $dir/CMakeLists.txt
    touch $dir/res/README.md

    echo "${cmake_lists_txt}">> $dir/CMakeLists.txt
    echo "${readme_md}" >> $dir/README.md
    echo "${lib_h}" >> $dir/lib.h
    echo "${lib_inl_h}" >> $dir/lib-inl.h
    echo "${lib_cpp}" >> $dir/lib.cpp
    echo "${main_cpp}" >> $dir/main.cpp
else
    echo "> error, ${dir} already exists" >&2; exit 2
fi