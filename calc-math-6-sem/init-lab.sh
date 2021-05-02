is_number='^[0-9]+$'
dir=lab-"$1"

if ! [[ $1 =~ $is_number ]] ; then
    echo "> error, expected lab's number" >&2; exit 1
elif mkdir $dir ; then
    mkdir $dir/res
    touch $dir/main.cpp
    touch $dir/lib.h
    touch $dir/lib-inl.h
    touch $dir/lib.cpp
    touch $dir/README.md
    touch $dir/CMakeLists.txt
    echo "project(${dir})" >> $dir/CMakeLists.txt
    echo "add_executable(main main.cpp)" >> $dir/CMakeLists.txt
    echo "add_library(lib lib.cpp)" >> $dir/CMakeLists.txt
    echo "target_link_libraries(main lib)" >> $dir/CMakeLists.txt
    echo "# ${dir}" >> $dir/README.md
else
    echo "> error, ${dir} already exists" >&2; exit 2
fi