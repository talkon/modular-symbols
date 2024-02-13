# usage:
#   ./build.sh [d|D] [c] [T]
#   release build by default
#   d: debug build (-O3)
#   D: debug build (-O0)
#   c: clean
#   T: run tests

if [[ $* == *D* ]]; then
  : ${CMAKE_DIR:=cmake-build-debug}
  cmake -D CMAKE_BUILD_TYPE=Debug -S . -B $CMAKE_DIR
elif [[ $* == *d* ]]; then
  : ${CMAKE_DIR:=cmake-build-relwithdebinfo}
  cmake -D CMAKE_BUILD_TYPE=RelWithDebInfo -S . -B $CMAKE_DIR
else
  : ${CMAKE_DIR:=cmake-build-release}
  cmake -D CMAKE_BUILD_TYPE=Release -S . -B $CMAKE_DIR
fi

if [[ $* == *c* ]]; then
  cmake --build $CMAKE_DIR --target clean
else
  cmake --build $CMAKE_DIR
fi

if [[ $* == *T* ]]; then
  ctest --test-dir $CMAKE_DIR
fi