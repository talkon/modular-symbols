# usage:
#   ./build.sh [d|D] [c] [T]
#   release build by default
#   d: debug build (-O3)
#   D: debug build (-O0)
#   c: clean
#   C: delete build repo
#   T: run tests

if [[ $* == *D* ]]; then
  : ${CMAKE_DIR:=cmake-build-debug}
  : ${CMAKE_BUILD_TYPE:=Debug}
elif [[ $* == *d* ]]; then
  : ${CMAKE_DIR:=cmake-build-relwithdebinfo}
  : ${CMAKE_BUILD_TYPE:=RelWithDebInfo}
else
  : ${CMAKE_DIR:=cmake-build-release}
  : ${CMAKE_BUILD_TYPE:=Release}
fi

if [[ $* == *C* ]]; then
  rm -rf $CMAKE_DIR
elif [[ $* == *c* ]]; then
  cmake -D CMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE -S . -B $CMAKE_DIR
  cmake --build $CMAKE_DIR --target clean
else
  cmake -D CMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE -S . -B $CMAKE_DIR
  cmake --build $CMAKE_DIR
fi

if [[ $* == *T* ]]; then
  ctest --test-dir $CMAKE_DIR
fi