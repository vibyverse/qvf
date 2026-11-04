# License: Apache-2.0.
# Cross toolchain for a big-endian check of the C++ library on s390x Linux,
# run under QEMU user-mode emulation. On Debian/Ubuntu:
#   apt-get install g++-s390x-linux-gnu qemu-user
# CTest runs every test through the emulator (CMAKE_CROSSCOMPILING_EMULATOR).
set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_SYSTEM_PROCESSOR s390x)

set(QVF_CROSS_TRIPLE s390x-linux-gnu)
set(CMAKE_C_COMPILER ${QVF_CROSS_TRIPLE}-gcc)
set(CMAKE_CXX_COMPILER ${QVF_CROSS_TRIPLE}-g++)
set(CMAKE_Fortran_COMPILER ${QVF_CROSS_TRIPLE}-gfortran)

set(CMAKE_FIND_ROOT_PATH /usr/${QVF_CROSS_TRIPLE})
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)

find_program(QVF_QEMU NAMES qemu-s390x qemu-s390x-static REQUIRED)
set(CMAKE_CROSSCOMPILING_EMULATOR ${QVF_QEMU} -L /usr/${QVF_CROSS_TRIPLE})
