The Fortran wrappers now return an optional `message` with exactly the copied
error text on gfortran 11 and 12, which previously returned empty or
garbage-padded text. `qvf_gto` refuses NaN exponents without raising an IEEE
invalid exception, so callers with floating-point trapping enabled get a
refusal instead of a signal. Focused development checks now cover gfortran
11, 12, 13 and 14 and LLVM flang 18, 19 and 20 on Ubuntu 24.04 x86_64. The C
ABI and the archive format are unchanged.
