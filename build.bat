@echo off
set build_dir=build
set bin_dir=bin

if not exist %build_dir% (
    mkdir %build_dir%
)

if not exist %bin_dir% (
    mkdir %bin_dir%
)

cd %build_dir%
call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Professional\VC\Auxiliary\Build\vcvars64.bat"
cl /Ox /EHsc /std:c++17 /Zi ..\main.cpp /Fe:..\%bin_dir%\test.exe
echo C++ Compilation complete.
nvcc -o ..\%bin_dir%\test_cuda.exe ..\test_cuda.cu
echo CUDA Compilation complete.