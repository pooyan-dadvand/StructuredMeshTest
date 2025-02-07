@echo off
call "C:\Program Files (x86)\Microsoft Visual Studio\2019\Professional\VC\Auxiliary\Build\vcvars64.bat"
cl /Ox /EHsc /std:c++17 /Zi main.cpp /Fe:test.exe
echo C++ Compilation complete.
nvcc -o test_cuda.exe test_cuda.cu
echo CUDA Compilation complete.