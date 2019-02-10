@echo off
set MATLAB=C:\PROGRA~1\MATLAB\R2018a
set MATLAB_ARCH=win64
set MATLAB_BIN="C:\Program Files\MATLAB\R2018a\bin"
set ENTRYPOINT=mexFunction
set OUTDIR=.\
set LIB_NAME=minSum_mex
set MEX_NAME=minSum_mex
set MEX_EXT=.mexw64
call setEnv.bat
echo # Make settings for minSum > minSum_mex.mki
echo CC=%COMPILER%>> minSum_mex.mki
echo CXX=%CXXCOMPILER%>> minSum_mex.mki
echo CFLAGS=%COMPFLAGS%>> minSum_mex.mki
echo CXXFLAGS=%CXXCOMPFLAGS%>> minSum_mex.mki
echo LINKER=%LINKER%>> minSum_mex.mki
echo OPTIMFLAGS=%OPTIMFLAGS%>> minSum_mex.mki
echo DEBUGFLAGS=%DEBUGFLAGS%>> minSum_mex.mki
echo LINKFLAGS=%LINKFLAGS%>> minSum_mex.mki
echo LINKOPTIMFLAGS=%LINKOPTIMFLAGS%>> minSum_mex.mki
echo LINKDEBUGFLAGS=%LINKDEBUGFLAGS%>> minSum_mex.mki
echo MATLAB_ARCH=%MATLAB_ARCH%>> minSum_mex.mki
echo OMPFLAGS= >> minSum_mex.mki
echo OMPLINKFLAGS= >> minSum_mex.mki
echo EMC_COMPILER=mingw64>> minSum_mex.mki
echo EMC_CONFIG=optim>> minSum_mex.mki
"C:\Program Files\MATLAB\R2018a\bin\win64\gmake" -j 1 -B -f minSum_mex.mk
