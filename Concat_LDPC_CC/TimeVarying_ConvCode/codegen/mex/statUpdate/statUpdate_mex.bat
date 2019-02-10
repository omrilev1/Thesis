@echo off
set MATLAB=C:\PROGRA~1\MATLAB\R2018a
set MATLAB_ARCH=win64
set MATLAB_BIN="C:\Program Files\MATLAB\R2018a\bin"
set ENTRYPOINT=mexFunction
set OUTDIR=.\
set LIB_NAME=statUpdate_mex
set MEX_NAME=statUpdate_mex
set MEX_EXT=.mexw64
call setEnv.bat
echo # Make settings for statUpdate > statUpdate_mex.mki
echo CC=%COMPILER%>> statUpdate_mex.mki
echo CXX=%CXXCOMPILER%>> statUpdate_mex.mki
echo CFLAGS=%COMPFLAGS%>> statUpdate_mex.mki
echo CXXFLAGS=%CXXCOMPFLAGS%>> statUpdate_mex.mki
echo LINKER=%LINKER%>> statUpdate_mex.mki
echo OPTIMFLAGS=%OPTIMFLAGS%>> statUpdate_mex.mki
echo DEBUGFLAGS=%DEBUGFLAGS%>> statUpdate_mex.mki
echo LINKFLAGS=%LINKFLAGS%>> statUpdate_mex.mki
echo LINKOPTIMFLAGS=%LINKOPTIMFLAGS%>> statUpdate_mex.mki
echo LINKDEBUGFLAGS=%LINKDEBUGFLAGS%>> statUpdate_mex.mki
echo MATLAB_ARCH=%MATLAB_ARCH%>> statUpdate_mex.mki
echo OMPFLAGS= >> statUpdate_mex.mki
echo OMPLINKFLAGS= >> statUpdate_mex.mki
echo EMC_COMPILER=mingw64>> statUpdate_mex.mki
echo EMC_CONFIG=optim>> statUpdate_mex.mki
"C:\Program Files\MATLAB\R2018a\bin\win64\gmake" -j 1 -B -f statUpdate_mex.mk
