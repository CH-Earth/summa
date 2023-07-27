# from {$maindir}/sundials/builddir, run 
# cp ../../summa/build/makefiles/Sundials/build_cmakeSundials.sh build_cmake
# run script from the builddir directory with ./build_cmake
# Note, using -DEXAMPLES_ENABLE_C=OFF -DEXAMPLES_ENABLE_F2003=OFF, if want to run examples should change

cmake ../../sundials-6.3.0/ -DEXAMPLES_ENABLE_C=OFF -DEXAMPLES_ENABLE_F2003=OFF -DBUILD_FORTRAN_MODULE_INTERFACE=ON -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_INSTALL_PREFIX=../../sundials/instdir -DEXAMPLES_INSTALL_PATH=../../sundials/instdir/examples


