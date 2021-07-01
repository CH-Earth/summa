

Extract ida-5.4.0.tar.gz and read INSTALL_GUIDE.pdf

Notes:

-In SUMMA we use the Fortran/C interface (fida)
You need to edit CMakeLists.txt a bit.
Search the keyword "Fortan" in CMakeLists.txt . Enable the appropriate options. Note that we are using Fortran 2003
The most important one is:
   # Fortran 2003 interface is disabled by default
   set(DOCSTR "Enable Fortran 2003 interfaces")
   option(F2003_INTERFACE_ENABLE "${DOCSTR}" OFF)
You should enable the above option

-The default compiler is gfortan. To change it you can add the following optoin to camke
   -DCMAKE_Fortran_COMPILER=your_favorite_compiler
   
-Make sure SUMMA can access to installed modules of sundials 

