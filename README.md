RABBIT Installation Instructions

1. Request access to the RABBIT source code on Gitlab
2. Clone the RABBIT repository into the directory of your choice using: *git clone git@gitlab.mpcdf.mpg.de:markusw/rabbit.git* 
3. Within the newly-created rabbit directory, create a new directory called “build”: *mkdir build* 
4. Navigate into the build directory: *cd build*
5. Generate cmake files, setting flags appropriately (details below): *cmake .. -DCMAKE_BUILD_TYPE=release -DCMAKE_Fortran_COMPILER=gfortran -DNETCDF_HOME=/opt/local -DOpenMP_Fortran_FLAGS=-fopenmp* 
   - DCMAKE_Fortran_COMPILER instructs cmake about which Fortran compiler it should use to compile RABBIT. Most likely, this should be set to gfortran. You can check if you have gfortran installed by running *gfortran --version*. If it returns something like this: GNU Fortran (Homebrew GCC 13.2.0) 13.2.0 , then set DCMAKE_Fortran_COMPILER=gfortran. 
   - DNETCDF_HOME tells cmake where NetCDF and its associated libraries can be found on your machine. If you have previously installed NetCDF using MacPorts, then this should probably be set to DNETCDF_HOME=/opt/local. If instead you have installed it with Homebrew, then it should probably be set to DNETCDF_HOME=/opt/homebrew. 
      - In case you don’t have NetCDF installed at all, you can use Homebrew to install it by running brew install netcdf, and then set DNETCDF_HOME=/opt/homebrew 
6. Once cmake files have been successfully generated, run make from inside the build directory: *make*

Troubleshooting: 
- On step 5, I run into this error: CMake Error at CMakeLists.txt:62 (message):  NETCDF_HOME variable seems wrong (set with -DNETCDF_HOME)
   - You may have both netcdf and netcdf-fortran installed, and cmake therefore won’t know which one you want to link to. You can distinguish between them by adding an extra flag (-DNETCDFF_HOME) to the cmake command from step 5. For example, cmake .. -DCMAKE_BUILD_TYPE=release -DCMAKE_Fortran_COMPILER=gfortran -DNETCDF_HOME=/opt/homebrew/Cellar/netcdf/4.9.2_1 -DNETCDFF_HOME=/opt/homebrew/Cellar/netcdf-fortran/4.6.1 -DOpenMP_Fortran_FLAGS=-fopenmp 
      - Navigate to your /opt/homebrew/Cellar in order to check which version numbers of each are appropriate for your flags 
 - On step 5, I run into this error: CMake Error at CMakeLists.txt:27 (project): The CMAKE_Fortran_COMPILER: ifort is not a full path and was not found in the PATH. 
   - You have either pointed cmake to a compiler that doesn’t exist on your machine or isn’t in your path. If you’re sure that you have the compiler that you’ve chosen installed, then add its location to your $PATH in .zshrc or .bashrc. If you’re not sure whether you have that compiler, run, for example, *ifort --version*. If it returns something like “command not found”, it means you don’t have the selected compiler and you’ll have to choose a different one. 

