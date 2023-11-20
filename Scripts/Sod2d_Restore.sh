# Remove modified SOD2D files
rm /home/apps/sod2d/src/lib_sod2d/sources/elem_convec.f90
rm /home/apps/sod2d/src/lib_sod2d/sources/elem_diffu.f90
rm /home/apps/sod2d/src/lib_sod2d/sources/mod_analysis.f90
rm /home/apps/sod2d/src/lib_sod2d/sources/mod_entropy_viscosity.f90
rm /home/apps/sod2d/src/lib_mainBaseClass/sources/TGVSolver.f90

# Rename original SOD2D files to the original names (remove *_old) 
mv /home/apps/sod2d/src/lib_sod2d/sources/elem_convec_old.f90 /home/apps/sod2d/src/lib_sod2d/sources/elem_convec.f90
mv /home/apps/sod2d/src/lib_sod2d/sources/elem_diffu_old.f90 /home/apps/sod2d/src/lib_sod2d/sources/elem_diffu.f90
mv /home/apps/sod2d/src/lib_sod2d/sources/mod_analysis_old.f90 /home/apps/sod2d/src/lib_sod2d/sources/mod_analysis.f90
mv /home/apps/sod2d/src/lib_sod2d/sources/mod_entropy_viscosity_old.f90 /home/apps/sod2d/src/lib_sod2d/sources/mod_entropy_viscosity.f90
mv /home/apps/sod2d/src/lib_mainBaseClass/sources/TGVSolver_old.f90 /home/apps/sod2d/src/lib_mainBaseClass/sources/TGVSolver.f90

# Rebuild SOD2D
cd /home/apps/sod2d/build
make clean
make