# Rename original SOD2D files to *_old 
mv /home/apps/sod2d/src/lib_sod2d/sources/elem_convec.f90 /home/apps/sod2d/src/lib_sod2d/sources/elem_convec_old.f90
mv /home/apps/sod2d/src/lib_sod2d/sources/elem_diffu.f90 /home/apps/sod2d/src/lib_sod2d/sources/elem_diffu_old.f90
mv /home/apps/sod2d/src/lib_sod2d/sources/mod_analysis.f90 /home/apps/sod2d/src/lib_sod2d/sources/mod_analysis_old.f90
mv /home/apps/sod2d/src/lib_sod2d/sources/mod_entropy_viscosity.f90 /home/apps/sod2d/src/lib_sod2d/sources/mod_entropy_viscosity_old.f90
mv /home/apps/sod2d/src/lib_mainBaseClass/sources/TGVSolver.f90 /home/apps/sod2d/src/lib_mainBaseClass/sources/TGVSolver_old.f90

# Copy mofified SOD2D files to appropriate folders
cp /home/apps/Refmap/SOD2D_Accel/Modified_Sod2d_Files/Blackbox_Analysis/elem_convec.f90 /home/apps/sod2d/src/lib_sod2d/sources/
cp /home/apps/Refmap/SOD2D_Accel/Modified_Sod2d_Files/Blackbox_Analysis/elem_diffu.f90 /home/apps/sod2d/src/lib_sod2d/sources/
cp /home/apps/Refmap/SOD2D_Accel/Modified_Sod2d_Files/Blackbox_Analysis/mod_analysis.f90 /home/apps/sod2d/src/lib_sod2d/sources/
cp /home/apps/Refmap/SOD2D_Accel/Modified_Sod2d_Files/Blackbox_Analysis/mod_entropy_viscosity.f90 /home/apps/sod2d/src/lib_sod2d/sources/
cp /home/apps/Refmap/SOD2D_Accel/Modified_Sod2d_Files/Blackbox_Analysis/TGVSolver.f90 /home/apps/sod2d/src/lib_mainBaseClass/sources/

# Rebuild SOD2D
cd /home/apps/sod2d/build
make clean
make
cd /home/apps/Refmap/SOD2D_Accel