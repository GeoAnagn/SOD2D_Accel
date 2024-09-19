
cp -r ../SOD2D_Versions/sod2d_gitlab ../SOD2D_Versions/sod2d_profiling
cp -r ../Modified_Sod2d_Files/Nsight_Modified_Files/lib_sod2d/sources/ ../SOD2D_Versions/sod2d_profiling/src/lib_sod2d
cp -r ../Modified_Sod2d_Files/Nsight_Modified_Files/lib_mainBaseClass/sources/ ../SOD2D_Versions/sod2d_profiling/src/lib_mainBaseClass
cp -r ../Modified_Sod2d_Files/Nsight_Modified_Files/lib_sod2d_incomp/sources/ ../SOD2D_Versions/sod2d_profiling/src/lib_sod2d_incomp
cp -r ../Modified_Sod2d_Files/Nsight_Modified_Files/tool_meshConversorPar/CMakeLists.txt ../SOD2D_Versions/sod2d_profiling/tool_meshConversorPar/CMakeLists.txt
cp -r ../Modified_Sod2d_Files/Nsight_Modified_Files/tool_commsPerformance/CMakeLists.txt ../SOD2D_Versions/sod2d_profiling/tool_commsPerformance/CMakeLists.txt
cd ../SOD2D_Versions/sod2d_profiling
mkdir build && cd build
cmake -DUSE_GPU=ON -DUSE_MEM_MANAGED=ON ..
make -j 8
cd ../../../Scripts

cp -r ../SOD2D_Versions/sod2d_gitlab ../SOD2D_Versions/sod2d_blackbox_coupled
cp -r ../Modified_Sod2d_Files/Blackbox_Analysis/Coupled/TGV/lib_sod2d/sources ../SOD2D_Versions/sod2d_blackbox_coupled/src/lib_sod2d
cd ../SOD2D_Versions/sod2d_blackbox_coupled
mkdir build && cd build
cmake -DUSE_GPU=ON -DUSE_MEM_MANAGED=ON ..
make -j 8
cd ../../../Scripts

cp -r ../SOD2D_Versions/sod2d_gitlab ../SOD2D_Versions/sod2d_blackbox_decoupled
cp -r ../Modified_Sod2d_Files/Blackbox_Analysis/Decoupled/TGV/lib_sod2d/sources ../SOD2D_Versions/sod2d_blackbox_decoupled/src/lib_sod2d
cd ../SOD2D_Versions/sod2d_blackbox_decoupled
mkdir build && cd build
cmake -DUSE_GPU=ON -DUSE_MEM_MANAGED=ON ..
make -j 8
cd ../../../Scripts

cp -r ../SOD2D_Versions/sod2d_gitlab ../SOD2D_Versions/sod2d_cpu
rm ../SOD2D_Versions/sod2d_cpu/cmake/compilerOps.cmake
cp ../Modified_Sod2d_Files/Cpu_Files/cmake/compilerOps.cmake ../SOD2D_Versions/sod2d_cpu/cmake/compilerOps.cmake
cd ../SOD2D_Versions/sod2d_cpu
mkdir build && cd build
cmake -DUSE_GPU=OFF ..
make -j 8
cd ../../../Scripts