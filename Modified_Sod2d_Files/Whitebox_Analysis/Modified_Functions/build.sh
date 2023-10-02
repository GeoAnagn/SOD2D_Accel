#!/bin/bash
set -x
target=$1
#0 CPU, 1 GPU

if [ $target -eq 0 ]
then
	flags="-DNOACC -fast -DNOACC"
else
	flags="-gpu=cuda11.7,managed,lineinfo -cuda -acc -fast"
fi

/opt/nvidia/hpc_sdk/Linux_x86_64/22.9/comm_libs/hpcx/latest/ompi/bin/mpif90 -DNOPRED -cpp -lstdc++ -D_USE_NVTX -lnvToolsExt -cuda -Minfo=accel ${flags} -c /home/apps/sod2d/src/lib_sod2d/sources/mod_constants.f90 -o ./mod_constants.f90.o
/opt/nvidia/hpc_sdk/Linux_x86_64/22.9/comm_libs/hpcx/latest/ompi/bin/mpif90 -DNOPRED -cpp -lstdc++ -D_USE_NVTX -lnvToolsExt -cuda -Minfo=accel ${flags} -c /home/apps/sod2d/src/lib_sod2d/sources/mod_nvtx.f90 -o ./mod_nvtx.f90.o
/opt/nvidia/hpc_sdk/Linux_x86_64/22.9/comm_libs/hpcx/latest/ompi/bin/mpif90 -DNOPRED -cpp -lstdc++ -D_USE_NVTX -lnvToolsExt -cuda -Minfo=accel ${flags} -c full_convec_ijk.f90 -o full_convec_ijk.o
/opt/nvidia/hpc_sdk/Linux_x86_64/22.9/comm_libs/hpcx/latest/ompi/bin/mpif90 -DNOPRED -cpp -lstdc++ -D_USE_NVTX -lnvToolsExt -cuda -Minfo=accel ${flags} -c full_diffu_ijk.f90 -o full_diffu_ijk.o

/opt/nvidia/hpc_sdk/Linux_x86_64/22.9/comm_libs/hpcx/hpcx-2.12/ompi/bin/mpif90 -DNOPRED -cpp -lstdc++ -D_USE_NVTX -lnvToolsExt -cuda -Minfo=accel ${flags} ./mod_constants.f90.o ./mod_nvtx.f90.o ./full_convec_ijk.o -o full_convec_ijk -Wl, -rpath, /usr/lib/x86_64-linux-gnu/libpthread.so /usr/lib/x86_64-linux-gnu/libz.so /usr/lib/x86_64-linux-gnu/libdl.so -lm /usr/lib/x86_64-linux-gnu/libpthread.so /usr/lib/x86_64-linux-gnu/libz.so /usr/lib/x86_64-linux-gnu/libdl.so -lm
/opt/nvidia/hpc_sdk/Linux_x86_64/22.9/comm_libs/hpcx/hpcx-2.12/ompi/bin/mpif90 -DNOPRED -cpp -lstdc++ -D_USE_NVTX -lnvToolsExt -cuda -Minfo=accel ${flags} ./mod_constants.f90.o ./mod_nvtx.f90.o ./full_diffu_ijk.o -o full_diffu_ijk -Wl, -rpath, /usr/lib/x86_64-linux-gnu/libpthread.so /usr/lib/x86_64-linux-gnu/libz.so /usr/lib/x86_64-linux-gnu/libdl.so -lm /usr/lib/x86_64-linux-gnu/libpthread.so /usr/lib/x86_64-linux-gnu/libz.so /usr/lib/x86_64-linux-gnu/libdl.so -lm