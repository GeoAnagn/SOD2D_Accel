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
/opt/nvidia/hpc_sdk/Linux_x86_64/22.9/comm_libs/hpcx/latest/ompi/bin/mpif90 -DNOPRED -cpp -lstdc++ -D_USE_NVTX -lnvToolsExt -cuda -Minfo=accel ${flags} -c elem_convec_solo_og.f90 -o elem_convec_og.o
/opt/nvidia/hpc_sdk/Linux_x86_64/22.9/comm_libs/hpcx/latest/ompi/bin/mpif90 -DNOPRED -cpp -lstdc++ -D_USE_NVTX -lnvToolsExt -cuda -Minfo=accel ${flags} -c elem_diffu_solo_og.f90 -o elem_diffu_og.o

/opt/nvidia/hpc_sdk/Linux_x86_64/22.9/comm_libs/hpcx/hpcx-2.12/ompi/bin/mpif90 -DNOPRED -cpp -lstdc++ -D_USE_NVTX -lnvToolsExt -cuda -Minfo=accel ${flags} ./mod_constants.f90.o ./mod_nvtx.f90.o ./elem_convec_og.o -o run_elem_convec -Wl, -rpath, /usr/lib/x86_64-linux-gnu/libpthread.so /usr/lib/x86_64-linux-gnu/libz.so /usr/lib/x86_64-linux-gnu/libdl.so -lm /usr/lib/x86_64-linux-gnu/libpthread.so /usr/lib/x86_64-linux-gnu/libz.so /usr/lib/x86_64-linux-gnu/libdl.so -lm
/opt/nvidia/hpc_sdk/Linux_x86_64/22.9/comm_libs/hpcx/hpcx-2.12/ompi/bin/mpif90 -DNOPRED -cpp -lstdc++ -D_USE_NVTX -lnvToolsExt -cuda -Minfo=accel ${flags} ./mod_constants.f90.o ./mod_nvtx.f90.o ./elem_diffu_og.o -o run_elem_diffu -Wl, -rpath, /usr/lib/x86_64-linux-gnu/libpthread.so /usr/lib/x86_64-linux-gnu/libz.so /usr/lib/x86_64-linux-gnu/libdl.so -lm /usr/lib/x86_64-linux-gnu/libpthread.so /usr/lib/x86_64-linux-gnu/libz.so /usr/lib/x86_64-linux-gnu/libdl.so -lm