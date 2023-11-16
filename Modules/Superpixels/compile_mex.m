mex -setup C++ -v

%% make SLIC
mex ./SLIC/fuzzyslic_mex.c -outdir ./
mex ./SLIC/slic_mex.c -outdir ./
mex ./SLIC/slico_mex.c -outdir ./
mex ./SLIC/slicsupervoxel_mex.c -outdir ./

%% make SNIC
mex ./SNIC/snic_mex.cpp -outdir ./

%% make ERS
mex -v -c ERS/MERCCInput.cpp
mex -v -c ERS/MERCOutput.cpp
mex -v -c ERS/MERCDisjointSet.cpp
mex -v -c ERS/MERCFunctions.cpp
mex -v -c ERS/MERCLazyGreedy.cpp
mex ERS/ers_mex.cpp MERCCInput.obj MERCOutput.obj MERCDisjointSet.obj MERCFunctions.obj MERCLazyGreedy.obj
delete *.obj

%% make scalp
mex -O CFLAGS="\$CFLAGS -Wall -Wextra -W -std=c99" ./SCALP/scalp_mex.cpp -outdir ./

%% make gGMMSP
copyfile ./gGMMSP/cudart64_80.dll ./
copyfile ./gGMMSP/mx_gGMMSP.mexa64 ./
copyfile ./gGMMSP/mx_gGMMSP.mexw64 ./

%% make gGMMSP
copyfile ./GMMSP/cudart64_80.dll ./
copyfile ./GMMSP/mx_GMMSP.mexw64 ./
