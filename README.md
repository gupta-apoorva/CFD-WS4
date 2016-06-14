To run the code please use 

mpirun -np num_proc ./sim problem.dat

here num_proc is for the number of processors and its value is equal to iproc*jproc i.e. the number of parts in which you want to break the system.

for eg if iproc is 2 and jproc is 1 then num_proc is 2.
