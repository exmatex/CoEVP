
The subdirectories of this directory contain the files needed to build the CoEVP proxy
application described in the specification document CoEVP.pdf.  The constitutive models
are contained in a library (cmlib.a) which is statically linked with LULESH.  To build
the library, link with LULESH and run, beginning in the directory containing this README
file:

     make 

To build a minimal set of features one can use:

     make REDIS=no FLANN=no SILO=no

which will enable building of the redis backend, the libflann based neighbor search and
silo output files. 

To use a system redis, libflann or libsilo use:

     make REDIS_LOC=/path/to/hiredis FLANN_LOC=/path/to/flann SILO_LOC=/path/to/silo

To build an mpi version of lulesh use:

     make COEVP_MPI=yes

To run lulesh use:

     LULESH/luslesh

and see use the "--help" option for help.
