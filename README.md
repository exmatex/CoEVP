
The subdirectories of this directory contain the files needed to build the CoEVP proxy
application described in the specification document CoEVP.pdf.  The constitutive models
are contained in a library (cmlib.a) which is statically linked with LULESH.  To build
the library, link with LULESH and run, beginning in the directory containing this README
file:

     make 

To build a minimal set of features one can use:

     make REDIS=no FLANN=no SILO-no

which will enable building of the redis backend, the libflann based neighbor search and
silo output files.

