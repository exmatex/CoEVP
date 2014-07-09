
The subdirectories of this directory contain the files needed to build the CoEVP proxy
application described in the specification document CoEVP.pdf.  The constitutive models
are contained in a library (cmlib.a) which is statically linked with LULESH.  To build
the library, link with LULESH and run, beginning in the directory containing this README
file:

     cd CM/exec

     gmake (builds and installs cmlib.a and header files needed by LULESH)

     cd ../../LULESH

     gmake   (builds LULESH)

     lulesh  (runs LULESH on the Taylor cylinder problem with the specified constitutive model)

