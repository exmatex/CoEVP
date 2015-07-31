
The subdirectories of this directory contain the files needed to build the CoEVP proxy
application described in the specification document CoEVP.pdf.  The constitutive models
are contained in a library (cmlib.a) which is statically linked with LULESH.  To build
the library, link with LULESH and run, beginning in the directory containing this README
file:

     cd CM/exec

     gmake (builds and installs cmlib.a and header files needed by LULESH)

     cd ../../LULESH

     gmake   (builds LULESH)

     lulesh  (runs LULESH on the Taylor cylinder problem with the specified
     constitutive model)

Additional instructions for running with adaptive sampling and the redis
database.  To enable adaptive sampling you must edit `lulesh.cc` and define
`USE_ADAPTIVE_SAMPLING`. There is also a slight change to the build instructions:

``` sh
cd CM/exec
make REDIS=yes -j10
cd ../../LULESH
make clean;make
./lulesh
```
Before running LULESH, make sure and start a redis-server (instructions for
doing this depend on your system environment). CoEVP will connect to redis on
the standard address (`localhost:6379`). If CoEVP detects any data in the redis databse
(from a previous run) it will be deleted. This is a (temporary) convenience so
that you don't have to find and delete the `dump.rdb` file. At  the end of the
run, CoEVP will print a summary of database statistics.


