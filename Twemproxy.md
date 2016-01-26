To utilize distributed redis, the following steps must be taken

1. Ensure that Twemproxy support is built in by running make with TWEMPROXY=yes
2. Prior to running LULESH, run twemproxy/buildNutcracker.py from the directory in which you will run LULESH.
3. Run lulesh with -s, -r, and -R flags enabled

If all goes well, each node will spawn a redis instance and a nutcracker (twemproxy) instance and the databases will be linked.

NOTE: Currently twemproxy support has very limited testing and is designed more with the charm++ version in mind than the MPI version. However, limited testing suggests that MPI will work so long as there is one redis/nutcracker instance per physical node.
