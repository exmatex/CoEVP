Using the Logging Library
---

This library, when linked in with CoEVP (and eventually Tabasco), will provide
limited facilities for logging performance and execution data during a run of
the code. The logging library currently uses REDIS to store the logged
information. Python code is used to process the collected logs.

### Building the Library

Logging is _on_ by default
(`LOGGER=yes` in
[`CoEVP/Makefile`](https://github.com/exmatex/CoEVP/blob/adbd900521b4651a7daa9782d695320999f7fb0f/Makefile#L25)). Also,
`lulesh` has been made dependent of the logging library (need link when merged
with master) so that it
it rebuilt if the library is updated. Finally, since we're logging to REDIS, it
must be turned on (`REDIS=yes`) in the build too.

A minimal build (from the `CoEVP` directory) with logging enabled looks like:
```sh
make REDIS=yes LOGGER=yes
```

The logging library also works with MPI:
```sh
make REDIS=yes LOGGER=yes COEVP_MPI=yes
```

### Instrumenting CoEVP for Logging

#### Types of log events

#### Formats of currently collected log events (from CoEVP)

### Running CoEVP with Logging

### Analyzing the Collected Logs

### Bugs and TODOs

* Startup and shutdown of REDIS is broken when running with MPI.
* Does it make sense to connect on multiple ports?
* How much of a hassle is it going to be to set up a Python environment?
* What is the correct balance between key length and number of values per key?
* Will caching be necessary when the runs get big?
* Timers need to be mapped to a keyword so many can run at once.
* Do we need to ensure thread safety?

