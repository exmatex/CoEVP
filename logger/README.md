Using the Logging Library
---

This library, when linked in with CoEVP (and eventually Tabasco), will provide
limited facilities for logging performance and execution data during a run of
the code. The logging library currently uses REDIS to store the logged
information. Python code is used to process the collected logs.

### Building the Library

Logging is _on_ by default
(`LOGGER=yes` in [`CoEVP/Makefile`](https://github.com/exmatex/CoEVP/blob/adbd900521b4651a7daa9782d695320999f7fb0f/Makefile#L25). 

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
