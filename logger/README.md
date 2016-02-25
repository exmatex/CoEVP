Simple Logging Library for CoEVP
---

This library, when linked in with CoEVP (and eventually Tabasco), will provide
limited facilities for logging performance and execution data during a run of
the code. The logging library currently uses REDIS to store the logged
information. Python code is used to process the collected logs. No attempt
has been made to build a general purpose logging library--this code is pretty
much CoEVP-specific.

You may ask why use the REDIS database instead of something simpler, like
files? The answer is mostly that I think we should eat our own dog food by
utilizing the same technology that we're using (and promoting) to build the
entire adaptive sampling application. Also, REDIS makes it easier to collect
all the logs in a central "location". Of course, that also means a central
point of interface to the database which could potentially cause performance
problems. However, in the dog food vein, using REDIS for this functionality may
also help us identify and address problems that might eventually crop up in the
full CoEVP/Tabasco application.

Finally, I don't expect anyone else on the project to really use this. I suspect
that I, as author of this stuff, will become the defacto logging and analysis
person--especially since most of the real analysis work is in the Python
scripts anyway.

### Building the Library

Logging is integrated into the CoEVP build system just like Flann, SILO,
etc. It is linkled into the `libcm` library so that it can be used by `lulesh`
and potentially any other code built with `libcm` (e.g. Tabasco).

Logging depends on REDIS, so if you turn on logging (`LOGGER=yes`) you must also turn on
REDIS (`REDIS=yes`). If not, the build will still succeed but logging functions will act as
`NO OP`s.  Note that logging uses its own REDIS server--a different instance of REDIS than CoEVP.

A minimal build (from the `CoEVP` directory) with logging enabled looks like:
```sh
make REDIS=yes LOGGER=yes
```

The logging library also works with MPI:
```sh
make REDIS=yes LOGGER=yes COEVP_MPI=yes
```

### Instrumenting CoEVP and LULESH for Logging

The library logs three types of events:

 * Info: Just logs an arbitrary text string. This might be useful for logging
   the parameters of the current run, its machine configuration, time and date,
   etc.
 * Count: A key/value pair indicating a count for a particular operation.
 * Timer: A key/value pair indicating the elapsed time of a particular
   operation.

The type of log entry to generate is selected using enums (`LOG_INFO`,
`LOG_TIMER`, `LOG_COUNT`). Keys for counts and timers (info events are not
key/value pairs) look like:
```
EVENT,node,id,timestep,func
```
where:
 * `EVENT`: either `COUNT` or `TIMER`
 * `node`: physical node that code was running on when logged (e.g. cn112)
 * `id`: thread/process/rank of code that was logging the event (e.g. MPI rank)
 * `timestep`: the active timestep (as defined by a main loop instrument) when
 an event was logged
 * `func`: the function (string key) that was being counted or timed
 (e.g. `fs_eval`)

Values for count and timer keys differ slightly:

 * `TIMER`: elspsed time in seconds, number of `funcs` timed
 * `COUNT`: number of `funcs` counted

#### Enabling logging

Once in each porcess (logging is not thread safe at this point) logging must be
enabled. Presumind that you compiled logging (with `LOGGER=yes) you'll need to
first include the header files. Obviously, header files must be includedin any
file that uses logging.
```c++
#if defined(LOGGER)      // CoEVP Makefile enforces assert LOGGER=REDIS=yes
#include "LoggerDB.h"    // Includes Logger base class too
#include "Locator.h"
#endif
```

Then, at **one and only one** location in a process (`main` would be locical) you
need to initialize the logger:
```c++
#if defined(LOGGER)
// Initialize logging to REDIS database
Locator::initialize();                   // *** 1 ***
if (logging) {                           // *** 2 ***
#if defined(COEVP_MPI)                   // *** 3 ***
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  char my_node[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(my_node, &name_len);
  LoggerDB  *logger_db = new LoggerDB(logdb, std::string(my_node), my_rank);
#else
  LoggerDB  *logger_db = new LoggerDB(logdb);
#endif
  Locator::provide(logger_db);           // *** 4 ***
}
Logger  &logger = Locator::getLogger();  // *** 5 ***
#endif
```

Some notes on the above. Refer to the numbers in the comments:

 1. The `initialize` call is required to provide a _dummy_ logger. If an attempt
 is made to use a logger that hasn't been enabled (by calling `new LoggerDB`)
 the dummy logger will be substitued and logging will esssentially become a `NO
 OP`.
 2. Even though logging may have been compiled in, the user (at least in CoEVP)
 must direct its usage by specifying a command line flag (here `logging`). As
 shown in this code snippet, the `Locator` is initialized even if we're not
 directed to log events. That's to provide a _dummy_ logger so that logging code
 continues to work (if only as a `NO OP`).
 3. Obviously this is MPI specific code to assign an appropriate `node` and
 `id` (rank in the case of MPI). If we use logging with some other distributed
 computing solution (e.g. Charm++) we'll have to add more specific code for
 that.
 4. This line will provide a _real_ logger, essentially stepping on the default
 _dummy_ logger.
 5. Here we grab a reference to the enabled (or potentially _dummy_) loogger
 for immediate use.

#### Timing the main loop

Since one almost always wants to time the entire code, and look at it on a
per-timestep basis, it is recommended that the user always time the main
loop. In addition to timing the main loop, the user should be sure to always
increment the logger's internal timestep counter so that it can be accurately
tied back to logged events. The main loop is the best place to do this.
Here's an example of you you might instrument the `lulesh` main loop:

```c++
void Lulesh::go(...) {

#if defined(LOGGER)   // did I mention that I hate this define stuff?
  Logger  &logger = Locator::getLogger();
#endif

  /* timestep to solution */
  while(domain.time() < domain.stoptime() ) {

#if defined(LOGGER)
    logger.logStartTimer("outer");
#endif

    //  guts of the inner loop

#if defined(LOGGER)
    logger.logIncrTimer("outer");
    logger.incrTimeStep();
#endif

  }  /* while */
}  /* go  /*
```

#### Logging other events

Now that you've initialized and enabled logging (e.g. at the top of `main) and
instrumented the main loop for timing and, most importantly, incrementing the
timestep, you can use the logger freely anywhere else in the process. Here's an
example of how you might log timers and counts in a critial part of the CoEVP
code:

```c++
#if defined(LOGGER)      // CoEVP Makefile enforces assert LOGGER=REDIS=yes
#include "LoggerDB.h"    // Includes Logger base class too
#include "Locator.h"
#endif

void
AdaptiveSampler::sample(...) {

//  blah blah blah

#if defined(LOGGER)
  Logger  &logger = Locator::getLogger();  // *** 1 ***
  logger.logStartTimer("interpolate");     // *** 2 ***
#endif
  bool interpolationSuccess = m_interp->interpolate(...);
#if defined(LOGGER)
  logger.logIncrTimer("interpolate");      // *** 3 ***
#endif

  if (interpolationSuccess == false) {     // *** interp check FAIL ***
#if defined(LOGGER)
    logger.logIncrCount("interp_fail");    // *** 4 ***
#endif

#if defined(LOGGER)
    logger.logStartTimer("fs_eval");
#endif
     fineScaleModel.evaluate(point, value);
#if defined(LOGGER)
    logger.logIncrTimer("fs_eval");
#endif

#if defined(LOGGER)
    logger.logStartTimer("fs_insert");
#endif
    m_interp->insert(...);
#if defined(LOGGER)
    logger.logIncrTimer("fs_insert");
#endif
  }
  else {                                   // *** interp check SUCCEED ***
#if defined(LOGGER)
    logger.logIncrCount("interp_succeed");
#endif

#if defined(LOGGER)
    logger.logStartTimer("second_interp");
#endif
    interpolationSuccess = m_interp->interpolate(...);
#if defined(LOGGER)
    logger.logIncrTimer("second_interp");
#endif
}                                          // *** end interp check if ***
```

Some notes on the above code:

 1. After including the required header files, get a reference to the enable
 logger. You only need to do this once in a file and the logger is global to
 the process.
 2. Just killing two birds with one stone within the ugly `define`. The logger
 is designed to **NOT** require defines. If logging were to be compiled in by
 default (and hence never not compiled in) the logging functions will always
 execute successfully. Either logging was enabled by the user on the commmand
 line, and an appropriate logger was constructed and enabled, or the dummy
 logger's `NO OP` methods will be called. End of rant.
 3. The `logIncrTimer` call actually does two things: it stops a running timer
 (producing and elapsed time since the `logStartTimer` call) and adds it to the
 cumulative execution time for this `node`, `id`, `timestep`, and `func`, and
 increments an internal count of the number of times this timer has been called
 on this `func`. The rationale is that many functions are called so often, and
 execute so quickly, that it only makes sense to time them in aggregate for a
 given timestep. The internal count is provided so that average time can be
 calculated (and can actually function as an event counter as well). For
 example, a _key_ for this `func` might look like `TIMER,cn117,3,73,interpolate`
 where `cn117` is the node it's executing on, `3` is the MPI rank, `73` is the
 timestep, and `interpolate` is the func. It's corresponding _value_ might be
 `0.035102,672` where 0.035102 seconds were consumed in 672 executions of this
 `func` (at the end of the timestep).
 4. The `logIncrCount` function simply increments the count of `interp_fail`
 `func`s. It's _key_ might look like `COUNT,cn115,1,274,interp_fail` abd its
 _value_ `950`.
 
#### Terminating the logger

Finally, at the end of the main loop (and presumably the code), the user needs
to manually delete the logger to flush remaining data and shut down the
database. The logger's destructor should be called when it falls out for scope,
but for some reason it isn't, and hence this step is ncessary.
```c++
#if defined(LOGGER)
if (logging) {
  delete(&logger);
}
#endif
```

### Running CoEVP with Logging

To run CoEVP with loging the user must:

 1. Compile logging via the `LOGGER=yes` flag
 2. Start an instance of REDIS for logging
 3. Enable logging via the command line and pass in the RESIS database location

On Darwin, the user should ssh  to a node that is part of a SLURM allocation
and start a REDIS server:
```sh
$ ssh cn111
$[cn111] module load redis
$[cn111] redis-server&
```

Note that while `lulesh` will automatically start its required instance(s) of
REDIS (for adaptive sampling) the logger does not yet do that.

To run `lulesh` with logging, adaptive sampling, and using MPI:
```sh
mpirun --map-by node -np 4 ./lulesh -l cn118:6379 -s
```

### Analyzing the Collected Logs

When the job terminates REDIS will leave a dump.rdb file on the node wher it
was started (e.g. `cn111` above). This file can be copied anywhere that REDIS
is available (e.g. your personal MacBook) and processed there. I intend to use
Python to process and analyze logs. The fors step is to start REDIS and convert
the dump.rdb file to a CSV (comma separated value) file that Python can
use. Please visit the `scripts` subdirectory for mor information on the Python
log processing pipeline.


