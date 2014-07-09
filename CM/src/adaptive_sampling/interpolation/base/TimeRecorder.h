//
// File:        TimeRecorder.h
// Package:     kriging algorithm
// 
// 
// 
// Description: Class implementing kriging model including derivative 
// Description: information.

#if !defined(included_krigalg_TimeRecorder)
#define included_krigalg_TimeRecorder

#include <ctime>

namespace krigalg {
 
    class TimeRecorder {
    public:
      /*!
       * Default constructor. Records the time upon creation.
       */
      TimeRecorder();

      /*! 
       * Default destructor.
       */
      ~TimeRecorder();
      
      /*!
       * Update time.
       */
      void update();

      /*!
       * Compute difference in seconds between this and other time.
       *
       * @param timeRecorder Other recorded time.
       *
       * @return Time difference in seconds.
       */
      int diff(const TimeRecorder & timeRecorder) const;

    private:
      //
      // copy construction/assignment
      //
      // TimeRecorder(const TimeRecorder &);
      // const TimeRecorder & operator=(const TimeRecorder &);
      
    private:
      time_t _time;
      
    };

}

#ifndef DEBUG_NO_INLINE
#include "TimeRecorder.I"
#endif

#endif // included_krigalg_TimeRecorder
