from __future__ import print_function
import pandas as pd
import numpy as np


def main():
    df = pd.read_csv('dump.csv')
    print(df.info())
    print("\n\n")
    print(df.dtypes)
    print("\n\n")
    print(df.describe())
    print("\n\n")
    timers =   df[(df.Mode == 'TIMER')]
    counters = df[(df.Mode == 'COUNT')]
    print(timers.head())
    print("\n\n")
    print(counters.head())

    counter_ops = counters.groupby(['Timestep', 'Operation'])
    print(counter_ops.aggregate(np.sum))
    print("\n\n")

    timer_ops = timers.groupby(['Timestep', 'Rank', 'Operation'])
    print(timer_ops.aggregate(np.sum))
    print("\n\n")

    
if __name__ == "__main__":
    main()
