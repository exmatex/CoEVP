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
    print("\n\n")

    outer = df[(df.Operation == 'outer')]
    print(outer.head())
    print("\n\n")

    print(outer.groupby('Rank').groups.keys())
    print("\n\n")

    print(outer.groupby('Rank')['Value'].quantile(.99))
    print("\n\n")

    
if __name__ == "__main__":
    main()
