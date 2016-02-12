from __future__ import print_function
import pandas as pd


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
    print(timers)
    print("\n\n")
    print(counters)
    print("\n\n")

    
if __name__ == "__main__":
    main()
