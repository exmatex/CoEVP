from __future__ import print_function
import pandas as pd


def main():
    df = pd.read_csv('dump.csv')
    print(df)
    

if __name__ == "__main__":
    main()
