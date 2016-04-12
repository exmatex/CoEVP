import sys
import os
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt

def genTimings(filePtr):
    scale = 0
    totalTime = 0.0
    retList = [[0, 0.0]]
    for line in filePtr:
        if line.startswith('Timer Output Frequency is'):
            #We are just starting out and need to get the frequency (1)
            strArr = line.split()
            scale = int(strArr[4])
        elif line.startswith('Changing Timer Output Frequency to'):
            #We are changing the frequency
            strArr = line.split()
            scale = int(strArr[5])
        elif not line.startswith('-1'):
            #This is a data line
            strArr = line.split()
            if len(strArr) == 4:
                #This is not the final line of a successful run
                index = int(strArr[0])
                if index == 0:
                    #This is a scale change line
                    index = int(strArr[2].replace(":",""))
                    totalTime = float(strArr[3])
                    retList.append([index, totalTime])
                else:
                    #This is a normal data line
                    index = int(strArr[2].replace(":",""))
                    totalTime += float(strArr[3])
                    retList.append([index, totalTime])
    retArr = np.array(retList)
    return retArr

def main():
    if len(sys.argv) != 2:
        print("Error: Pass in a file")
        return
    inFile = open(sys.argv[1], 'r')
    fileArr = genTimings(inFile)
    return

if __name__ == "__main__":
    main()

