import sys
import os
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt


def main():
    argList = ["sf", "sgf", "s", "sg", "sgrf", "sgr"]
    cwd = os.getcwd()
    f, axarr = plt.subplots(2, sharex=True)
    for arg in argList:
        subdir = os.path.join(cwd, arg)
        fPath = os.path.join(subdir, "runTimes.dat")
        data = np.loadtxt(fPath)
        xVals = data[:, 0]
        timeVals = data[:, 3]
        callVals = data[:, 4]
        axarr[0].plot(xVals, timeVals, label="-" + arg)
        axarr[1].plot(xVals, callVals, label="-" + arg)
    plt.xlabel("Iteration")
    axarr[0].set_ylabel("Time (s)")
    #axTime.title("Time per Iteration for Different Arguments")
    #axTime.legend()
    #axCalls.xlabel("Iteration")
    axarr[1].set_ylabel("Fine Scale Calls")
    #axCalls.title("Fine Scale Calls  per Iteration for Different Arguments")
    #axCalls.legend()
    axarr[0].legend(loc='upper right')
    plt.title("Time and Calls per Iteration for Different Arguments")
    outFile = os.path.join(cwd, "plot.pdf")
    plt.savefig(outFile)
    return


if __name__ == "__main__":
    main()
