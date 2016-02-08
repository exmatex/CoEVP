import sys
import os
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt


def main():
    argList = [ ("sf", "HashMap + FLANN"),
                ("sgf", "Global HashMap + FLANN"),
                ("s", "Hashmap + MTree"),
                ("sgrf", "Global Redis + FLANN"),
                ("sg", "Global HashMap + MTree"),
                ("sgr", "Global Redis + MTree"),
                ("brute", "No Adaptive Sampling")]
    #argList = ["sf", "sgf", "s", "sgrf", "sg", "sgr", "brute"]
    #argList = ["sf", "sgf", "s", "sgrf", "brute"]
    cwd = os.getcwd()
    f, axarr = plt.subplots(2, sharex=True)
    for (arg, wordies) in argList:
        subdir = os.path.join(cwd, arg)
        fPath = os.path.join(subdir, "runTimes.dat")
        data = np.loadtxt(fPath)
        xVals = data[:, 0]
        timeVals = data[:, 3]
        callVals = data[:, 5]
        axarr[0].plot(xVals, timeVals, label=wordies)
        axarr[1].plot(xVals, callVals, label= wordies)
    plt.xlabel("Iteration")
    axarr[0].set_ylabel("Time (s)")
    #axTime.title("Time per Iteration for Different Arguments")
    #axTime.legend()
    #axCalls.xlabel("Iteration")
    axarr[1].set_ylabel("Fine Scale Calls")
    #axCalls.title("Fine Scale Calls  per Iteration for Different Arguments")
    #axCalls.legend()
    lgd = plt.legend(loc='lower center', bbox_to_anchor=(0.5, -0.83), ncol=2)
    #axarr[0].legend(loc='upper right')
    plt.title("Time and Calls per Iteration for Different Arguments")
    outFile = os.path.join(cwd, "plot.pdf")
    plt.savefig(outFile, bbox_extra_artists=(lgd,), bbox_inches='tight')
    return


if __name__ == "__main__":
    main()
