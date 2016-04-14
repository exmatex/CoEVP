import sys
import os
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt

lineStyles = ['--', ':', '-', '-.']

markerStyles = {}
markerStyles[1] = 'bo'
markerStyles[16] = 'gD'
markerStyles[128] = 'rv'
markerStyles[256] = 'c^'
markerStyles[512] = 'bx'
markerStyles[1024] = 'mh'
markerStyles[2048] = 'kh'
markerStyles[416] = 'm*'
markerStyles[208] = 'm+'

def genAbsTimings(filePtr):
    scale = 0
    totalTime = 0.0
    retList = [[0, 0.0]]
    for line in filePtr:
        strArr = line.split()
        if strArr[0][0].isdigit():
            #This is a data line
            if len(strArr) == 5:
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
        elif line.startswith("Total Cycles:"):
            index = int(strArr[2])
            totalTime = float(strArr[4])
            retList.append([index, totalTime])
    retArr = np.array(retList)
    return retArr

def genStepTimings(filePtr):
    scale = 0
    totalTime = 0.0
    retList = [[0, 0.0]]
    for line in filePtr:
        strArr = line.split()
        if strArr[0][0].isdigit():
            #This is a data line
            if len(strArr) == 5:
                startDex = int(strArr[0])
                endDex = int(strArr[2].replace(":",""))
                time = float(strArr[3]) / int(endDex - startDex)
                retList.append([endDex, time])
        elif line.startswith("Total Cycles:"):
            endDex = int(strArr[2])
            startDex = retList[len(retList)-1][0]
            time = float(strArr[4]) / int(endDex - startDex)
            retList.append([endDex, time])
    retArr = np.array(retList)
    return retArr

def plotAbsTime(prefixList, absList, maxSteps):
    for (name, data) in absList:
        prefCnt = 0
        for (prefix, tag) in prefixList:
            if name.startswith(prefix) and name[len(prefix)].isdigit():
                arrLen = 0
                if maxSteps != 0:
                    for i in xrange(data.shape[0]):
                        if data[i,0] <= maxSteps:
                            arrLen = arrLen + 1
                else:
                    arrLen = data.shape[0]
                strArr = name.split(prefix)
                nodeCount = int(strArr[1].split('.timer')[0])
                lineLabel = tag + ": " + str(nodeCount) + " nodes"
                lineStyle = lineStyles[prefCnt] + markerStyles[int(nodeCount)]
                plt.plot(data[:arrLen, 0], data[:arrLen, 1], lineStyle, label=lineLabel)
            prefCnt = prefCnt + 1
    lgd = plt.legend(loc='lower center', prop={'size':6}, ncol=len(prefixList),
            bbox_to_anchor=(float(1.0 / len(prefixList)),-0.3))
    plt.xlabel('Iterations')
    plt.ylabel('Execution Time (s)')
    plt.savefig("absplot.pdf", format='pdf',  bbox_extra_artists=(lgd,),
            bbox_inches='tight')
    plt.figure()
    return

def plotStepTime(prefixList, stepList, maxSteps):
    for (name, data) in stepList:
        prefCnt = 0
        for (prefix, tag) in prefixList:
            if name.startswith(prefix) and name[len(prefix)].isdigit():
                arrLen = 0
                if maxSteps != 0:
                    for i in xrange(data.shape[0]):
                        if data[i,0] <= maxSteps:
                            arrLen = arrLen + 1
                else:
                    arrLen = data.shape[0]
                strArr = name.split(prefix)
                nodeCount = int(strArr[1].split('.timer')[0])
                lineLabel = tag + ": " + str(nodeCount) + " nodes"
                lineStyle = lineStyles[prefCnt] + markerStyles[int(nodeCount)]
                plt.plot(data[:arrLen, 0], data[:arrLen, 1], lineStyle, label=lineLabel)
            prefCnt = prefCnt + 1
    lgd = plt.legend(loc='lower center', prop={'size':6}, ncol=len(prefixList),
            bbox_to_anchor=(float(1.0 / len(prefixList)),-0.3))
    plt.yscale('log')
    plt.xlabel('Iterations')
    plt.ylabel('Execution Time Per Step (s)')
    plt.savefig("stepplot.pdf", format='pdf',  bbox_extra_artists=(lgd,),
            bbox_inches='tight')
    plt.figure()
    return

def main():
    if len(sys.argv) != 2:
        print("Error: Pass in a directory")
        return
    inPath = sys.argv[1]
    absList = []
    stepList = []
    for root, dirs, files in os.walk(inPath):
        for fName in files:
            if fName.endswith('.timer'):
                fPath =  os.path.join(root, fName)
                inFile = open(fPath, 'r')
                absArr = genAbsTimings(inFile)
                inFile.close()
                if absArr.shape[0] > 1:
                    absList.append((fName, absArr))
                inFile = open(fPath, 'r')
                stepArr = genStepTimings(inFile)
                inFile.close()
                if stepArr.shape[0] > 1:
                    stepList.append((fName, stepArr))
    #print [i[0] for i in absList]
    prefixList = []
    #prefixList.append(("_strong1_charmpp_adapt_128x416:416_0.1000:0_mtree_local:hashmap_",
    #    "Strong1"))
    prefixList.append(("_taylDBSweep_charmpp_adapt_300x256:256_1.0000:0_mtree_local:hashmap_",
        "Adaptive Hashmap"))
    prefixList.append(("_taylDBSweep_charmpp_brute_300x256:256_1.0000:0_mtree_local:hashmap_",
            "Brute Force"))
    plotAbsTime(prefixList, absList, 200)
    plotStepTime(prefixList, stepList, 200)
    return

if __name__ == "__main__":
    main()

