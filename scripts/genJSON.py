import sys
import os

#Vars we won't be changing
visInterval = 20
remoteDB = 0
resultPath = "./results"

def numFineScaleCalls(edge, height):
    core = edge / 4
    wing = edge - core
    nElem = (core * edge * height) + (core * height * wing)
    return nElem

def genName(eElems, hElems, sTime, domains, sample, remoteDB, db, nn, dw, nSteps, nNodes, rt):
    nameStr = ""
    if(rt == 0):
        nameStr += "charmpp_"
    elif(rt == 1):
        nameStr += "mpi_"
    else:
        nameStr += "error_"
    if (sample == 1):
        nameStr += "adapt"
    else:
        nameStr += "brute"
    nameStr += "_"
    nameStr += str(eElems) + "x" + str(hElems) + ":"
    nameStr += str(domains) + "_"
    nameStr += '{0:.4f}'.format(sTime) + ":"
    nameStr += str(nSteps) + "_"
    if (nn == 0):
        nameStr += "flann"
    else:
        nameStr += "mtree"
    nameStr += "_"
    if (remoteDB == 0):
        nameStr += "local"
    else:
        nameStr += "remote"
    nameStr += ":"
    if (db ==  0):
        nameStr += "redis"
    elif (db == 1):
        nameStr += "hashmap"
    elif (db == 2):
        nameStr += "posix"
    elif (db == 3):
        nameStr += "hio"
    elif (db == 4):
        nameStr += "twemproxy"
    else:
        nameStr += "DB_ERROR"
    if (dw == 1):
        nameStr += "_DW"
    nameStr += "_" + str(nNodes)
    return nameStr

def genJSON(tag, eElems, hElems, simTime, domains, sample, remoteDB, db, nn, dw, nSteps):
    fileName = tag + '.json'
    numFS = numFineScaleCalls(eElems, hElems)
    #Boilerplate
    jsonFile = open(fileName, 'w')
    jsonFile.write('{"parameter": {\n')
    jsonFile.write('"header": "TabaSCo Test Input File",\n')
    #Coarse Scale Model
    jsonFile.write('"CoarseScaleModel": [\n')
    jsonFile.write('{"id": "type", "value": "0"},\n')
    jsonFile.write('{"id": "count", "value": "' + str(domains) + '"},\n')
    jsonFile.write('{"id": "use adaptive sampling", "value": "' + str(sample) + '"},\n')
    jsonFile.write('{"id": "stop time", "value": "' + str(simTime) + '"},\n')
    jsonFile.write('{"id": "max steps", "value": "' + str(nSteps) + '"},\n')
    jsonFile.write('{"id": "visit data interval",   "value": "' + str(visInterval) + '"},\n')
    jsonFile.write('{"id": "file parts","value": "4"},\n')
    jsonFile.write('{"id": "edge elems","value": "' + str(eElems) + '"},\n')
    jsonFile.write('{"id": "height elems", "value": "' + str(hElems) + '"},\n')
    jsonFile.write('{"id": "timer sampling rate", "value": "1"}\n')
    jsonFile.write('],\n')
    #Fine Scale Model
    jsonFile.write('"FineScaleModel": [\n')
    jsonFile.write('{"id": "type", "value": "0"}\n')
    jsonFile.write('],\n')
    #Nearest Neighbor Search
    jsonFile.write('"NearestNeighborSearch": [\n')
    jsonFile.write('{"id": "type", "value": "' + str(nn) + '"},\n')
    jsonFile.write('{"id": "count", "value": "' + str(numFS) + '"},\n')
    jsonFile.write('{"id": "point dimension", "value": "6"},\n')
    jsonFile.write('{"id": "number of trees", "value": "1"}\n')
    jsonFile.write('],\n')
    #Interpolate (?)
    jsonFile.write('"Interpolate": [\n')
    jsonFile.write('{"id": "type", "value": "0"},\n')
    jsonFile.write('{"id": "count", "value": "10"}\n')
    jsonFile.write('],\n')
    #Evaluate
    jsonFile.write('"Evaluate": [\n')
    jsonFile.write('{"id": "type", "value": "0"},\n')
    jsonFile.write('{"id": "count", "value": "' + str(numFS) + '"}\n')
    jsonFile.write('],\n')
    #DBInterface
    jsonFile.write('"DBInterface": [\n')
    jsonFile.write('{"id": "type", "value": "' + str(db) + '"},\n')
    jsonFile.write('{"id": "count", "value": "10"},\n')
    jsonFile.write('{"id": "remote", "value": "' + str(remoteDB) + '"}\n')
    jsonFile.write(']\n')
    #Boilerplate
    jsonFile.write('}\n}\n')
    jsonFile.close()
    return

def main():
    dbType = 2 #0 = redis, 1 = hashmap, 2=posix, 3=hio, 4=twemproxy 
    nnType = 1 #mtree? VERIFY
    sampling = 0 #Default to brute force
    domains = 1 #One domain for now
    simTime = 1e-01 #This will be fun to set from the command line...
    edgeElems = 32
    heightElems = 52
    nNodes = 64
    rt = 0 #0 = charm, 1 = mpi. 2 = circle or liblouis, whichever is first
    dw = 0 #No datawarp by default
    nSteps = 100 #Don't limit by steps, by default
    #TODO: Add in command line args... or just do a predefined sweep
    if not os.path.exists(resultPath):
        os.mkdir(resultPath)
    #TODO: Put the for loops here, probably
    tag = genName(edgeElems, heightElems, simTime, domains, sampling, remoteDB, dbType,
            nnType, dw, nSteps, nNodes, rt)
    tagPath = os.path.join(resultPath, tag)
    if not os.path.exists(tagPath):
        os.mkdir(tagPath)
    os.chdir(tagPath)
    if rt == 0:
        genJSON(tag, edgeElems, heightElems, simTime, domains, sampling, remoteDB,
                dbType, nnType, dw, nSteps)
    return



if __name__ == "__main__":
    main()
