import sys
import os

REDIS_PORT = 6379
NUTCRACKER_PORT=6380

def buildConfigs(nodeList):
    cwd = os.getcwd();
    for node in nodeList:
        #mkdir for node
        nPath = os.path.join(cwd, node)
        if not os.path.exists(nPath):
            os.mkdir(nPath)
        #$(NODE)/nutcracker.yml file
        yPath = os.path.join(nPath, "nutcracker.yml")
        yFile = open(yPath, 'w')
        #write yaml file
        yFile.write(node + "_nut:\n")
        yFile.write("  listen: " + node + ":"+ str(NUTCRACKER_PORT) + "\n")
        yFile.write("  hash: fnv1a_64\n")
        yFile.write("  distribution: ketama\n")
        yFile.write("  auto_eject_hosts: true\n")
        yFile.write("  redis: true\n")
        yFile.write("  server_retry_timeout: 2000\n")
        yFile.write("  server_failure_limit: 1\n")
        yFile.write("  servers:\n")
        for rServer in nodeList:
            rHost = rServer + ":" + str(REDIS_PORT) + ":1"
            yFile.write("   - " + rHost + "\n")
        yFile.write("\n")
        yFile.close()


def genSlurmList(nodePrefix):
    # Get environmental variable
    nodeString = os.environ['SLURM_JOB_NODELIST']
    # For debug/devel
    # Remove leading 'prefix'
    nodeString = nodeString.lstrip(nodePrefix)
    nodeList = []
    if (nodeString.startswith('[') == False):
        # Just the one
        nodeList.append(int(nodeString))
    else:
        # Strip the brackets, we no longer need them
        nodeString = nodeString.strip('[]')
        # Tokenize the string
        tokenList = nodeString.rsplit(',')
        # Iterate over tokens
        for token in tokenList:
            # Two cases: Single or range
            if '-' in token:
                inclusiveRange = token.rsplit('-')
                start = int(inclusiveRange[0])
                end = int(inclusiveRange[1]) + 1
                for node in xrange(start, end):
                    nodeList.append(node)
            else:
                nodeList.append(int(token))
    hostList = []
    for node in nodeList:
        host = nodePrefix + str(node)
        hostList.append(host)
    return hostList


#This is for cray environments
def genCrayList(nodePrefix):
    nFilePath = os.environ['PBS_NODEFILE']
    nFile = open(nFilePath, 'r')
    nodeList = []
    for line in nFile:
        nodeList.append(int(line))
    hostList = []
    for node in nodeList:
        host = nodePrefix + str(node)
        hostList.append(host)
    return hostList


def selectServerNodes(nodeList, numServers):
    numNodes = len(nodeList)
    offset = 0
    if numServers >= numNodes:
        offset = 1
    else:
        offset = numNodes / numServers
    serverList = []
    if offset == 1:
        for i in xrange(0, min(numNodes, numServers)):
            serverList.append(nodeList[i])
    else:
        for i in xrange(0, numNodes, offset):
            serverList.append(nodeList[i])
    return serverList


def main():
    if(len(sys.argv) != 4):
        print("python buildNutcracker.py $(ENVIRONMENT) $(PREFIX) $(NUMBER_OF_SERVERS)")
        exit(1)
    #Get nodelist
    environ = sys.argv[1]
    prefix = sys.argv[2]
    numberOfServers = int(sys.argv[3])
    fullList = []
    #fullList = genSlurmList("cn")
    #fullList = genCrayList("nid000"i)
    if environ == "SLURM":
        fullList = genSlurmList(prefix)
    elif environ == "cray":
        fullList = genCrayList(prefix)
    else:
        print("We only support SLURM and cray currently")
        exit(1)
    #print(fullList)
    #Select which nodes are servers
    serverList = selectServerNodes(fullList, numberOfServers)
    #print(serverList)
    #Make configs
    buildConfigs(serverList)


if __name__ == "__main__":
    main()
