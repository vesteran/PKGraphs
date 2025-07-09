import os

#Gets all Graphs from the paths that are listed in a file. and returns them as a tuple
##The first part of the tuple is the Segment Graphs for the PKs and the second is the PKGraphs
##In these lists the indices will match up if they are pointing to the same RNA structure
def GetGraphsFromFile(file):
    try:
        GraphList = open(file, "r")
    except:
        print("invalid file path")
        return
    curline = GraphList.readline()
    PKSegmentsGraphs = []
    PKGraphs = []
    FileNames = []
    while curline :
        curline = curline.strip()
        (curSegmentGraph, curPKGraph)  = makeGraphsFromFile(curline)
        if(curPKGraph):
            PKGraphs.append(curPKGraph)
            curSegmentGraph = trimSegments(curSegmentGraph)
            PKSegmentsGraphs.append(curSegmentGraph)
            FileNames.append(curline)
        curline = GraphList.readline()
    
    return (PKSegmentsGraphs, PKGraphs, FileNames)

def makeGraphsFromFile(filePath):
    pkGraph = []
    pkSegmentGraph = []

    try:
        file = open(filePath, "r")
    except:
        print("invalid file path for Graph" + filePath)
        return None
    
    lines = file.readlines()
    i = findLineOfString("#PK-graph", lines)
    curline = lines[i]
    while not curline.startswith("#"):
        i += 1
        curedge = curline.split()
        pkGraph.append((int(curedge[0]),int(curedge[1])))
        curline = lines[i]

    i = findLineOfString("#All-segments",lines)
    curline = lines[i]
    while not curline.startswith("#"):
        i += 1
        pk = curline.split()
        bp = (pk[0].split(".."))
        bp1 = bp[0]
        bp = (pk[1].split(".."))
        bp2 = bp[1]
        pkSegmentGraph.append((int(bp1),int(bp2)))
        curline = lines[i]
    return (pkSegmentGraph, pkGraph)

####no longer used###
# def makePKGraphFromFile(filePath):
#     pkGraph = []
#     try:
#         file = open(filePath, "r")
#     except:
#         print("invalid file path for PK Graph")
#         return None
    
#     lines = file.readlines()
#     i = findLineOfString("#PK-graph", lines)
#     curline = lines[i]
#     while not curline.startswith("#"):
#         i += 1
#         curedge = curline.split()
#         pkGraph.append((int(curedge[0]),int(curedge[1])))
#         curline = lines[i]
#     return pkGraph

# def makePKSegmentGraphFromFile(filePath):
#     pkSegmentGraph = []
#     try:
#         file = open(filePath, "r")
#     except:
#         print("invalid file path for PK Segment Graph")
#         return None
#     lines = file.readlines()
#     HasPK = -5
#     for line in lines:
#         if ("PK" in line):
#             HasPK = 1
    
#     if HasPK == -5:
#         return None

#     i = findLineOfString("PK",lines) - 1
#     curline = lines[i]
#     while not curline.startswith("PK1."):
#         i += 1
#         pk = curline.split()
#         bp = (pk[2].split(".."))
#         bp1 = bp[0]
#         bp = (pk[3].split(".."))
#         bp2 = bp[1]
#         pkSegmentGraph.append((int(bp1),int(bp2)))
#         curline = lines[i]

#     segmentGraph = []
#     i = findLineOfString("segment",lines) - 1
#     curline = lines[i]
#     while curline.startswith("segment"):
#         pk = curline.split()
#         bp = (pk[2].split(".."))
#         bp1 = bp[0]
#         bp = (pk[4].split(".."))
#         bp2 = bp[1]
#         segmentGraph.append((int(bp1),int(bp2)))
#         i = i + 1
#         if(i < len(lines)):
#             curline = lines[i]
#         else:
#             curline = ""

#     pkSegmentGraph = trimSegments(pkSegmentGraph, segmentGraph)
#     return pkSegmentGraph

#a segment is trimable if it is not crossed by a pseudoknot
#Therefore if a segment crosses a pseudoknot it is added to the graph to calculate genus. 
#You can check if a PK crosses a segment by seeing if either one but not both of the endpoints of a pk are inside the segment
def trimSegments(SegGraph):
    pkSegGraph = []
    for seg in SegGraph:
        SegisPk = False
        for checkseg in SegGraph:
            if (seg[0] < checkseg[0] < seg[1] < checkseg[1]) or (checkseg[0] < seg[0] < checkseg[1] < seg[1]):
                SegisPk = True
        if SegisPk == True:
            pkSegGraph.append(seg)
    return pkSegGraph

def findLineOfString(string, fileLines):
    n = 0
    for line in fileLines:
        n += 1
        if line.startswith(string):
            return n

    return None

# make tsv output for parsed pks w/ genus, type and name, size of graph
