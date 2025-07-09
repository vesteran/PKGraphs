import sys
import os
from datetime import datetime
import ReadGraphs as reader
import numpy as np
import csv

def main():
    
    # Attempt to grab graphs from file provided and exits program if it fails. 
    if  len(sys.argv) <= 1:
        print("No file inputed")
        return
    try:
        (SegGraph, PKGraph, RNANames) = reader.GetGraphsFromFile(sys.argv[1])
    except:
        print("Unable to retrieve file paths from file provided.")
        return

    if(SegGraph == None or PKGraph == None):
        print("One or More paths in provided paths file are invalid")
        return
    
    OverallStructureTSV = []
    IndividualPsuedoKnotTSV = []

    for i in range(len(SegGraph)):
        curRNAName = RNANames[i]
        try: 
            # curGraphsPKs = [(curSegmentGraph,PKGraph[i])]
            totalGenus = calculatePKGenus(SegGraph[i])
            print(identifyPseudoKnots(PKGraph[i],totalGenus))

            OverallStructureTSV.append([curRNAName,totalGenus,len(SegGraph[i]),SegGraph])

            (splitSegGraph, splitPKGraph) = seperatePseudoknots(SegGraph[i],PKGraph[i])
            for i in range(len(splitSegGraph)):

                #Sorts segment graph by first BP so we can corrolate it to the PK after seperating the knots
                SegmentGraphSortedIndices = np.argsort(splitSegGraph[i],0)
                curSegmentGraph = []
                for ind in SegmentGraphSortedIndices:
                    curSegmentGraph.append(splitSegGraph[i][ind[0]])
                splitSegGraph[i] = curSegmentGraph

                (splitSegGraph[i],splitPKGraph[i]) = mergeStems(splitSegGraph[i],splitPKGraph[i])
                print(splitSegGraph[i])
                print(splitPKGraph[i])

                genus = calculatePKGenus(splitSegGraph[i])
                print(identifyPseudoKnots(splitPKGraph[i], genus))
                IndividualPsuedoKnotTSV.append([curRNAName, i, genus, len(splitPKGraph[i]),splitSegGraph[i]])
        except:
            print("An Error occured while classifying the pseudoknot from file " + curRNAName)

    makeTSVs(OverallStructureTSV, IndividualPsuedoKnotTSV)
    return



###############################################################################
# Calculates the Genus of an inputed RNA segment graph. 
# Input: Segment Graph
# Output: Calculated Genus (Int)
# Notes: Genus Formula is G = (P-L)/2 where P is the number of stems and L is the number of internal Loops
###############################################################################

def calculatePKGenus(graph):
    
    maxbp = np.max(graph)
    Visitedbps = []
    internalLoops = -1
    for edge in graph:
        curedge = edge
        curbp = curedge[0]
        startingbp = curbp
        if(curbp not in Visitedbps):
            while (curbp not in Visitedbps): #Travels in a path along the edges and the bottom of the edge graph untill it makes a loop or encounters a previously encountered edge
                Visitedbps.append(curbp)
                if (curbp == curedge[0]):
                    curbp = curedge[1]
                else:
                    curbp = curedge[0]

                if curbp >= maxbp:          #If the bp is the highest of all the bps then we loop it back to the beginning. 
                    curbp = -1

                testEdge = (-5,-5)
                while curbp not in testEdge:    #Iterate along the line untill we get to another edge then we set it as our next edge
                    curbp += 1
                    for edge_ in graph:
                        if curbp in edge_:
                            curedge = edge_
                            testEdge = curedge
                    
                
            if(startingbp == curbp):
                internalLoops += 1


    print("genus = (" + str(len(graph)) + " - " + str(internalLoops) + ") / 2")
        
    if((len(graph) - internalLoops)%2 == 0):
        return (len(graph) - internalLoops)/2
    else:
        print("Genus calculated as non whole number")
        print(graph)
        exit()

###############################################################################
# Outputs either the denoting letter for pks of genus 1 or it outputs the genus and number of crossings to the terminal. 
# Input: PK Graph and Genus
# Output: None (Outputs to terminal)
###############################################################################

def identifyPseudoKnots(graph, genus):
    if(len(graph) == 1):
        return "H"
    elif(len(graph) == 2 and genus == 1):
        return "K"
    elif(len(graph) == 3 and genus == 1):
        return "L"
    elif genus == 1:
        return 'M'
    else:
        return ("Pseudo Knot of Genus " + str(genus) + " consisting of " + str(len(graph)) + " crossings")


###############################################################################
# Seperates the different pseudoknot structures contained in an RNA structure.  
# Input: Segment Graph, PK Graph
# Output: A list of segments graphs for each seperate component and the same for the PK graph. 
# Notes: A structure is splittable into seperate pseudo knots if there is more then one individual component in the pkgraph
#        The outputs lists will be ordered such such that segmentgraph[1] will corrospond to pkgraph[1]
###############################################################################
def seperatePseudoknots(SegGraph, PKGraph):
    SplitPKGraphs = []
    SplitSegGraphs = []
    visitedSegs = []
    segs_to_visit = SegGraph.copy()
    
    while segs_to_visit:
        cur_pk = []
        curseg = segs_to_visit.pop()
        if curseg not in visitedSegs:
            visitedSegs.append(curseg)
            segs_to_visit2 = get_seg_crossings(SegGraph,curseg)
            if segs_to_visit2:
                cur_pk = [curseg]
            

            while segs_to_visit2:
                curseg = segs_to_visit2.pop()
                if curseg not in cur_pk:
                    cur_pk.append(curseg)
                    visitedSegs.append(curseg)
                    segs_to_visit2.extend(get_seg_crossings(SegGraph,curseg))

            if cur_pk:
                SplitSegGraphs.append(cur_pk)

    # Recreate the pk graphs from the split segment graphs 
    for subseggraph in SplitSegGraphs:
        subpkgraph = []
        for i in range(len(subseggraph)):
            for j in range(i + 1, len(subseggraph)):
                seg1 = subseggraph[i]
                seg2 = subseggraph[j]
                if(seg2[0] <= seg1[0] <= seg2[1] <= seg1[1]) and ((j,i) not in subpkgraph):
                    subpkgraph.append((j,i))
                if(seg1[0] <= seg2[0] <= seg1[1] <= seg2[1]) and ((i,j) not in subpkgraph):
                    subpkgraph.append((i,j))
        SplitPKGraphs.append(subpkgraph)
    return (SplitSegGraphs,SplitPKGraphs)

### No longer used
#Takes a graph in the form of a list of edges [(1,2),(2,3)...] and returns a True if it has a cycle and a false otherwise. 
# def checkForCycles(Graph):
#     visitedEdges = []
#     visitedEdges = []
#     visitedEdges.append(Graph[0][0])
#     frontier = getConnectedEdgesFromNode(Graph,Graph[0][0])
#     while frontier:
#         curEdge = frontier[0]
#         frontier.pop(0)
#         visitedEdges.append(curEdge)
#         if curEdge[0] in visitedEdges and curEdge[1] in visitedEdges:
#             return True
#         elif not curEdge[0] in visitedEdges:
#             visitedEdges.append(curEdge[0])
#             for edge in getConnectedEdgesFromNode(Graph,curEdge[0]):
#                 if not edge in visitedEdges:
#                     frontier.append(edge)
#         elif not curEdge[1] in visitedEdges:
#             visitedEdges.append(curEdge[1])
#             for edge in getConnectedEdgesFromNode(Graph,curEdge[1]):
#                 if not edge in visitedEdges:
#                     frontier.append(edge)
#     return False



###############################################################################
# Merges any stems that can be merged while maintaining genus
# Input: Segment Graph, PK Graph
# Output: A trimmed segment graph and pk graph. 
# Notes: 2 segments can be merged while preserving the genus if and only if the cross over the exact same other segments and are nested. 
#        Example: (( [ )) ] -> ([)] but not ([)(]) or [(]( [ )) ]
###############################################################################    
def mergeStems(SegGraph, pkGraph):
    curgraph = SegGraph
    trimmedPKGraph = []
    trimmedSegGraph = []
    SegsToNotAdd = []
    cur_len = len(SegGraph)
    prev_len = cur_len + 1
    
    while(cur_len != prev_len):
        trimmedSegGraph = []
        SegsToNotAdd = []
        for seg in curgraph: 
            for seg2 in curgraph:
                if(seg2 not in SegsToNotAdd and seg not in SegsToNotAdd and segsAreStacked(curgraph,seg2,seg)):
                    SegsToNotAdd.append(seg2)

        for seg in curgraph:
            if seg not in SegsToNotAdd:
                trimmedSegGraph.append(seg)
        prev_len = cur_len
        cur_len = len(trimmedSegGraph)
        curgraph = trimmedSegGraph
        
    #reconstruct pkGraph from new seg graph as that is easier then trying to make sure we only remove the correct nodes from pk as well
    for i in range(len(trimmedSegGraph)):
        for j in range(i + 1, len(trimmedSegGraph)):
            seg1 = trimmedSegGraph[i]
            seg2 = trimmedSegGraph[j]
            if(seg2[0] <= seg1[0] <= seg2[1] <= seg1[1]) and ((j,i) not in trimmedPKGraph):
                trimmedPKGraph.append((j,i))
            if(seg1[0] <= seg2[0] <= seg1[1] <= seg2[1]) and ((i,j) not in trimmedPKGraph):
                trimmedPKGraph.append((i,j))

    return (trimmedSegGraph, trimmedPKGraph)
    

###############################################################################
# Tells us if any two specefic segments are stacked
# Input: segment graph and the two segments to check.
# Output: True or False
# Notes: 2 segments can be merged while preserving the genus (are stacked) if and only if the cross over the exact same other segments and are nested. 
#        Example: (( [ )) ] -> ([)] but not ([)(]) or [(]( [ )) ]
###############################################################################    
def segsAreStacked(graph,seg1,seg2):
    #orders segments such that seg1 alwas starts at a lower bp then seg 2
    if(seg1[0] > seg2[0]):
        temp = seg1
        seg1 = seg2
        seg2 = temp

    #checks to see if they cross each other
    if (seg1[0] < seg2[0] < seg1[1] < seg2[1]):
        return False

    #Makes sure that they are nested. 
    if not(seg1[0] < seg2[0] and seg2[1] < seg1[1]):
        return False

    #checks to see if there are any other segments splitting the two we are checking. 
    for curseg in graph:
        if(curseg != seg1 and curseg != seg2):
            if (seg1[0] <= curseg[0] <= seg2[0]) or (seg1[0] <= curseg[1] <= seg2[0]):
                return False
            elif (seg2[1] <= curseg[0] <= seg1[1]) or (seg2[1] <= curseg[1] <= seg1[1]):
                return False
    return True

###############################################################################
# Tells us all of the segments that are crossing a specified other segment. 
# Input: segment graph and the segments to find crossing for.
# Output: a list of the segments that cross
# Notes: 2 segments cross if one of them starts outside of the other but ends inside or vice versa
###############################################################################   
def get_seg_crossings(SegGraph, seg):
    crossingSegs = []
    for seg2 in SegGraph:
        if seg2 != seg and (seg2[0] <= seg[0] <= seg2[1] <= seg[1] or seg[0] <= seg2[0] <= seg[1] <= seg2[1]):
            crossingSegs.append(seg2)
    return crossingSegs

### No longer Used

#Grabs all nodes that are connected by an edge to a node in a graph of the form [edge1, edge2] where each edge is of the form (node1,node2)
# def getConnectedNodesFromNode(graph, node):
#     connectedNodes = []       
#     for edge in graph:
#         if node == edge[0] :
#             connectedNodes.append(edge[1])
#         elif node == edge[1] :
#             connectedNodes.append(edge[0])
#     return connectedNodes
# 
#Grabs all connected edges to a node in a graph of the form [edge1, edge2] where each edge is of the form (node1,node2)
# def getConnectedEdgesFromNode(graph, node):
#     connectedEdges = []       
#     for edge in graph:
#         if node in edge :
#             connectedEdges.append(edge)
#     return connectedEdges
# 
# def listbps(graph):
#     bps = []
#     for edge in graph:
#         if not edge[0] in bps:
#             bps.append(edge[0])
#         if not edge[1] in bps:
#             bps.append(edge[1])
#     return   bps

###############################################################################
# Creates TSV's with the data we calculated. 
# Input: A list of overall genus data [[RNA_Name, Genus, Number of Segments], ...]
#        and a list of individual genus data for each individual pseudoknot [[RNA_Name, PK ID,Genus, Number of Crossings], ...]
# Output: None(Creates 2 files in Output_Files labelled with the current time)
# Notes: Files created will be tab delimited tsv files 
###############################################################################  

def makeTSVs(Overall, Individual):
    curtime = str(datetime.now()).replace(" ","_").replace(":","-")
    os.makedirs("Output_Files/" + curtime)

    Fields = ["File Name", "Genus", "Number of Stems","Seg Graph"]
    FileName = "Output_Files/"+ curtime + "/Overall_Genus_" + curtime + ".tsv"
    with open(FileName, 'w') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        writer.writerow(Fields)
        writer.writerows(Overall)

    Fields = ["File Name", "PK" ,"Genus", "Number of Crossings","Seg Graph"]
    FileName = "Output_Files/"+ curtime + "/individual_Genus_" + curtime + ".tsv"
    with open(FileName, 'w') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        writer.writerow(Fields)
        writer.writerows(Individual)

if __name__ == '__main__':
    main()
