# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Description: this code snippet is developed to construct a multi-floor space network and then analyze, visualize the betweenness of the vertices. 
# Version: 0.1
# Last Updated: 6-April-2013
# Supervised by: Martin Tomko     tomkom@unimelb.edu.au
# Created by: Yiqun(Benny) Chen     chen1@unimelb.edu.au
#
# Read Me:
# (1) The code is only tested on Mac (10.8.2).
# (2) rlg pacakge is requred, which needs XQuartz library for 3D visualization.
# (3) change the working directory properly to run on your machine.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# required packages

library(igraph);
library(utils);
library(rgl);
library(Matrix);

# setup R working directory 
setwd("/Users/yiqunc/UniProjects/ABPBuilding")

# set data directory path
dataPath <- "./data/OldABPBuilding/120601_enrichedAxial"

# floor fullname vector, for loading individual floor plan (.net)
floorNames <- c("Ground", "FirstFloor", "SecondFloor", "ThirdFloor", "FourthFloor", "FifthFloor", "SixthFloor", "SeventhFloor")

# floor shortname vector, for displaying in 3D scene. 
floorShortNames <- c("ground", "1st", "2nd", "3rd", "4th", "5th", "6th", "7th")

# a list to contain graph objects for each floor
floors <- list()

# constants declaration
# edge default weight
CONST_EDGE_WEIGHT_DEFAULT = 1.0
# edge penalty weight
CONST_EDGE_WEIGHT_HUGE = 1000.0
# floor text height interval in 3D scene
CONST_FLOOR_TEXT_HEIGHT_RATIO = 2/(length(floorNames)-1)
# render colorramp size
CONST_COLRAMP_SIZE = 5
CONST_COLRAMP_2ENDS = c("white","red")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Step 1: construct a multi-floor network based on individual floor plans
# load individual .net network data, all network will be unioned by vertex name. assign a proper name for every vertex.
for(i in 1:length(floorNames)){
  floors[[i]] <- read.graph(file=sprintf("%s/0%i_%s/FABP_0%i.net", dataPath, i-1,floorNames[i],i-1), format="pajek")
  V(floors[[i]])$name <-paste(sprintf("0%i_",i-1), as.integer(V(floors[[i]])$id)+1, sep="")
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Step 2: Build a bridge to connect floors
# In the original ".net" files, there is no connection between each floor, A "bridgeIdx.txt" file is built to link them all together, the structure is:
# FN,SA1,SA2,L1,L2,SB1,SB2
# "00",11,12,14,13,5,3
# "01",21,20,1,0,84,83
# "02",8,18,1,0,32,-1
# "03",4,6,1,0,30,-1
# "04",6,4,1,0,13,-1
# "05",40,38,1,0,-1,-1
# "06",1,6,-1,-1,-1,-1
# "07",2,6,-1,-1,-1,-1
# Where FN stands for "floor name"; SA1,SA2,SB1,SB2 are four marked stairs on the "BuildingBinder.pdf" file: L1, L2 are lifts. They are all joints. 
# The number refers to the original vertex index in each floor's ".net" file: 
# E.g, 11 stands for the "SA1" vertex index in "00" .net file. 
# -1 means the joint does not exist on that floor. E.g. SB2 is not shown on "02" floor.
# Ideally, existing joints in the same position of each floor should connect vertically. E.g. for SA1, 11-21-8-4-6-40-1-2 forms a continous edge; for SB2, 3-83 forms a continous edge.
# The bridge graph is created for the ideal scenario, and When exception exists, the "brokenedge" will handle it.

# load the bridge data.
bridgeDF <- read.table(file=sprintf("%s/bridgeIdx.txt", dataPath),header=TRUE,sep=",")
bridgeEdgeSeq <- c()
bridgeVertexName <- c()
bridgeVertexCounter = 1
for(j in 2:ncol(bridgeDF)){
  vtmp <- c()
  for(i in 1:nrow(bridgeDF)){
    if(bridgeDF[i,j] > -1){
      vtmp <- c(vtmp,bridgeVertexCounter)
      bridgeVertexName = c(bridgeVertexName, sprintf("0%i_%i",bridgeDF[i,1],bridgeDF[i,j]+1))
      bridgeVertexCounter = bridgeVertexCounter +1
      vtmp <- c(vtmp,bridgeVertexCounter)
    }
  }
  if (length(vtmp)>2){
    vtmp <- vtmp[1:(length(vtmp)-2)]
  }
  bridgeEdgeSeq <- c(bridgeEdgeSeq, vtmp)
}

# create a graph for the bridge
gbridge <- graph(bridgeEdgeSeq, directed=FALSE)
V(gbridge)$name <- bridgeVertexName

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Step 3: create the entire multi-floor network
# init gunion
gunion <- graph.union.by.name(floors[[1]], floors[[2]])
# merge floors
for (fidx in 3:length(floors)){
  gunion <- graph.union.by.name(gunion, floors[[fidx]])
}
# merge with bridge
gunion <- graph.union.by.name(gunion, gbridge)

# copy coords (x,y,z) from original graph based on name
for(i in 1:length(V(gunion)$name)){
  vname = V(gunion)$name[i]
  for(fidx in 1:length(floors)){
    if(vname %in% V(floors[[fidx]])$name){
      V(gunion)[i]$x = V(floors[[fidx]])[V(floors[[fidx]])$name==vname]$x
      V(gunion)[i]$y = V(floors[[fidx]])[V(floors[[fidx]])$name==vname]$y
      V(gunion)[i]$z = fidx - 1
      break;
    }
  }
}

# init betweenness value and vetex color
V(gunion)$btwn = -1

# init edge weight
E(gunion)$weight = CONST_EDGE_WEIGHT_DEFAULT

# load spaceinfo database
spaceInfoDF <- read.csv(file=sprintf("%s/spaceinfo.txt", dataPath),header=TRUE,sep=",")

# create a unique space index
spaceInfoDF[,"UNISIDX"]=sprintf("0%i_%i",spaceInfoDF[,"FN"], spaceInfoDF[,"SIDX"]+1)

# attach space info to graph
for(i in 1:length(V(gunion)$name)){
  vname = V(gunion)$name[i]
  if (vname %in% spaceInfoDF[,"UNISIDX"]){
    filter = vname==spaceInfoDF[,"UNISIDX"]
    V(gunion)[i]$spacename = as.character(spaceInfoDF[filter,"SN"])
    V(gunion)[i]$capacity = spaceInfoDF[filter,"CAPACITY"]
    V(gunion)[i]$accesslevel = spaceInfoDF[filter,"ACCESSLVL"]
    V(gunion)[i]$owners = as.character(spaceInfoDF[filter,"OWNER"])
  }
}

# back up new network before modification
gunion_original <- gunion

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Step 4: use brokenedge and brokenvertice to modify network.

brokenVerticesDF <- read.table(file=sprintf("%s/brokenvertices.txt", dataPath),header=TRUE,sep=",")
brokenEdgesDF <- read.table(file=sprintf("%s/brokenedges.txt", dataPath),header=TRUE,sep=",")

brokenEdgesIdx <- c()

# find broken edge indices
if(nrow(brokenEdgesDF) > 0){
  # check if vertex1 names are valid in network
  filterV1 = brokenEdgesDF[,1] %in% V(gunion)$name
  # check if vertex2 names are valid in network
  filterV2 = brokenEdgesDF[,2] %in% V(gunion)$name
  # v1 v2 both must be valid
  filter = filterV1 & filterV2
  # remove those invalid ones
  brokenEdgesDF <- subset(brokenEdgesDF,filter)
  
  if(nrow(brokenEdgesDF) > 0){
    for(i in 1:nrow(brokenEdgesDF)){
      edgeIdx = gunion[as.character(brokenEdgesDF[i,1]), as.character(brokenEdgesDF[i,2]), edges=TRUE]
      # check if broken edge exists in current network
      if(edgeIdx > 0){
        brokenEdgesIdx <- c(brokenEdgesIdx, edgeIdx)
      }
    }
  }
}

# find all edge indices which are connected to broken vertices
if(nrow(brokenVerticesDF) > 0){
  # test if broken vertices exist in network
  filter = brokenVerticesDF[,1] %in% V(gunion)$name
  # remove those invalid ones
  brokenVerticesDF <- subset(brokenVerticesDF,filter)
  # then find connected edges to remove
  if(nrow(brokenVerticesDF) > 0){
    for(i in 1:nrow(brokenVerticesDF)){
      brokenEdgesIdx <- c(brokenEdgesIdx, gunion[[as.character(brokenVerticesDF[i,1]),edges=TRUE]][[1]])
    }
  }
}

# remove duplicated edge indice before remove them
brokenEdgesIdx <- brokenEdgesIdx[!duplicated(brokenEdgesIdx)]

# remove brokenEdges from network
if(length(brokenEdgesIdx) > 0){
  gunion <- delete.edges(gunion,E(gunion)[brokenEdgesIdx])
}

# remove brokenVertices from network
if(nrow(brokenVerticesDF) > 0){
  gunion <- delete.vertices(gunion, as.character(brokenVerticesDF[,1]))
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Step 4: calc the betweeness for vertex
V(gunion)$btwn <- betweenness(gunion, directed=FALSE, normalized=TRUE)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Step 5: save gunion as the final result
write.graph(gunion,file="./outputs/gunion.xml",format="graphml")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Step 6: for visualization
# copy betweeness to the original network
for(vname in V(gunion)$name){
    if(vname %in% V(gunion_original)$name){
      V(gunion_original)[vname]$btwn = V(gunion)[vname]$btwn
  }
}

# use gunion_original for visualization test
colorramp = colorRampPalette(CONST_COLRAMP_2ENDS)(CONST_COLRAMP_SIZE)
vecBtwn = V(gunion_original)$btwn
itvValue = (max(vecBtwn)-min(vecBtwn)) / (CONST_COLRAMP_SIZE-1)
V(gunion_original)$colorIdx = as.integer(V(gunion_original)$btwn / itvValue) + 1
V(gunion_original)$color = colorramp[V(gunion_original)$colorIdx]

#V(gunion_original)[btwn>0.1]$color = "red"
V(gunion_original)$size = 3
if(nrow(brokenVerticesDF) > 0){
  V(gunion_original)[as.character(brokenVerticesDF[,1])]$color = "black"
  V(gunion_original)[as.character(brokenVerticesDF[,1])]$size = 5
}

# get a middle color between to colors : colorRampPalette(c("#FFFFFF","#FF0000"))(3)[2]
# the edge color is calculated based on its two vertices color
edgeVerticePair = get.edges(gunion_original,E(gunion_original))
for(i in 1:nrow(edgeVerticePair)){
  sColor = V(gunion_original)[edgeVerticePair[i,1]]$color
  eColor = V(gunion_original)[edgeVerticePair[i,2]]$color
  E(gunion_original)[i]$color = colorRampPalette(c(sColor,eColor))(3)[2]
}

# have to set a default edge width value for all edges first 
E(gunion_original)$width = 1
if(length(brokenEdgesIdx) > 0){
  E(gunion_original)[brokenEdgesIdx]$color = "black"
  E(gunion_original)[brokenEdgesIdx]$width = 2 
}

# plot 2D network
#plot.igraph(gunion_original,vertex.color=V(gunion_original)$color,vertex.label=NA,vertex.size=3)

mat <- do.call(cbind, list(V(gunion_original)$x, V(gunion_original)$y, V(gunion_original)$z))
# plot 3D network
rglplot(gunion_original,vertex.label=NA,layout=mat)
rgl.viewpoint(0,-90);

# to mathc with floor position, the z value for text has to be clamped to [-1,1]. The interval between floor text is CONST_FLOOR_TEXT_HEIGHT_RATIO
for(i in 1:length(floorShortNames)){
  rgl.texts(x=-1.2, y=-1.2, z=-1+(i-1)*CONST_FLOOR_TEXT_HEIGHT_RATIO, text=floorShortNames[i])
}

# auto rotate the 3D network for video recording
#start <- proc.time()[3]
#while ((i <- 36*(proc.time()[3]-start)) < 3360) {
#  rgl.viewpoint(i/20,i/30); 
#}



