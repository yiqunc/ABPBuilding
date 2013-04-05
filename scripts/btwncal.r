library(igraph);
library(utils);
library(rgl);
library(Matrix);

setwd("/Users/yiqunc/UniProjects/ABPBuilding")

floorNames <- c("Ground", "FirstFloor", "SecondFloor", "ThirdFloor", "FourthFloor", "FifthFloor", "SixthFloor", "SeventhFloor")

floorShortNames <- c("ground", "1st", "2nd", "3rd", "4th", "5th", "6th", "7th")

floors <- list()

CONST_EDGE_WEIGHT_DEFAULT = 1.0
CONST_EDGE_WEIGHT_HUGE = 1000.0
CONST_FLOOR_TEXT_HEIGHT_RATIO = 2/(length(floorNames)-1)


# load .net network data
for(i in 1:length(floorNames)){
  floors[[i]] <- read.graph(file=sprintf("./data/OldABPBuilding/120601_enrichedAxial/0%i_%s/FABP_0%i.net",i-1,floorNames[i],i-1), format="pajek")
  V(floors[[i]])$name <-paste(sprintf("0%i_",i-1), as.integer(V(floors[[i]])$id)+1, sep="")
}

# build a bridge to connect adjacent floors
bridgeDF = read.table(file="./data/OldABPBuilding/120601_enrichedAxial/bridgeIdx.txt",header=TRUE,sep=",")
bridgeEdgeSeq = c()
bridgeVertexName = c()
bridgeVertexCounter = 1
for(j in 2:ncol(bridgeDF)){
  vtmp =c()
  for(i in 1:nrow(bridgeDF)){
    if(bridgeDF[i,j] > -1){
      vtmp = c(vtmp,bridgeVertexCounter)
      bridgeVertexName = c(bridgeVertexName, sprintf("0%i_%i",bridgeDF[i,1],bridgeDF[i,j]+1))
      #print(sprintf("[%i]:%i",bridgeVertexCounter,bridgeDF[i,j]))
      bridgeVertexCounter = bridgeVertexCounter +1
      vtmp = c(vtmp,bridgeVertexCounter)
    }
  }
  if (length(vtmp)>2){
    vtmp = vtmp[1:(length(vtmp)-2)]
  }
  #print(vtmp)
  bridgeEdgeSeq = c(bridgeEdgeSeq, vtmp)
}
gbridge = graph(bridgeEdgeSeq, directed=FALSE)
V(gbridge)$name = bridgeVertexName
#gbridge = graph(c(1,2,2,3,3,4,4,5,5,6,7,8,8,9,9,10,10,11,11,12),directed=FALSE)
#V(gbridge)$name = c("00_15","01_2","02_2","03_2","04_2","05_2", "00_14","01_1","02_1","03_1","04_1","05_1")

# init gunion
gunion = graph.union.by.name(floors[[1]], floors[[2]])

# merge floors
for (fidx in 3:length(floors)){
  gunion = graph.union.by.name(gunion, floors[[fidx]])
}
# merge bridge
gunion = graph.union.by.name(gunion, gbridge)

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

# back up new network before modification
gunion_original = gunion

# have a chance to modify the new network before caculation.
# brokenVerticeNames: a broken vertices vector, its connected edges index will be copied to brokenEdgesIdx.
# brokenEdgesSeq: a broken edge sequences vector, corresponding edge index will be copied to brokenEdgesIdx.
# brokenEdgesIdx: a broken edge index vector, which is used to create new network

brokenEdgesIdx = c()

# (1) define a broken vertices vector, all edges connected to these vertices will get a HUGE weight
brokenVerticeNames = c()
#if(length(brokenVerticeNames)>0){
#  for(vname in brokenVerticeNames){
    #E(gunion)[gunion[[vname,edges=TRUE]][[1]]]$weight = CONST_EDGE_WEIGHT_HUGE
#  }
#}
# (2) define a broken edges vector. all edges in this vector will get a HUGE weight
brokenEdgesSeq = c()
#for(i in 1:(length(brokenEdgesSeq)/2)){
#  edgeIdx = gunion[brokenEdgesSeq[(i-1)*2+1], brokenEdgesSeq[(i-1)*2+2], edges=TRUE]
#  if(edgeIdx > 0){
    #E(gunion)[edgeIdx]$weight = CONST_EDGE_WEIGHT_HUGE
#  }
#}

if(length(brokenEdgesSeq) > 0){
  for(i in 1:(length(brokenEdgesSeq)/2)){
    edgeIdx = gunion[brokenEdgesSeq[(i-1)*2+1], brokenEdgesSeq[(i-1)*2+2], edges=TRUE]
    if(edgeIdx > 0){
      brokenEdgesIdx = c(brokenEdgesIdx,edgeIdx)
    }
  }
}

if(length(brokenVerticeNames) > 0){
  for(vname in brokenVerticeNames){
    brokenEdgesIdx = c(brokenEdgesIdx, gunion[[vname,edges=TRUE]][[1]])
  }
}

brokenEdgesIdx = brokenEdgesIdx[!duplicated(brokenEdgesIdx)]

# (3) remove brokenEdges from network
#gunion = delete.edges(gunion,E(gunion, P=brokenEdges))
if(length(brokenEdgesIdx) > 0){
  gunion = delete.edges(gunion,E(gunion)[brokenEdgesIdx])
}

# (4) remove brokenVerticeNames from network
if(length(brokenVerticeNames) > 0){
  gunion = delete.vertices(gunion, brokenVerticeNames)
}

# calc the betweeness for vertex
V(gunion)$btwn=betweenness(gunion, directed=FALSE, normalized=TRUE)

# copy betweeness to the original network
for(vname in V(gunion)$name){
    if(vname %in% V(gunion_original)$name){
      V(gunion_original)[vname]$btwn = V(gunion)[vname]$btwn
  }
}

# save gunion as the final result
write.graph(gunion,file="./outputs/gunion.xml",format="graphml")

# use gunion_original for visualization test 
V(gunion_original)[btwn>0.1]$color = "red"
V(gunion_original)$size = 3
if(length(brokenVerticeNames) > 0){
  V(gunion_original)[brokenVerticeNames]$color = "black"
  V(gunion_original)[brokenVerticeNames]$size = 5
}

# have to set a default edge width value for all edges first 
E(gunion_original)$width = 1
if(length(brokenEdgesIdx) > 0){
  E(gunion_original)[brokenEdgesIdx]$color = "black"
  E(gunion_original)[brokenEdgesIdx]$width = 2 
}

# plot 2D network
#plot.igraph(gunion_original,vertex.color=V(gunion_original)$color,vertex.label=NA,vertex.size=3)


mat = do.call(cbind, list(V(gunion_original)$x, V(gunion_original)$y, V(gunion_original)$z))
# plot 3D network
#rglplot(gunion_original,vertex.color=V(gunion_original)$color,vertex.label=NA,vertex.size=3,layout=mat)
rglplot(gunion_original,vertex.label=NA,layout=mat)
rgl.viewpoint(0,-90);

# to mathc with floor position, the z value for text has to be clamped to [-1,1]. The interval between floor text is CONST_FLOOR_TEXT_HEIGHT_RATIO
for(i in 1:length(floorShortNames)){
  rgl.texts(x=-1.2, y=-1.2, z=-1+(i-1)*CONST_FLOOR_TEXT_HEIGHT_RATIO, text=floorShortNames[i])
}

#start <- proc.time()[3]
#while ((i <- 36*(proc.time()[3]-start)) < 3360) {
#  rgl.viewpoint(i/20,i/30); 
#}


