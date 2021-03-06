### Edgelist gives case ids of citing and cited cases
edgelist <- read.table("allcites.txt",header=F)
## Need to add characters to ids so that they are not used as integers by network functions
edgelist[,1] <- paste("c",edgelist[,1],sep="")
edgelist[,2] <- paste("c",edgelist[,2],sep="")

### Read in node attributes
node_attributes <- read.csv("judicial.csv",stringsAsFactors=F)
node_attributes[,1] <- paste("c",node_attributes$caseid,sep="")

### Add additional attributes from the SCDB
## See data documentation, http://scdb.wustl.edu/_brickFiles/2016_01/SCDB_2016_01_codebook.pdf
SCDB <- read.csv("SCDB_2016_01_caseCentered_Citation.csv",stringsAsFactors=F)
SCDBinFowler <- match(SCDB$usCite,node_attributes$usid)
SCDB <- SCDB[!is.na(SCDBinFowler),]
SCDBinFowler <- match(SCDB$usCite,node_attributes$usid)
node_attributes$issueArea <- NA
node_attributes$issueArea[SCDBinFowler] <- SCDB$issueArea
node_attributes$lawType <- NA
node_attributes$lawType[SCDBinFowler] <- SCDB$lawType

### Create a netowrk object
## The following sequence assures that edges are matched to the correct nodes
library(network)
## this is a directed acyclic graph, so it makes most senese to analyze it as undirected (i.e., since there can be no cycling)
scnetwork1 <- network.initialize(nrow(node_attributes),dir=F)
## set name to case ids, which will make it easy to add edges
network.vertex.names(scnetwork) <- node_attributes[,1]

# convert edgelist to a matrix
edgelist <- as.matrix(edgelist)

# match node indices to numeric place
edge_row <- match(edgelist[,1],network.vertex.names(scnetwork))
edge_col <- match(edgelist[,2],network.vertex.names(scnetwork))

# estimate timing
test_time <- system.time( for(i in 2:200){ scnetwork[cbind(edge_row[(i-1):i], edge_col[(i-1):i])] <- 1} )

# convert to hours
estimated_hours <- test_time/199*length(edge_row)/60/60

# view estimated timing
print(estimated_hours)

# run all rows to store edges
for(i in 98500:length(edge_row)){ #stopped at 98500
  if(i %% 250 == 0) cat("Starting iteration", i, "\n")
  scnetwork[cbind(edge_row[(i-1):i], edge_col[(i-1):i])] <- 1}

### Lets use the age of the case and the salience indicator
set.vertex.attribute(scnetwork,c("year","salience", "area", "type"),node_attributes[,c("year","oxford", "issueArea", "lawType")])


# only use vertices created after 1953
sc<- scnetwork

scnetwork54<- delete.vertices(sc, 1:21065)
rm(sc)

save.image(file="supreme54.RData")

# only use vertices where issua area and law type is known

sc<- scnetwork

nodes.to.delete <- which(is.na(node_attributes$issueArea) | is.na(node_attributes$lawType), arr.ind=TRUE)

scnetwork_additional_data <- delete.vertices(sc, nodes.to.delete)
rm(sc)


save.image(file="supreme_additional_data.RData")
