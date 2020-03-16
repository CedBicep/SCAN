#-- Program SCAN : compares differents assemblies method using similarity networks --#
options(warn = -1)
library(igraph)
library(yaml)
library(RSvgDevice)
require(plotrix)
library(snow)

#############################################read the yaml file to upload the required parameters#################################################################

args <- (commandArgs(TRUE))
data <- yaml.load_file(args[1])
n <- data$nbars
alpha <- data$alpha
method <- attributes(data$metas)$name
tagref <- ""
tagtest <- ""
centralities <- c("Porportion","ArticulationLocal","ArticulationGlobal","DegreeOne","Geodesic")

############################################################declaration of functions##############################################################################

#----------------------------------------------------------------#
#- Function HasTag : tells nodes having a specific tag		-#
#----------------------------------------------------------------#

#-----------------------------------
#- cc: graph object
#- tag: a given string of characters
#-----------------------------------

HasTag <- function(cc,tag)
{
	has.this.tag <- sapply(strsplit(V(cc)$name, "_"), "[",2)==tag
	has.this.tag
}


#------------------------------------------------------------------------------------------------------------------#
#- Function GetNumLocalArtPoint : gives the number of specific nodes that are local articulation point in a graph -#
#------------------------------------------------------------------------------------------------------------------#

#-----------------------------------
#- cc: graph object
#- tag: a given string of characters
#-----------------------------------

GetNumLocalArtPoint <- function(cc,tag)
{
	all <- which(HasTag(cc,tag)==TRUE)
	ret <- 0
	if(length(all) > 0)
	{
		g <- graph.neighborhood(cc, 1, all)
		ret <- sum(unlist(sapply(g,Is.Articulationpoint)) & unlist(sapply(g,HasTag,tag = tag)))
	}
	ret
}


#------------------------------------------------------------------------------------------------#
#- Function Is.Articulationpoint : tells which node are articulation point in a graph		-#
#------------------------------------------------------------------------------------------------#

#-----------------------------------
#- cc: graph object
#-----------------------------------

Is.Articulationpoint <- function(cc)
{
	is.articulation <- logical(length(V(cc)))
	is.articulation[articulation.points(cc)] <- TRUE
	is.articulation
}


#------------------------------------------------------------------------#
#- Function Distance : computes the jaccard index of two samples	-#
#------------------------------------------------------------------------#

#-----------------------------------
#- x: numeric vector
#- y: numeric vector
#-----------------------------------

Distance <- function(x,y)
{
	ret <- length(intersect(x,y))/length(union(x,y))
	ret
}


#--------------------------------------------------------------------------------#
#- Function PlotTwoDistribution : plotting function of two distributions	-#
#--------------------------------------------------------------------------------#

#-----------------------------------
#- dist1: numeric vector
#- dist2: numeric vector
#- tag1: a given string of characters
#- tag2: a given string of characters
#- bar: numeric
#- sim: numeric
#- meth: a given string of characters
#- xlab: a given string of characters
#- ylab: a given string of characters
#- main: a given string of characters
#-----------------------------------

PlotTwoDistribution <- function(dist1,dist2,tag1,tag2,bar,sim,meth,xlab,ylab,main)
{	
	if(length(unique(c(dist1,dist2))) > 1)
	{
		cols <- c("azure3","black")
		l <- list(dist1,dist2)
		mh <- multhist(l,breaks = bar,args.legend = list(x = "topleft"),plot.it = FALSE)
		devSVG(file = file.path("images",paste(sim,tag1,tag2,meth,bar,"plot.svg", sep = "_"), fsep = .Platform$file.sep), bg = "white",width = 7,height = 7)
		multhist(l,ylim = c(0,max(mh[[2]])*1.4), xlab = xlab,ylab = ylab,main = main, breaks = bar)
		legend("topleft",legend = c(tag1, tag2),lwd = 4,bty = "n",col = cols,y.intersp = 2.5,inset = 0.095)
		dev.off()
	}
}


#--------------------------------------------------------------------------------#
#- Function PlotOneDistribution : plotting function of one distribution 	-#
#--------------------------------------------------------------------------------#

#-----------------------------------
#- dist: numeric vector
#- tag1: a given string of characters
#- tag2: a given string of characters
#- bar: numeric
#- sim: numeric
#- meth: a given string of characters
#- xlab: a given string of characters
#- ylab: a given string of characters
#- main: a given string of characters
#-----------------------------------

PlotOneDistribution <- function(dist,tag1,tag2,bar,sim,meth,xlab,ylab,main)
{
	#cols <- c("azure3","black")
	mh <- hist(dist,nclass = bar,plot = FALSE)
	devSVG(file = file.path("images",paste(sim,tag1,tag2,meth,bar,"plot.svg", sep = "_"), fsep = .Platform$file.sep), bg = "white",width = 7,height = 7)
	hist(dist,ylim = c(0,max(mh$count)*1.4),col = c("azure3"),xlab = xlab,ylab = ylab,main = main,nclass = bar)
	dev.off()
}


#--------------------------------------------------------------------------------#
#- Function WriteCC : writes a connected component in a cytoscape file format	-#
#--------------------------------------------------------------------------------#

#-----------------------------------
#- cc: graph object
#- filename: a given string of characters
#- sim: a given string of characters
#- tag: a given string of characters
#-----------------------------------

WriteCC <- function(cc,filename,sim,tag)
{
	x <- data.frame(get.edgelist(cc, names = TRUE))	
	write.table(x, file = file.path("data",paste(sim,"_",tag,"_",filename,".txt", sep = ""), fsep = .Platform$file.sep), sep = "\t",row.names = FALSE,col.names = FALSE, quote = FALSE, qmethod = "escape",append = TRUE)		
}


#----------------------------------------------------------------#
#- Function WriteSeq : writes contigs in text file		-#
#----------------------------------------------------------------#

#-----------------------------------
#- cc: graph object
#- filename: a given string of characters
#- sim: a given string of characters
#- tag: a given string of characters
#-----------------------------------

WriteSeq <- function(cc,filename,sim,tag)
{	
	prefixe <- sapply(strsplit(V(cc)$name, "_"), "[",2)
	contig <- V( cc)$name[which(prefixe==tag)]
	x <- data.frame(contig)
	write.table(x, file = file.path("data",paste(sim,"_",tag,"_",filename,".txt", sep = ""), fsep = .Platform$file.sep), sep = "\t",row.names = FALSE,col.names = FALSE, quote = FALSE, qmethod = "escape",append = TRUE)
}


#----------------------------------------------------------------#
#- Function SumUpResults : writes a summary of all centralities	-#
#----------------------------------------------------------------#

#-----------------------------------
#- nParam: numeric
#- nbrContig: numeric
#- totalContig: numeric
#- tag: a given string of characters
#- sim: a given string of characters
#-----------------------------------

SumUpResults <- function(nParam,nbrContig,totalContig,tag,sim)
{
	FileConn <- file(file.path("data",paste(tag,sim,"Summary.txt",sep = "_"),fsep = .Platform$file.sep),"w")
	if(nParam==0)
	{
		writeLines(paste("No suspicious assembly. ",nbrContig," contig(s) seem successful",sep = ""), FileConn)
	}
	else
	{
		writeLines(paste("Suspicious assembly according to ",nParam," parameter(s).",nbrContig," contig(s) ",round((nbrContig/totalContig)*100),"% seem successful",sep = ""), FileConn)
	}
	close(FileConn)
}


#------------------------------------------------------------------------#
#- Function MakeGoodContigList : writes sequences name from graph 	-#
#------------------------------------------------------------------------#

#-----------------------------------
#- filename: a given string of characters
#- tag: a given string of characters
#-----------------------------------

MakeGoodContigList <- function(contigs,filename,tag)
{	
	write.table(data.frame(contigs), file = paste("GoodContigList_",filename, sep = ""), sep = "\t",row.names = FALSE,col.names = FALSE, quote = FALSE, qmethod = "escape",append = FALSE)
}


#----------------------------------------------------------------------------------------#
#- Function MakeSummaryOfThisMethod : writes the summary for one method	 		-#
#----------------------------------------------------------------------------------------#

#-----------------------------------
#- method: a given string of characters
#- centrality: a given vector of string of characters
#- similarity: numeric
#- pvalue: a given numeric vector
#- distance: a given numeric vector
#-----------------------------------

MakeSummaryOfThisMethod <- function(method,reference,centrality,similarity,pvalue,distance)
{
	FileConn <- file("SummaryByMethod.txt","a")
	writeLines(paste(method,reference,centrality,similarity,unlist(pvalue),unlist(distance),sep = "\t"), FileConn)
	close(FileConn)
}


#----------------------------------------------------------------------------------------#
#- Function MakeSummaryOfCentrality : writes the summary for one centrality		-#
#----------------------------------------------------------------------------------------#

#-----------------------------------
#- fic
#- method: a given string of characters
#- centrality: a given vector of string of characters
#- similarity: numeric
#- pvalue: a given numeric vector
#- distance: a given numeric vector
#-----------------------------------

MakeSummaryOfCentrality <- function(fic,tag1,tag2,method,centrality,similarity,pvalue,distance)
{
	FileConn<-file(file.path("data",paste(similarity,method,tag1,tag2,fic,".txt",sep = "_"), fsep = .Platform$file.sep),"w")
	writeLines(paste("Method","Centrality","Network Similarity","pvalue","Distance",sep = "\t"), FileConn)
	writeLines(paste(method,centrality,similarity,pvalue,distance,sep = "\t"), FileConn)
	close(FileConn)
}

############################################################################################################################################

#----------------------------------------------------------------------------------------#
#- Function GetProportion : computes the proportion of different nodes in a graph	-#
#----------------------------------------------------------------------------------------#

#----------------------------------- 
#- cc: graph object
#- tag1: a given string of characters
#- tag2: a given string of characters
#-----------------------------------
GetProportion <- function(cc,tag1,tag2)
{
	vc <- vcount(cc)
	ntag1 <- sum(HasTag(cc,tag1))
	ntag2 <- sum(HasTag(cc,tag2))
	if(sum(ntag1,ntag2)==0)
	{
		return(c(1,ntag1/vc,ntag2/vc))
	}
	return(c(prop.test(c(ntag1,ntag2),c(vc,vc))$p.value,ntag1/vc,ntag2/vc))
}


#--------------------------------------------------------------------------------#
#- Function GetlonguestGeod : computes the longuest geodesic in a graph		-#
#--------------------------------------------------------------------------------#

#----------------------------------- 
#- cc: graph object
#- tag1: a given string of characters
#- tag2: a given string of characters
#-----------------------------------

GetlonguestGeod <- function(cc,tag1,tag2) 
{
	tg2 <- delete.vertices(cc, which(sapply(strsplit(V(cc)$name, "_"), "[",2)!=tag2))
	tg1 <- delete.vertices(cc, which(sapply(strsplit(V(cc)$name, "_"), "[",2)!=tag1))
	return(c( diameter(tg1, weights = rep(1, ecount(tg1))),diameter(tg2, weights = rep(1, ecount(tg2))) ))
}


#------------------------------------------------------------------------------------------------#
#- Function GetPropLocalArtPoint : computes the proportion of local articulation point		-#
#------------------------------------------------------------------------------------------------#

#----------------------------------- 
#- cc: graph object
#- tag1: a given string of characters
#- tag2: a given string of characters
#-----------------------------------
GetPropLocalArtPoint <- function(cc,tag1,tag2)
{
	vc <- vcount(cc)
	ntag1 <- GetNumLocalArtPoint(cc,tag1)
	ntag2 <- GetNumLocalArtPoint(cc,tag2)
	if(sum(ntag1,ntag2)==0)
	{
		return(c(1,ntag1/vc,ntag2/vc))
	}
	return(c(prop.test(c(ntag1,ntag2),c(vc,vc))$p.value,ntag1/vc,ntag2/vc))	
}


#------------------------------------------------------------------------------------------------#
#- Function GetPropGlobalArtPoint : computes the proportion of global articulation point	-#
#------------------------------------------------------------------------------------------------#

#----------------------------------- 
#- cc: graph object
#- tag1: a given string of characters
#- tag2: a given string of characters
#-----------------------------------

GetPropGlobalArtPoint <- function(cc,tag1,tag2)
{
	vc <- vcount(cc)
	ntag1 <- sum((Is.Articulationpoint(cc)) & (HasTag(cc,tag1)))
	ntag2 <- sum((Is.Articulationpoint(cc)) & (HasTag(cc,tag2)))
	if(sum(ntag1,ntag2)==0)
	{
		return(c(1,ntag1/vc,ntag2/vc))
	}
	return(c(prop.test(c(ntag1,ntag2),c(vc,vc))$p.value,ntag1/vc,ntag2/vc))	
}


#--------------------------------------------------------------------------------#
#- Function GetPropDegree1 : computes the proportion of node with degree=1 	-#
#--------------------------------------------------------------------------------#

#----------------------------------- 
#- cc: graph object
#- tag1: a given string of characters
#- tag2: a given string of characters
#-----------------------------------

GetPropDegree1 <-function(cc,tag1,tag2)
{
	vc <- vcount(cc)
	ntag1 <- sum((degree(cc)==1) & (HasTag(cc,tag1)))
	ntag2 <- sum((degree(cc)==1) & (HasTag(cc,tag2)))
	if(sum(ntag1,ntag2)==0)
	{
		return(c(1,ntag1/vc,ntag2/vc))
	}
	return(c(prop.test(c(ntag1,ntag2),c(vc,vc))$p.value,ntag1/vc,ntag2/vc))	
}


#------------------------------------------------------------------------------------------------#
#- Function HasGoodJaccard : tells if a two tagged nodes have a high jaccard in a graph 	-#
#------------------------------------------------------------------------------------------------#

#-----------------------------------
#- cc: graph object
#- tag1: a given string of characters
#- tag2: a given string of characters
#-----------------------------------

HasGoodJaccard <- function(cc,tag1,tag2)
{
	adj <- (get.adjlist(cc))
	gel <- get.edgelist(cc,names = FALSE)
	tgRef <- which(HasTag(cc,tag = tag1))
	tgTest <- which(HasTag(cc,tag = tag2))
	AllDist <- numeric(1)
	ret <- list()
	ret[["bool"]] <- FALSE
	ret[["dist"]] <- 0
	if( (vcount(cc)==2) & (length(tgRef) > 0) & (length(tgTest) > 0) )
	{
		ret[["bool"]] <- TRUE
		ret[["dist"]] <- 1
		return(ret)
	}
	lns <- which( ( ((gel[,2] %in% tgTest) & (gel[,1] %in% tgRef)) | ((gel[,1] %in% tgTest) & (gel[,2] %in% tgRef)) ) == TRUE)
	if(length(lns) > 0)
	{
		AllDist <- numeric(length(lns))
		for(i in 1:length(lns))
		{
			from <- as.integer(unlist(adj[ gel[lns[i],2] ]))
			to <- as.integer(unlist(adj[ gel[lns[i],1] ]))
			from <- from[from != gel[lns[i],2]]
			from <- from[from != gel[lns[i],1]]
			to <- to[to != gel[lns[i],2]]
			to <- to[to != gel[lns[i],1]]
			dst <- Distance(from,to)
			AllDist[i] <- dst
			if(dst > 0.9)
			{
				ret[["bool"]] <- TRUE
			}
		}
		ret[["dist"]] <- AllDist
	}
	return(ret)
}


#------------------------------------------------------------------------------------------------#
#- Function decomposeGraph : finds connected components in a graph 				-#
#------------------------------------------------------------------------------------------------#

#-----------------------------------
#- file: a given string of characters
#-----------------------------------

decomposeGraph<-function(file)
{
	sif <- read.csv(file, head = FALSE, sep = "\t")
	from <- as.character(sif[,1])
	to <- as.character(sif[,2])
	
	graph <- graph.data.frame(data.frame(from = from,to = to), directed = FALSE)	#overall grpah
	CC <- decompose.graph(graph, min.vertices = 2)				#find all connected components
	return(CC)
}

############################################################################################################################################

#----------------------------------------#
#- Function Main : Start of the script	-#
#----------------------------------------#

Main <- function()
{
	print("Loading input files ...",quote=FALSE)
	unlink("images", recursive = TRUE)
	unlink("data", recursive = TRUE)
	
	dir.create("images", showWarnings = TRUE, recursive = FALSE, mode = "0777")
	dir.create("data", showWarnings = TRUE, recursive = FALSE, mode = "0777")
	
	FileConn <- file("SummaryByMethod.txt","w")
	writeLines(paste("Method","Reference","Centrality","NetworkSimilarity","pvalue","Distance",sep="\t"), FileConn)
	close(FileConn)	
	
	
	description=vector()
	for(l in 1:length(data$test))
	{
		for(m in 1:length(data$ref[[l]]))
		{
			#print(paste(data$ref[[l]][m],"_",data$test[l],sep=""))
			description=append(description,paste(data$test[l],"_",data$ref[[l]][m],sep=""))
		}	
	}
	
	
	FileConn <- file("CentralitySummary.txt","w")
	writeLines(paste(c("Similarity",description),collapse="\t"), FileConn)#attributes(data$metas)$names)
	close(FileConn)
	
	FileConn <- file("NetworksSummary.txt","w")
	writeLines(paste(c("Method","Reference","Similarity","Network having the most good contigs","# of Good contigs"),collapse="\t"), FileConn)
	close(FileConn)
	
	similarities <- attributes(data$metas[[1]])$names
	centralitiesSummary <- data.frame(similarities)
	
	print("SCAN is computing, please wait ...",quote=FALSE)
	
	for(i in 1:length(data$metas))	#foreach method
	{
		
		tagtest <- data$test[i]
		
		AllPvalForThisMeth <- list()	#liste of list of pvalue for each similarity and each centrality
		AllDistForThisMeth <- list()	#liste of list of distance for each similarity and each centrality
		numGoodPvalForThisMeth <- integer(length(data$metas[[i]]))		#number of good pvalue according the similarity	
		numBadPvalForThisMeth <- integer(length(data$metas[[i]]))
		numGoodContigForThisMeth <- integer(length(data$metas[[i]]))
		GoodContigForThisMeth <- list()
		
		
		cl <- makeSOCKcluster(rep("localhost",data$nproc),type = "SOCK")
		clusterEvalQ(cl, library(igraph))
		clusterExport(cl,"decomposeGraph")
		AllCC<-parSapply(cl, data$metas[[i]],decomposeGraph,simplify = FALSE )
		stopCluster(cl)

		for(k in 1:length(data$ref[[i]]))
		{
			tagref <- data$ref[[i]][k]
				
			for(j in 1:length(data$metas[[i]]))	#for each similarity
			{				
				CC <- AllCC[[j]][]
				CC <- CC[ which((sapply(sapply(CC,HasTag,tag=tagref),sum) > 0) & (sapply(sapply(CC,HasTag,tag=tagtest),sum) > 0)) ]

				numOfContig <- sum(unlist(sapply(CC,HasTag,tag=tagtest)))
				
				cl <- makeSOCKcluster(rep("localhost",data$nproc),type = "SOCK")
				
				clusterEvalQ(cl, library(igraph))
				clusterExport(cl,"tagref")
				clusterExport(cl,"tagtest")
				clusterExport(cl,"Distance")
				clusterExport(cl,"HasTag")
				clusterExport(cl,"GetNumLocalArtPoint")
				clusterExport(cl,"Is.Articulationpoint")
				
				jaccard <- parSapply(cl,CC,HasGoodJaccard,tag1 = tagref,tag2 = tagtest)		#
				geodesics <- parSapply(cl,CC,GetlonguestGeod,tag1 = tagref,tag2 = tagtest)		#
				degree1 <- parSapply(cl,CC,GetPropDegree1,tag1 = tagref,tag2 = tagtest)		#
				artpointglobal <- parSapply(cl,CC,GetPropGlobalArtPoint,tag1 = tagref,tag2 = tagtest)	#
				artpointlocal <- parSapply(cl,CC,GetPropLocalArtPoint,tag1 = tagref,tag2 = tagtest)	#
				proportion <- parSapply(cl,CC,GetProportion,tag1 = tagref,tag2 = tagtest)		#	
				
				stopCluster(cl)
				
				ks1 <- ks.test(proportion[2,],proportion[3,])				
				ks2 <- ks.test(artpointlocal[2,],artpointlocal[3,])
				ks3 <- ks.test(artpointglobal[2,],artpointglobal[3,])
				ks4 <- ks.test(degree1[2,],degree1[3,])
				ks5 <- ks.test(geodesics[2,],geodesics[1,])
				
				Allresult <- matrix(rep(0,length(CC)*6), nrow = 6, ncol = length(CC))	
				Allresult[1,which(unlist(jaccard[1,])==TRUE)] <- 1
				if(ks5$p.value > alpha )						#a special case
				{
					Allresult[2,] <- 1
				}
				Allresult[3,which(degree1[1,] > alpha)] <- 1
				Allresult[4,which(artpointglobal[1,] > alpha)] <- 1
				Allresult[5,which(artpointlocal[1,] > alpha)] <- 1
				Allresult[6,which(proportion[1,] > alpha)] <- 1
				goodCC <- which(colSums(Allresult) > as.numeric(data$centhr))
				s <- sapply(CC[goodCC],WriteSeq,sim = names(data$metas[[i]][j]),filename = "SafeContigs4all",tag = tagtest)
				
				numGoodContigForThisMeth[j] <- sum(unlist(sapply(CC[goodCC], HasTag,tag = tagtest)))
				GoodContigForThisMeth[[j]] <- unlist(sapply(CC[goodCC],get.vertex.attribute,name="name"))[unlist(sapply(CC[goodCC],HasTag,tag=tagtest))]
				
				numGoodPvalForThisMeth[j] <- sum(c(ks1$p.value,ks2$p.value,ks3$p.value,ks4$p.value,ks5$p.value) > alpha)
				numBadPvalForThisMeth[j] <- sum(c(ks1$p.value,ks2$p.value,ks3$p.value,ks4$p.value,ks5$p.value) <= alpha)
				AllPvalForThisMeth[[j]] <- c(ks1$p.value,ks2$p.value,ks3$p.value,ks4$p.value,ks5$p.value)
				AllDistForThisMeth[[j]] <- c(ks1$statistic,ks2$statistic,ks3$statistic,ks4$statistic,ks5$statistic)
				
				PlotTwoDistribution(proportion[2,],proportion[3,],tagref,tagtest,10,names(data$metas[[i]][j]),"Proportion",xlab = "Proportion of Contigs/Reference node per graph",ylab = "Number of connected components",main = "Porportion")
				PlotTwoDistribution(artpointlocal[2,],artpointlocal[3,],tagref,tagtest,10,names(data$metas[[i]][j]),"ArticulationPoint_local",xlab = "Proportion of Contigs/Reference Nodes which are articulation point per component",ylab = "Number of connected components",main = "Articulation Local")
				PlotTwoDistribution(artpointglobal[2,],artpointglobal[3,],tagref,tagtest,10,names(data$metas[[i]][j]),"ArticulationPoint_global",xlab = "Proportion of Contigs/Reference Nodes which are articulation point per component",ylab = "Number of connected components",main = "Articulation Global")
				PlotTwoDistribution(degree1[2,],degree1[3,],tagref,tagtest,10,names(data$metas[[i]][j]),"DegreeOne",xlab = "Density of connected components according to the proportion of Reference/Testdata",ylab = "Number of connected components",main = "Degree One")
				PlotTwoDistribution(geodesics[2,],geodesics[1,],tagref,tagtest,10,names(data$metas[[i]][j]),"Longest_chain",xlab = "Longest distance between two Testdata points or two Reference Nodes",ylab = "Number of connected components",main = "Geodisc")
				PlotOneDistribution(unlist(jaccard[2,]),tagref,tagtest,10,names(data$metas[[i]][j]),"Jaccard",xlab = paste("Distance between",tagref,"/",tagtest,"and reference nodes per graph", sep = " " ),ylab = "Number of Pair Nodes",main = "Jaccard")
				
				MakeSummaryOfCentrality(fic = "ArticulationPoint_global_KS",tag1 = tagref,tag2 = tagtest,method = attributes(data$metas[i])$names,centrality = "ArticulationGlobal",similarity = names(data$metas[[i]][j]),pvalue = ks3$p.value,distance = ks3$statistic)
				MakeSummaryOfCentrality(fic = "ArticulationPoint_local_KS",tag1 = tagref,tag2 = tagtest,method = attributes(data$metas[i])$names,centrality = "ArticulationLocal",similarity = names(data$metas[[i]][j]),pvalue = ks2$p.value,distance = ks2$statistic)
				MakeSummaryOfCentrality(fic = "DegreeOne_KS",tag1 = tagref,tag2 = tagtest,method = attributes(data$metas[i])$names,centrality = "Degree1",similarity = names(data$metas[[i]][j]),pvalue = ks4$p.value,distance = ks4$statistic)
				MakeSummaryOfCentrality(fic = "Longest_chain_KS",tag1 = tagref,tag2 = tagtest,method = attributes(data$metas[i])$names,centrality = "Geodesic",similarity = names(data$metas[[i]][j]),pvalue = ks5$p.value,distance = ks5$statistic)
				MakeSummaryOfCentrality(fic = "Proportion_KS",tag1 = tagref,tag2 = tagtest,method = attributes(data$metas[i])$names,centrality = "Porportion",similarity = names(data$metas[[i]][j]),pvalue = ks1$p.value,distance = ks1$statistic)
				
				SumUpResults( numBadPvalForThisMeth[j], sum(unlist(sapply(CC[goodCC],HasTag,tag = tagtest))), numOfContig, tag = tagtest, sim = names(data$metas[[i]][j]))
				
				
				cl <- makeSOCKcluster(rep("localhost",data$nproc),type = "SOCK")
				clusterEvalQ(cl, library(igraph))
				
				clusterExport(cl, "WriteCC")
				clusterExport(cl, "WriteSeq")
				clusterExport(cl, "tagtest")
				clusterExport(cl, "data")
				tasks <- list(
				s <- function () sapply(CC[which(Allresult[6,]==1)],WriteCC,sim = names(data$metas[[i]][j]),filename = "SafeGeneFamilies_proportion",tag = tagtest),
				s <- function () sapply(CC[which(Allresult[6,]==1)],WriteSeq,sim = names(data$metas[[i]][j]),filename = "SafeContigs_proportion",tag = tagtest),
				s <- function () sapply(CC[which(Allresult[6,]==0)],WriteCC,sim = names(data$metas[[i]][j]),filename = "SuspiciousGeneFamilies_proportion",tag = tagtest),
				s <- function () sapply(CC[which(Allresult[6,]==0)],WriteSeq,sim = names(data$metas[[i]][j]),filename = "SuspiciousContigs_proportion",tag = tagtest),
				
				s <- function () sapply(CC[which(Allresult[5,]==1)],WriteCC,sim = names(data$metas[[i]][j]),filename = "SafeGeneFamilies_articulation_local",tag = tagtest),
				s <- function () sapply(CC[which(Allresult[5,]==1)],WriteSeq,sim = names(data$metas[[i]][j]),filename = "SafeContigs_articulation_local",tag = tagtest),
				s <- function () sapply(CC[which(Allresult[5,]==0)],WriteCC,sim = names(data$metas[[i]][j]),filename = "SuspiciousGeneFamilies_articulation_local",tag = tagtest),
				s <- function () sapply(CC[which(Allresult[5,]==0)],WriteSeq,sim = names(data$metas[[i]][j]),filename = "SuspiciousContigs_articulation_local",tag = tagtest),
				
				s <- function () sapply(CC[which(Allresult[4,]==0)],WriteCC,sim = names(data$metas[[i]][j]),filename = "SuspiciousGeneFamilies_articulation_global",tag = tagtest),
				s <- function () sapply(CC[which(Allresult[4,]==0)],WriteSeq,sim = names(data$metas[[i]][j]),filename = "SuspiciousContigs_articulation_global",tag = tagtest),
				s <- function () sapply(CC[which(Allresult[4,]==1)],WriteSeq,sim = names(data$metas[[i]][j]),filename = "SafeContigs_articulation_global",tag = tagtest),
				s <- function () sapply(CC[which(Allresult[4,]==1)],WriteCC,sim = names(data$metas[[i]][j]),filename = "SafeGeneFamilies_articulation_global",tag = tagtest),
				
				s <- function () sapply(CC[which(Allresult[3,]==1)],WriteCC,sim = names(data$metas[[i]][j]),filename = "SafeGeneFamilies_degree1",tag = tagtest),
				s <- function () sapply(CC[which(Allresult[3,]==1)],WriteSeq,sim = names(data$metas[[i]][j]),filename = "SafeContigs_degree1",tag = tagtest),
				s <- function () sapply(CC[which(Allresult[3,]==0)],WriteCC,sim = names(data$metas[[i]][j]),filename = "SuspiciousGeneFamilies_degree1",tag = tagtest),
				s <- function () sapply(CC[which(Allresult[3,]==0)],WriteSeq,sim = names(data$metas[[i]][j]),filename = "SuspiciousContigs_degree1",tag = tagtest),
				
				s <- function () sapply(CC[which(Allresult[2,]==1)],WriteCC,sim = names(data$metas[[i]][j]),filename = "SafeGeneFamilies_longest_chain",tag = tagtest),
				s <- function () sapply(CC[which(Allresult[2,]==1)],WriteSeq,sim = names(data$metas[[i]][j]),filename = "SafeContigs_longest_chain",tag = tagtest),
				s <- function () sapply(CC[which(Allresult[2,]==0)],WriteSeq,sim = names(data$metas[[i]][j]),filename = "SuspiciousContigs_longest_chain",tag = tagtest),
				s <- function () sapply(CC[which(Allresult[2,]==0)],WriteCC,sim = names(data$metas[[i]][j]),filename = "SuspiciousGeneFamilies_longest_chain",tag = tagtest),
				
				s <- function () sapply(CC[which(Allresult[1,]==1)],WriteCC,sim = names(data$metas[[i]][j]),filename = "SafeGeneFamilies_jaccard",tag = tagtest),
				s <- function () sapply(CC[which(Allresult[1,]==1)],WriteSeq,sim = names(data$metas[[i]][j]),filename = "SafeContigs_jaccard",tag = tagtest),
				s <- function () sapply(CC[which(Allresult[1,]==0)],WriteCC,sim = names(data$metas[[i]][j]),filename = "SuspiciousGeneFamilies_jaccard",tag = tagtest),
				s <- function () sapply(CC[which(Allresult[1,]==0)],WriteSeq,sim = names(data$metas[[i]][j]),filename = "SuspiciousContigs_jaccard",tag = tagtest)
				)
				out <- clusterApply( 
				  cl,
				  tasks,
				  function(f) f()
				)
				
				stopCluster(cl)
	
			}
		
			
			centralitiesSummary <- cbind(centralitiesSummary,numGoodPvalForThisMeth)
			
			GoodSim <- which( (numGoodContigForThisMeth==max(numGoodContigForThisMeth)) & (max(numGoodContigForThisMeth) > 0) )
	
			if(length(GoodSim) > 0)
			{
				if(length(GoodSim)==1)
				{
					s <- MakeGoodContigList(contigs = GoodContigForThisMeth[[GoodSim[1]]],filename = as.character(data$metas[[tagtest]][GoodSim[1]]), tag = tagtest)
					s <- MakeSummaryOfThisMethod(method = attributes(data$metas[i])$names, reference = tagref,centrality = centralities,similarity = attributes(data$metas[[tagtest]])$names[GoodSim[1]],pvalue = AllPvalForThisMeth[[GoodSim[1]]],distance = AllDistForThisMeth[[GoodSim[1]]])
					
					FileConn <- file("NetworksSummary.txt","a")
					writeLines(paste(c( attributes(data$metas[i])$names,tagref, attributes(data$metas[[tagtest]])$names[GoodSim[1]], as.character(data$metas[[tagtest]][GoodSim[1]]),	numGoodContigForThisMeth[GoodSim[1]]),collapse="\t"), FileConn)
					close(FileConn)
				}
				else
				{
					totake <- which(numGoodPvalForThisMeth[GoodSim]==max(numGoodPvalForThisMeth[GoodSim]))[1]
					s <- MakeGoodContigList(contigs = GoodContigForThisMeth[[GoodSim[totake]]], filename = as.character(data$metas[[tagtest]][GoodSim[totake]]), tag = tagtest)
					s <- MakeSummaryOfThisMethod(method = attributes(data$metas[i])$names, reference = tagref,centrality = centralities,similarity = attributes(data$metas[[tagtest]])$names[totake],pvalue = AllPvalForThisMeth[[totake]],distance = AllDistForThisMeth[[totake]])
					
					FileConn <- file("NetworksSummary.txt","a")
					writeLines(paste(c( attributes(data$metas[i])$names,tagref, attributes(data$metas[[tagtest]])$names[totake], as.character(data$metas[[tagtest]][GoodSim[totake]]),	numGoodContigForThisMeth[GoodSim[totake]]),collapse="\t"), FileConn)
					close(FileConn)
				}
			}
		
		}
		
	}
	
	write.table(centralitiesSummary, file = "CentralitySummary.txt", sep = "\t",row.names = FALSE,col.names = FALSE, quote = FALSE, qmethod = "escape",append = TRUE)
}
############################################################################################################################################


Main()

print("Computation Successful",quote=FALSE)
