gaphunter<-function(object,Beta=NULL,threshold=0.05,keepoutliers=FALSE,outcutoff=0.01,verbose=TRUE){
if (is.null(Beta)&is.null(object))
	stop("[gaphunter] At least one of 'object' or 'Beta' must be supplied.")
if ((threshold <= 0)|(threshold >= 1))
	stop("[gaphunter] 'threshold' must be between 0 and 1.") 
if ((outcutoff <= 0)|(outcutoff >= 0.5))
	stop("[gaphunter] 'outcutoff' must be between 0 and 0.5.") 
if (is.null(Beta)){
	if(verbose)
		message("[gaphunter] Calculating beta matrix.")
		Beta<-getBeta(object)
	}
if (verbose){
	message("[gaphunter] Using ",prettyNum(nrow(Beta),big.mark=",",scientific=FALSE)," probes and ",prettyNum(ncol(Beta),big.mark=",",scientific=FALSE)," samples.")
	message("[gaphunter] Searching for gap signals.")
}

sorting<-t(apply(Beta,1,sort.int,method="quick",index.return=TRUE))
sortedindices<-lapply(sorting,function(n)return(n$ix))
sortedbeta<-do.call("rbind",lapply(sorting,function(n)return(n$x)))
rownames(sortedbeta)<-rownames(Beta)
diffs<-rowDiffs(sortedbeta)
gapind<-rowSums(diffs>threshold)
sortedindices<-lapply(which(gapind>0),function(h){return(sortedindices[[h]])})
sortedbeta<-sortedbeta[which(gapind>0),]
diffs<-diffs[which(gapind>0),]
gapind<-gapind[which(gapind>0)]

breakpoints<-apply(diffs,1,function(x){return(which(x>threshold))})

returngroups<-function(x,y){
	template<-rep(1,ncol(sortedbeta))
	count<-2
	for (j in x){
	template[y[(j+1):length(y)]]<-count
	count<-count+1
	}
	return(template)
}	

groupanno<-t(mapply(returngroups,breakpoints,sortedindices))
gapanno<-apply(groupanno,1,function(k){
	tempgap<-rep(0,(length(unique(k))+1))
	tempgap[1]<-length(unique(k))
	tempgap[2:length(tempgap)]<-table(k)
	return(tempgap)})
gapanno<-do.call("rbind",lapply(gapanno,function(x){length(x)<-(max(gapind)+2); return(x)})) 
rownames(gapanno)<-rownames(sortedbeta) ##
rownames(groupanno)<-rownames(sortedbeta)

if(verbose)
message("[gaphunter] Found ",prettyNum(nrow(gapanno),big.mark=",",scientific=FALSE)," gap signals.")

if (keepoutliers==FALSE){
	if (verbose)
	message("[gaphunter] Filtering out gap signals driven by outliers.")
		markme<-unlist(lapply(1:nrow(gapanno),function(blah){
		numgroup<-gapanno[blah,1]
		analyze<-gapanno[blah,(2:(1+numgroup))]
		maxgroup<-which(analyze==max(analyze))
		if (length(maxgroup)==1){	
			if(sum(analyze[-maxgroup])<outcutoff*ncol(Beta))			
			{return(blah)}
		}}))
	if (length(markme)>0){
	gapanno<-gapanno[-markme,]
	groupanno<-groupanno[-markme,]}
	if (verbose)
	message("[gaphunter] Removed ",prettyNum(length(markme),big.mark=",",scientific=FALSE)," gap signals driven by outliers from results.")
}	
gapanno<-data.frame(gapanno)

labs<-c()
for (i in 1:(ncol(gapanno)-1)){
	temp<-paste("Group",i,sep="")
	labs<-c(labs,temp)
}
colnames(gapanno)<-c("Groups",labs)

algorithm<-list("threshold"=threshold,"outcutoff"=outcutoff,"keepoutliers"=keepoutliers)

return(list("proberesults"=gapanno,"sampleresults"=groupanno,"algorithm"=algorithm))
}