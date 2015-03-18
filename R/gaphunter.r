gaphunter<-function(object,Beta=NULL,threshold=0.05,keepoutliers=FALSE,outcutoff=0.01,verbose=TRUE){
if (is.null(Beta)&is.null(object))
	stop("[gaphunter] At least one of 'object' or 'Beta' must be supplied.")
if ((threshold <= 0)|(threshold >= 1))
	stop("[gaphunter] 'threshold' must be between 0 and 1.") 
if ((outcutoff <= 0)|(outcutoff >= 0.5))
	stop("[gaphunter] 'outcutoff' must be between 0 and 0.5.") 
if (is.null(Beta)){
	if(verbose)
		message("[gaphunter] Calculating beta matrix")
		Beta<-getBeta(object)
	}
if (verbose)
	message("[gaphunter] Using ",prettyNum(nrow(Beta),big.mark=",",scientific=FALSE)," probes and ",prettyNum(ncol(Beta),big.mark=",",scientific=FALSE)," samples.")
if (!getDoParRegistered())
        registerDoSEQ()
    workers <- getDoParWorkers()
    backend <- getDoParName()
    version <- getDoParVersion()
    if (verbose) {
        if (workers == 1) {
            mes <- "[gaphunter] Using a single core (backend: %s, version: %s)."
            message(sprintf(mes, backend, version))
        }
        else {
            mes <- "[gaphunter] Parallelizing using %s workers/cores (backend: %s, version: %s)."
            message(sprintf(mes, workers, backend, version))
        }
    }
if (verbose)
	message("[gaphunter] Searching for gap signals.")
	
chunksize <- ceiling(nrow(Beta)/workers) 

tmp <- foreach(probes = iter(Beta, by = "row", chunksize = chunksize)) %dorng%
	{apply(probes,1,function(k){
		atprobe<-sort(as.numeric(k))
		diffprobe<-as.numeric(lapply(2:ncol(Beta),function(x){atprobe[x]-atprobe[(x-1)]}))
		gaps<-(sum(diffprobe>threshold))
		if (gaps>0){
			grouptemp<-rep(1,ncol(Beta))
			gaptemp<-rep(0,(gaps+2))
			breakpoint<-which(diffprobe>threshold)
			for (j in 1:length(breakpoint)){
				sortpersonbreak<-breakpoint[j]+1
				personbreak<-which(k==atprobe[sortpersonbreak])[1]
				clust1<-which(k<k[personbreak])
				clust2<-which(k==k[personbreak]|k>k[personbreak])
				grouptemp[clust2]<-(j+1)
			}
			gaptemp[1]<-gaps+1
			gaptemp[2:(2+gaps)]<-table(grouptemp)
		return(list(gapanno = gaptemp,groupanno = grouptemp)) 
		}
	}
	)
	}
tmp<-unlist(tmp,recursive=FALSE)
names(tmp)<-rownames(Beta)
myout<-Filter(Negate(is.null),tmp)
groupanno<-do.call("rbind",lapply(myout,function(x)return(x$groupanno)))
pregaps<-lapply(myout,function(x)return(x$gapanno))
mymax<-max(unlist(lapply(pregaps,length)))
gapanno<-do.call("rbind",lapply(pregaps,function(x){length(x)<-mymax; return(x)}))
rownames(gapanno)<-names(myout)
rownames(groupanno)<-names(myout)

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
			if(sum(analyze[-maxgroup])<outcutoff*ncol(Beta)){
			return(blah)}
		}}))
	if (length(markme)>0){
	gapanno<-gapanno[-markme,]
	groupanno<-groupanno[-markme,]}
	if (verbose)
	message("[gaphunter] Removed ",prettyNum(length(markme),big.mark=",",scientific=FALSE)," gaps driven by outliers from results.")
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