prepVCF <- function(vcfdir,bcftools.dir=NULL,bcftools.loadpath=NULL,chromosomes,temp.dir){
    if (!is.null(bcftools.dir) && !is.null(bcftools.loadpath))
        stop("Please supply bcftools load path or install directory.")
    if (is.null(bcftools.dir) && is.null(bcftools.loadpath))
        stop("Please supply bcftools load path or install directory.")

    if (!is.null(bcftools.loadpath)) { # if it has to be preloaded then the command will just be
        bcftools.dir <- "bcftools"     # this
        bcftools.loadpath <-           # and the commands have to be executed together
            paste0(bcftools.loadpath," && ")
    } else {
        bcftools.loadpath <- ""        # else is just empty
    }

    for (c in chromosomes){
        cat("\nPreparing chromosome",c,'\n')
        commands <- paste(bcftools.loadpath,"./ld_easy_prepvcf.sh",bcftools.dir,vcfdir,c,temp.dir)
        ## cat("Executing:\n",commands,'\n')
        system(commands)
    }
    cat("Looks like we're done here.\n")
}

postVCF <- function(vcfdir,bcftools.dir=NULL,bcftools.loadpath=NULL,temp.dir){
    if (!is.null(bcftools.dir) && !is.null(bcftools.loadpath))
        stop("Please supply bcftools load path or install directory.")
    if (is.null(bcftools.dir) && is.null(bcftools.loadpath))
        stop("Please supply bcftools load path or install directory.")

    if (!is.null(bcftools.loadpath)) { # if it has to be preloaded then the command will just be
        bcftools.dir <- "bcftools"     # this
        bcftools.loadpath <-           # and the commands have to be executed together
            paste0(bcftools.loadpath," && ")
    } else {
        bcftools.loadpath <- ""        # else is just empty
    }

    commands <- paste(bcftools.loadpath,"./ld_easy_postvcf.sh",bcftools.dir,vcfdir,temp.dir)
    system(commands)
}

getCorr <- function(locus1,locus2){
    suppressWarnings(
        abs(cor(as.vector(unlist(locus1)),as.vector(unlist(locus2)),"complete.obs"))
    )
}

getCorrWindowed <- function(loci1,loci2){ # not used
    m1 <- matrix(t(loci1),ncol=nrow(loci1))
    m2 <- matrix(t(loci2),ncol=nrow(loci2))
    r <- try(suppressWarnings(
        mean(abs(cor(m1,m2,"complete.obs")**2))
    ),T)
    if(!is.numeric(r)) 1 else r
}

jumper <- function(focal,limit,rmax,direction,m=m){
    j <- focal                  # start/focal site
    jj <- j+direction           # next site
    o <- vector()
    while ( abs(jj-limit) > 1 ){
        rr <- getCorr(m[j,],m[jj,])
        if(rr<rmax){            # if corr is low, 
            o <- c(o,jj)        # add locus to output,
            j <- jj             # new focal is old next
            jj <- jj+direction  # and next site also by 1
        } else {
            jj <- jj+direction  # else just step next site and leave focal as is
        }
    }
    return(o)
}

runPruning <- function(rlimit,windowsize,get.loci=T,get.geno=F,l=l,m=m){
    wins <- c(seq(l[1],l[length(l)],windowsize),l[length(l)])
    s <- vector()
    pb <- txtProgressBar(min=0,max=length(wins),style=3)
    for (i in 1:(length(wins)-1)){
        start <- which.min(abs(l-wins[i]))
        stop <- which.min(abs(l-wins[i+1]))
        if (abs(start-stop)<=1) next 
        focal <- sample(start:(stop),1)
        s <- c(s,
               jumper(focal,start,rlimit,-1,m=m), # this is reverse jumping
               jumper(focal,stop,rlimit,1,m=m))   # this is forward jumping
        setTxtProgressBar(pb,i)
    }
    close(pb)

    s <- sort(s)
    d <- m[s,]
    l <- l[s]
    cat("Selected",nrow(d),"out of",nrow(m),"loci,",round(nrow(d)/nrow(m)*100),"%. \n")
    if (get.geno) d
    else if (get.loci) l
    else NULL
}

pruneChromosomes <- function(temp.dir,windowsize,rlimit){
    files <- list.files(paste0(temp.dir,'/ld_easy_temp'))
    mfiles <- files[grep('matrix',files)]
    run <- F

    for (f in mfiles){
        
        ## load files
        m <- fread(paste0(temp.dir,'/ld_easy_temp/',f))
        c <- gsub("matrix_|.txt","",f)
        l <- as.vector(fread(
            paste0(temp.dir,'/ld_easy_temp/positions_',c,'.txt'),
            data.table=F)[,1]
            )

        ## run pruning for chromosome
        cat("\nStarting with pruning chromosome",c,":\n")
        d <- runPruning(rlimit,windowsize,m=m,l=l)

        ## prepare regions file for bcftools
        cat("\nWriting sites to regions file.\n")
        fwrite(data.frame(c,d),paste0(temp.dir,'/ld_easy_temp/regions.txt'),sep='\t',quote=F,col.names=F,append=run)
        run <- T
    }
}

