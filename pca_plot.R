require(dplyr)
require(ggplot)
require(ggrepel)
require(data.table)
require(cowplot)
require(vcfR)
require(adegenet)


source("pca_funs.R")


vcffile <-
    "/nfs/nas22/fs2201/biol_impb_group_evo_euler//data/swiss_tetra/gatk3hc_300_tetra_combined.only4.ldeasy_pruned.vcf"



## data prep
vcf <- read.vcfR(vcffile)
gel <- vcfR2genlight.tetra(vcf)

## namefix
s <- fread('samples.tsv',data.table=F) %>% select(sample,newid) %>% unique
n <- data.frame(sample=indNames(gel))
nn <- left_join(n,s) %>% mutate(outname=ifelse(is.na(newid),sample,newid)) %>% select(outname) %>% unique %>% unlist %>% unname

indNames(gel) <- nn

## name populations
pop(gel) <- substr(indNames(gel),1,3)

## filter
f <- filter(quick,ploidy==4)$population
gel <- gel[gel@pop%in%f]


## checks
gel
indNames(gel)
ploidy(gel)

## run PCA
pca0 <- glPcaFast(gel, nf=5)

## save(pca0,file="pca0.RData") # remote
## load("pca0.RData") # local

## Varprop
(vp1 <- pca0$eig[1]/sum(pca0$eig)) # proportion of variation explained by 1st axis
(vp2 <- pca0$eig[2]/sum(pca0$eig)) # proportion of variation explained by 2nd axis
(vp3 <- pca0$eig[3]/sum(pca0$eig)) # proportion of variation explained by 3rd axis

## prep for plotting
pca0d <- makepca0(pca0)

pca0d <- left_join(pca0d,fread("300_quicklist.csv",data.table=F),by=c('pop'='population'))

pca0d$Lineage <- factor(pca0d$class,
                        levels=unique(pca0d$class),
                        labels=c("Hercynian","Ruderal","S. Carpathian","W. Carpathian"),
                        exclude="")

pca0d$pop <- gsub("MOS","MOU",pca0d$pop) # name correction
coords <- pca0d %>% group_by(pop) %>% summarise(m_ev1=mean(ev1),m_ev2=mean(ev2),m_ev3=mean(ev3))

colors <- c(
Hercynian="#A5839E",
Ruderal="#CCC59D",
`S. Carpathian`="#5F4E27",
`W. Carpathian`="#ABBA61",
GOS="deepskyblue4",
MOU="red2"
)

## plots

(pp2 <- pcaplot0(pca0d,1,2,vp1,vp2,x=ev1,y=ev2,
                 color=Lineage
                 )+
     stat_ellipse(data=filter(pca0d,pop%in%c("MOU","GOS")),mapping=aes(ev1,ev2,color=pop),
                  type="euclid",level=2.5,size=.7)+
     geom_text_repel(data = coords, aes(m_ev1,m_ev2,label=pop))+
     scale_color_manual(breaks=unique(pca0d$Lineage),values=colors)+
     theme(legend.position = c(.7,.7),
           legend.background = element_rect(fill=alpha('white',.6),
                                            size=.3,
                                            linetype='solid',
                                            color='black'),
           legend.margin = margin(r=.8,l=.8,t=.2,b=.7,unit='lines'),
           legend.title = element_blank(),
           legend.text = element_text(size=15))+
     guides(color=guide_legend(override.aes=list(alpha=1,shape=15,size=5))))

saveplot(pp2,'pca0_12_20211102')




(pp1 <- pcaplot0(pca0d,1,3,vp1,vp3,x=ev1,y=ev3,
                 color=Lineage
                 )+
     stat_ellipse(data=filter(pca0d,pop%in%c("MOU","GOS")),mapping=aes(ev1,ev3,color=pop),
                  type="euclid",level=2.5,size=.7)+
     geom_text_repel(data = coords, aes(m_ev1,m_ev3,label=pop),min.segment.length=1.5,force=2)+
     scale_color_manual(breaks=unique(pca0d$Lineage),values=colors)+
     theme(legend.position = c(.7,.7),
           legend.background = element_rect(fill=alpha('white',.6),
                                            size=.3,
                                            linetype='solid',
                                            color='black'),
           legend.margin = margin(r=.8,l=.8,t=.2,b=.7,unit='lines'),
           legend.title = element_blank(),
           legend.text = element_text(size=15))+
     guides(color=guide_legend(override.aes=list(alpha=1,shape=15,size=5))))

saveplot(pp1,'~/pro/300_analyses/plots/pca0_13_20211102')



