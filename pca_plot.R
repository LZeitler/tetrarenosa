require(tidyverse)
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

## Varprop
(vp1 <- pca0$eig[1]/sum(pca0$eig)) # proportion of variation explained by 1st axis
(vp2 <- pca0$eig[2]/sum(pca0$eig)) # proportion of variation explained by 2nd axis
(vp3 <- pca0$eig[3]/sum(pca0$eig)) # proportion of variation explained by 3rd axis

## prep for plotting
pca0d <- makepca0(pca0)
fwrite(pca0d,'/nfs/nas22/fs2201/biol_impb_group_evo_euler/lz/temp/pca0d-hc-all-easyldpruned-4x.txt')

pca0d <- fread('/mnt/evo_euler/lz/temp/pca0d-hc-all-easyldpruned-4x.txt',data.table = F)

pca0d$new <- sapply(seq_along(pca0d$sample.id), function(i)
    grepl('z',substr(pca0d$sample.id[i],nchar(pca0d$sample.id[i]),nchar(pca0d$sample.id[i]))))
pca0d$pop <- gsub("MOS","MOU",pca0d$pop)
coords <- pca0d %>% group_by(pop) %>% summarise(m_ev1=mean(ev1),m_ev2=mean(ev2),m_ev3=mean(ev3))

## vp1 <- 0.0409657
## vp2 <- 0.03231271
## vp3 <- 0.03072874
pca0d <- left_join(pca0d,quick,by=c('pop'='population'))
pca0d$Lineage <- factor(pca0d$class,
                        levels=unique(pca0d$class),
                        labels=c("Swiss","Hercynian","Ruderal","S. Carpathian","W. Carpathian"),
                        exclude="")


(pp1 <- pcaplot0(pca0d,1,3,vp1,vp3,x=ev1,y=ev3,
                 color=Lineage
                 )+
     geom_text_repel(data = coords, aes(m_ev1,m_ev3,label=pop)))

saveplot(pp1,'pca_hc_all_vcfR_tetra2_easyldpruned_vp13')


(pp1 <- pcaplot0(pca0d,1,2,vp1,vp2,x=ev1,y=ev2,
                 color=Lineage
                 )+
     geom_text_repel(data = coords, aes(m_ev1,m_ev2,label=pop)))

saveplot(pp1,'pca_hc_all_vcfR_tetra2_easyldpruned_vp12')



