require(qqman)
require(ggplot2)
require(dplyr)
require(ggrepel)

getSpeciesName <- function(mgsName, taxo) {
  species = taxo[match(mgsName, rownames(taxo)),"species"]
  return(species)
}

prepare_manhatan_ggplot2 <- function(es_df){
  res_df<- es_df %>% group_by(CHR) %>% summarise(chr_len=max(BP)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>% select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(es_df, ., by=c("CHR"="CHR")) %>%
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>% mutate( BPcum=BP+tot)
  return(res_df)
}
##load taxonomy dataframe
taxo<-read.csv("../taxo.csv",row.names=1)

#Define order and colors by cohort, info in cohort_color file
cohort_color<-read.csv("cohort_color.tsv", header=F)
cohort_color$newchr <-as.integer(rownames(cohort_color))

#Read ouptup from python script effect size
ESdata<-read.csv("Effect_size.Country_ctrl.tsv", sep="\t", row.name=1)
ESdata$species<-getSpeciesName(ESdata$msp,taxo)

efup<-ESdata[,c(1,2,3,5,6,8)]
names(efup)=c("SNP","BP","P","cohort","CHR","species")
#change CHR with final desired order
efup<-cohort_color %>% left_join(efup, ., by=c("CHR"="V2"))
efup$CHR<-efup$newchr

#effect size down column
efdw<-ESdata[,c(1,2,4,5,6,8)]
names(efdw)=c("SNP","BP","P","cohort", "CHR","species")
efdw<-cohort_color %>% left_join(efdw, ., by=c("CHR"="V2"))
efdw$CHR<-efdw$newchr

#plot effect size up all MSP
don <- prepare_manhatan_ggplot2(efup)
axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
#create labels based following final order
axisdf <- cohort_color  %>%  select(-V2) %>% left_join(axisdf, ., by=c("CHR"="newchr"))
#plot effec upp
#pdf(file= "EffectSize_enriched.pdf", width=6, height=4.5 )
jpeg(f="EffectSize_enriched.jpg", w=3*300, h=1.5*300,pointsize=10)
  ggplot(don, aes(x=BPcum, y=P )) +
      # Show all points
      # remove space between plot area and x axis
      # Custom the theme:
      geom_point( aes(color=as.factor(V3)), alpha=0.8, size=1.3) +
      labs(title="Enriched in disease", col="Category", x="Disease cohort", y="Effect Size")+
      geom_hline(yintercept=0.3, linetype="dashed",color = "blue", size=0.5)+
      geom_hline(yintercept=0.5, linetype="dashed",color = "red", size=0.5)+
      geom_text_repel(data = don[which(don$SNP %in% c('msp_1234') & don$P> 0.3) ,],
                      aes(x=BPcum,y=P,label=paste(species,SNP,sep="\n")), segment.color = 'grey', size=5, max.overlaps=10) +
      #scale_color_manual(values = rep(c("grey", "skyblue"), 39 )) +
      # custom X axis:
      scale_x_continuous( label = axisdf$V1, breaks= axisdf$center ) +
      scale_y_continuous(expand = c(0, 0) ) +
      ylim(0,1.05)+
      theme_bw() +
      theme(
        plot.title=element_text(hjust=0, size=20),
        legend.position = 'right',
        legend.text=element_text(size=15),
        legend.title=element_blank(),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x=element_text(angle = 70,hjust=0.99,vjust=0.95, size=13 )

      )
dev.off()

#Plot Effect size down all MSP
don <- prepare_manhatan_ggplot2(efdw)
axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
#create labels based following final order
axisdf <- cohort_color  %>%  select(-V2) %>% left_join(axisdf, ., by=c("CHR"="newchr"))
#plot effec upp
#pdf(file= "EffectSize_depeleted.pdf", width=6, height=4.5 )
jpeg(f="EffectSize_depleted.jpg", w=3*300, h=1.5*300,pointsize=10)
  ggplot(don, aes(x=BPcum, y=P )) +
    # Show all points
    # remove space between plot area and x axis
    # Custom the theme:
    geom_point( aes(color=as.factor(V3)), alpha=0.8, size=1.3) +
    labs(title="Depleated in disease", col="Category", x="Disease cohort", y="Effect Size")+
    geom_hline(yintercept=0.3, linetype="dashed",color = "blue", size=0.5)+
    geom_hline(yintercept=0.5, linetype="dashed",color = "red", size=0.5)+
    geom_text_repel(data = don[which(don$SNP %in% c('msp_0107') & don$P> 0.3),],
                  aes(x=BPcum,y=P,label=paste(species,SNP,sep="\n")), segment.color = 'grey', size=5, max.overlaps=Inf) +
    #scale_color_manual(values = rep(c("grey", "skyblue"), 39 )) +
    # custom X axis:
    scale_x_continuous( label = axisdf$V1, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +
    ylim(0,1.05)+
    theme_bw() +
    theme(
      plot.title=element_text(hjust=0, size=20),
      legend.text=element_text(size=15),
      legend.position = 'right',
      legend.title= element_blank(),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x=element_text(angle = 70,hjust=0.99,vjust=0.95, size=13 )
    )
dev.off()

#MSP entriched and depleted in >5 different cohorts
msp_up=c("msp_0610","msp_0009","msp_0086","msp_0056","msp_0020","msp_0213", "msp_1327", "msp_0331")
msp_dw=c("msp_0107","msp_0259","msp_1291","msp_0188","msp_0436")

#for cycle for each msp_up msp_dw
for (msp in msp_up) {
  #filter efup or efdw by MSP in msp_up/dw list
  mspdf=subset(efup, SNP == msp,)
  don <- prepare_manhatan_ggplot2(mspdf)
  axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  axisdf <- cohort_color  %>%  select(-V2) %>% left_join(axisdf, ., by=c("CHR"="newchr"))
  pdf(file= paste(msp,"enriched","pdf",sep=".") , width=6, height=4.5 )
  print(ggplot(don, aes(x=BPcum, y=P )) +
        # Show all points
        # remove space between plot area and x axis
        # Custom the theme:
        geom_point( aes(color=as.factor(V3)), alpha=0.8, size=1.3) +
        labs(title=msp, col="Category", x="Disease cohort", y="Effect Size")+
        geom_hline(yintercept=0.3, linetype="dashed",color = "blue", size=0.5)+
        geom_hline(yintercept=0.5, linetype="dashed",color = "red", size=0.5)+

        #scale_color_manual(values = rep(c("grey", "skyblue"), 39 )) +
        # custom X axis:
        scale_x_continuous( label = axisdf$V1, breaks= axisdf$center ) +
        scale_y_continuous(expand = c(0, 0) ) +
        ylim(0,1.05)+
        theme_bw() +
        theme(
          plot.title=element_text(hjust=0),
          legend.position = 'right',
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x=element_text(angle = 90,hjust=0.95,vjust=0.2, size=7 )
        )
      )
  dev.off()
}

for (msp in msp_dw) {
  #filter efup or efdw by MSP in msp_up/dw list
  mspdf=subset(efdw, SNP == msp,)
  don <- prepare_manhatan_ggplot2(mspdf)
  axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  axisdf <- cohort_color  %>%  select(-V2) %>% left_join(axisdf, ., by=c("CHR"="newchr"))
  pdf(file= paste(msp,"depelted","pdf",sep=".") , width=6, height=4.5 )
  print(ggplot(don, aes(x=BPcum, y=P )) +
        # Show all points
        # remove space between plot area and x axis
        # Custom the theme:
        geom_point( aes(color=as.factor(V3)), alpha=0.8, size=1.3) +
        labs(title=msp, col="Category", x="Disease cohort", y="Effect Size")+
        geom_hline(yintercept=0.3, linetype="dashed",color = "blue", size=0.5)+
        geom_hline(yintercept=0.5, linetype="dashed",color = "red", size=0.5)+

        #scale_color_manual(values = rep(c("grey", "skyblue"), 39 )) +
        # custom X axis:
        scale_x_continuous( label = axisdf$V1, breaks= axisdf$center ) +
        scale_y_continuous(expand = c(0, 0) ) +
        ylim(0,1.05)+
        theme_bw() +
        theme(
          plot.title=element_text(hjust=0),
          legend.position = 'right',
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x=element_text(angle = 90,hjust=0.95,vjust=0.2, size=7 )
        )
      )
  dev.off()
}


'''
manhattan(efup, genomewideline=0.5, suggestiveline = 0.3, ylim = c(0, 1.1), #col=cohort_cols
      logp=F, #chrlabs=label_cohort,
       xlab="Disease cohort", ylab="Effect Size",
       cex.axis=0.5, las=2, main="Enriched in disease")  # main = targetMain = "a. histaminiformans"
manhattan(efdw, genomewideline=0.5, suggestiveline = 0.3, ylim = c(0, 1.1), #col=cohort_cols
     logp=F, #chrlabs=label_cohort,
      xlab="Disease cohort", ylab="Effect Size",
      cex.axis=0.5, las=2, main="Depleted in disease")
'''
