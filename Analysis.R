###############################################Required packages####################################################################
library(biomformat)
library(tidyverse)
library(data.table)
library(fs)
library(phyloseq)
library(writexl)
library(phylostratr)
library(ape)  # For phylogenetic trees
library(taxize)  # To fetch taxonomic names from NCBI
######################Kaiju data transformation##################################################################################
Kaijubioms<-dir("Kaiju Taxonomy Abundance",pattern="*biom")
Kaijutaxa<-dir("Kaiju Taxonomy Abundance",pattern="kaiju_taxonomy*")
for (biom in Kaijubioms){
  for(taxa in Kaijutaxa){
    if(substring(biom,21,nchar(biom)-5)==substring(taxa,16,nchar(taxa)-4))
      name=substring(biom,21,nchar(biom)-5)
    OTUtable<-as.data.frame(as.matrix(biom_data(read_biom(paste("Kaiju Taxonomy Abundance",biom,sep="/")))))%>%
      rownames_to_column("Feature.ID")
    taxatable<-read.table(paste("Kaiju Taxonomy Abundance",taxa,sep="/"),sep='\t',header = TRUE)
    Abundance<-inner_join(OTUtable,taxatable)%>%
      select(-"Feature.ID")%>%
      relocate(Taxon)%>%
      mutate(Taxon=str_replace_all(Taxon,".__",""))%>%
      separate(Taxon,c("Domain","Phylum","Class","Order","Family","Genus","Species"),sep=";")%>%
      mutate(Phylum=if_else(is.na(Phylum),Domain,Phylum),
             Class=if_else(is.na(Class),Phylum,Class),
             Order=if_else(is.na(Order),Class,Order),
             Family=if_else(is.na(Family),Order,Family),
             Genus=if_else(is.na(Genus),Family,Genus),
             Species=if_else(is.na(Species),Genus,Species))
    write_xlsx(Abundance,paste("Kaiju Taxonomy Abundance/kaiju_abun-",name,".xlsx",sep=""))
  }
}

######################Kaiju phyloseq object##################################################################################

Kaijubioms<-dir("Kaiju Taxonomy Abundance",pattern="*biom")
Kaijutaxa<-dir("Kaiju Taxonomy Abundance",pattern="kaiju_taxonomy*")
for (biom in Kaijubioms){
  for(taxa in Kaijutaxa){
    if(substring(biom,21,nchar(biom)-5)==substring(taxa,16,nchar(taxa)-4))
      name=substring(biom,21,nchar(biom)-5)
    OTUtable<-as.matrix(biom_data(read_biom(paste("Kaiju Taxonomy Abundance",biom,sep="/"))))
    taxatable<-read.table(paste("Kaiju Taxonomy Abundance",taxa,sep="/"),sep='\t',header = TRUE)%>%
      mutate(Taxon=str_replace_all(Taxon,".__",""))%>%
      separate(Taxon,c("Domain","Phylum","Class","Order","Family","Genus","Species"),sep=";")%>%
      mutate(Phylum=if_else(is.na(Phylum),Domain,Phylum),
             Class=if_else(is.na(Class),Phylum,Class),
             Order=if_else(is.na(Order),Class,Order),
             Family=if_else(is.na(Family),Order,Family),
             Genus=if_else(is.na(Genus),Family,Genus),
             Species=if_else(is.na(Species),Genus,Species))%>%
      column_to_rownames("Feature.ID")%>%
      as.matrix()
    physeq <- phyloseq(otu_table(OTUtable, taxa_are_rows = TRUE), tax_table(taxatable))
    saveRDS(physeq,paste("Kaiju Taxonomy Abundance/phyloseq_kaiju-",name,".rds",sep=""))
  }
}



##################################Creating Bracken Taxonomic abundance Table##################################################### 

brackentables<-dir(path="BrackenOutputs",pattern="table-bracken*")
brackentables<-dir_info("BrackenOutputs",regexp = "table-bracken*")$path[dir_info("BrackenOutputs",regexp = "table-bracken*")$modification_time>as.POSIXct("2025-02-06")]

brackentaxonomy<-dir(path="Bracken Outputs",pattern="taxonomy-bracken*")
brackentaxonomy<-dir_info("BrackenOutputs",regexp = "taxonomy-bracken*")$path[dir_info("BrackenOutputs",regexp = "taxonomy-bracken*")$modification_time>as.POSIXct("2025-02-06")]

dir.create("Kraken-Bracken Taxonomy Abundance")
for(each in brackentables){
  name<-substring(each,30)
  for(each2 in brackentaxonomy){
    name2<-substring(each2,33)
    if(name==name2)
    {table<-as.data.frame(as.matrix(biom_data(read_biom(paste(each,"feature-table.biom",sep='/')))))%>%
      rownames_to_column("Feature.ID")
    taxonomy<-read.table(paste(each2,"taxonomy.tsv",sep='/'),sep="\t",header = TRUE)%>%
      mutate(Feature.ID=as.character(Feature.ID))
    taxatable<-inner_join(table,taxonomy)%>%
      select(-"Feature.ID")%>%
      relocate(Taxon)%>%
      mutate(Taxon=str_replace_all(Taxon,".__",""))%>%
      separate(Taxon,c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species"),sep=";")
    #write_tsv(taxatable,paste("Kraken-Bracken Taxonomy Abundance/Kraken-Bracken_abun-",name,".tsv",sep=""))
    }
  }
}

Allfiles<-NULL
path<-"Kraken Outputs/New/kraken_reports-demux-paired"
filenames<-list.files(path)
folder<-"Raw"
for(each in filenames){
  file<-substring(each,1,(nchar(each)-11))
  file1<-read.table(paste(path,each,sep="/"),sep='\t',header = FALSE)%>%
    as.data.frame()%>%
    setNames(c("percent","clade_reads", "direct_reads", "rank_code", "tax_id", "name"))%>%
    mutate(name=trimws(name),
           sample=file,
           datatype=folder)%>%
    filter(name %in% c("unclassified","cellular organisms","Bacteria","Viruses","Eukaryota"))%>%
    select(-c(3,4,5))%>%
    relocate(datatype,sample,name)
  Allfiles<-Allfiles%>%
    rbind(file1)
}
write_tsv(Allfiles,paste(folder,"report.tsv"))

Allfiles<-NULL
path<-"C:/Users/Ankan_Choudhury/Box/M_COLON/Stool Sample Sequencing/MRDNA/Mahmud/A_kraken_results"
filenames<-list.files(path,pattern="*kraken_report.txt")
folder<-"Trimmed-filtered-assembled"
for(each in filenames){
  file<-substring(each,1,(nchar(each)-26))
  file1<-read.table(paste(path,each,sep="/"),sep='\t',header = FALSE)%>%
    as.data.frame()%>%
    setNames(c("percent","clade_reads", "direct_reads", "rank_code", "tax_id", "name"))%>%
    mutate(name=trimws(name),
           sample=file,
           datatype=folder)%>%
    filter(name %in% c("unclassified","cellular organisms","Bacteria","Viruses","Eukaryota"))%>%
    select(-c(3,4,5))%>%
    relocate(datatype,sample,name)
  Allfiles<-Allfiles%>%
    rbind(file1)
}
write_tsv(Allfiles,paste(folder,"report_Mahmud.tsv"))


##################################Creating Bracken Phyloseq##################################################### 

brackentables<-dir(path="BrackenOutputs/Do",pattern="table-bracken*")
brackentables<-dir_info("BrackenOutputs/Do",regexp = "table-bracken*")$path[dir_info("BrackenOutputs",regexp = "table-bracken*")$modification_time>as.POSIXct("2025-02-06")]

brackentaxonomy<-dir(path="Bracken Outputs",pattern="taxonomy-bracken*")
brackentaxonomy<-dir_info("BrackenOutputs",regexp = "taxonomy-bracken*")$path[dir_info("BrackenOutputs",regexp = "taxonomy-bracken*")$modification_time>as.POSIXct("2025-02-06")]

brackentree<-dir(path="BrackenOutputs",pattern="bracken-reports*")
brackentree<-dir_info("BrackenOutputs",regexp = "bracken-reports*")$path[dir_info("BrackenOutputs",regexp = "bracken-reports*")$modification_time>as.POSIXct("2025-02-06")]


dir.create("BrackenOutputs/BrackenPhyloseq")
counter=0
for(each in brackentables){
  name<-substring(each,30)
  table<-as.matrix(biom_data(read_biom(paste(each,"feature-table.biom",sep='/'))))
  taxonomy<-read.table(paste("BrackenOutputs/taxonomy-bracken-",name,"/taxonomy.tsv",sep=''),sep="\t",header = TRUE)%>%
    mutate(Taxon=str_replace_all(Taxon,".__",""))%>%
    separate(Taxon,c("Domain","Phylum","Class","Order","Family","Genus","Species"),sep=";")%>%
    mutate(Phylum=if_else(is.na(Phylum),Domain,Phylum),
           Class=if_else(is.na(Class),Phylum,Class),
           Order=if_else(is.na(Order),Class,Order),
           Family=if_else(is.na(Family),Order,Family),
           Genus=if_else(is.na(Genus),Family,Genus),
           Species=if_else(is.na(Species),Genus,Species))%>%
    column_to_rownames("Feature.ID")%>%
    as.matrix()
  tree<-read_tree(paste("BrackenOutputs/bracken-reports-",name,"/output.tree",sep=''))
  counter=counter+1
  print(counter)
  physeqBracken<-phyloseq(otu_table(table,taxa_are_rows = TRUE),tax_table(taxonomy))
  physeqBracken<-merge_phyloseq(physeqBracken,tree)
  saveRDS(physeqBracken,paste("BrackenOutputs/BrackenPhyloseq/phyloseq_bracken-",name,".rds",sep=""))
  print(name)
}

#######################################Creating Phylum level abundance data from bracken reports (Mahmud)##########################
#Folder names and data-types for Mahmud, will switch the path and folder variable accordingly: bc_bracken_results=Raw,bracken_results=Trimmed-filtered,A_bracken_results=Trimmed-filtered-contigs

PhylumData_Mahmud<-NULL
path<-"C:/Users/Ankan_Choudhury/Box/M_COLON/Stool Sample Sequencing/MRDNA/Mahmud/bc_bracken_results" 
filenames<-list.files(path,pattern="*bracken_report.txt")
folder<-"Raw"
for(each in filenames){
  file<-substring(each,1,(nchar(each)-26))
  file1<-read.table(paste(path,each,sep="/"),sep='\t',header = FALSE)%>%
    as.data.frame()%>%
    setNames(c("percent","clade_reads", "direct_reads", "rank_code", "tax_id", "name"))%>%
    mutate(name=trimws(name),
           sample=file,
           datatype=folder)%>%
    filter(rank_code=="P")%>%
    select(-c(3,4,5))%>%
    relocate(datatype,sample,name)
  PhylumData_Mahmud<-PhylumData_Mahmud%>%
    rbind(file1)
}

PhylumData_Mahmud_summary<-PhylumData_Mahmud%>%    #summarizing across datatypes
  group_by(datatype,name)%>%
  summarize(mean_percent=round(mean(percent),2),
            SE=sd(percent)/sqrt(n()))%>%
  filter(mean_percent>0.01)

#######################################Creating Phylum level abundance data from bracken reports (Ankan)##########################
#Folder names and data-types for Ankan, will switch the path and folder variable accordingly: demux-paired=Raw,trimmed-filtered=Trimmed-filtered,trimmed-filtered-contigs=Trimmed-filtered-contigs


PhylumData_Ankan<-NULL
path<-"C:/Users/Ankan_Choudhury/Box/M_COLON/Stool Sample Sequencing/MRDNA/Ankan Output/BrackenOutputs/bracken-reports-demux-paired"
filenames<-list.files(path,pattern="*report.txt")
folder<-"Raw"
for(each in filenames){
  print(paste("Processing:",each))
  file<-substring(each,1,(nchar(each)-11))
  file1<-tryCatch({read.table(paste(path,each,sep="/"),sep='\t',header = FALSE)%>%
      as.data.frame()%>%
      setNames(c("percent","clade_reads", "direct_reads", "rank_code", "tax_id", "name"))%>%
      mutate(name=trimws(name),
             sample=file,
             datatype=folder)%>%
      filter(rank_code=="P")%>%
      select(-c(3,4,5))%>%
      relocate(datatype,sample,name)
  }, warning = function(w) {
    message("Warning in file: ", each, "\n", w)  # Print file name if warning occurs
    return(NULL)  # Return NULL so the loop continues
  }, error = function(e) {
    message("Error in file: ", each, "\n", e)  # Print file name if error occurs
    return(NULL)
  })
  PhylumData_Ankan<-PhylumData_Ankan%>%
    rbind(file1)
}

PhylumData_Ankan_summary<-PhylumData_Ankan%>%        #summarizing across datatypes
  group_by(datatype,name)%>%
  summarize(mean_percent=round(mean(percent),2),
            SE=sd(percent)/sqrt(n()))%>%
  filter(mean_percent>0.01)  

#############################################Color palette################################################################
library(Polychrome)
set.seed(123456)
P12 <- createPalette(100, c("#FFFF00","#FF0000", "#00FF00", "#0000FF"), range = c(0, 2000))    #creating color palette
swatch(P12)
P12 <- sortByHue(P12)
P12 <- as.vector(t(matrix(P12, ncol=5)))
names(P12) <- NULL

#####################################Bracken Stacked taxa bar plots##############################################################  
PhylumData_Ankan_summary%>%
  ggplot(aes(x=datatype,y=mean_percent,fill=name))+
  geom_bar(position="fill",stat="identity",color="black",width=1)+
  scale_y_continuous(labels = scales::percent_format())+
  labs(y="Relative Abundance",x="Data Type",fill="Phylum",title="Phylum Relative abundances: Ankan's Bracken Data")+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(values=P12)+theme_classic()+
  scale_x_discrete(expand = c(-2, -2),labels = c("Raw" = "Raw", "Trimmed-filtered" = "Trimmed\nfiltered", "Trimmed-filtered-contigs" = "Trimmed\nfiltered\ncontigs"))+
  coord_equal(9,expand = FALSE,xlim=c(0,4))+
  theme(
    panel.background = element_rect(fill="white"),
    text = element_text(size=22,color="black",face="bold"),
    panel.border = element_rect(size=1,fill="NA"),
    strip.text = element_text(size=22,color="black",face="bold"),
    legend.text = element_text(size=15,color="black",face="italic"),
    legend.margin = unit(0,"cm"),
    legend.title=element_text(size=20,face="bold",color="black"),
    legend.key.size= unit(1, 'cm'),
    plot.title = element_text(hjust=0, size=20),
    plot.margin = margin(0,0,0,0))

PhylumData_Ankan%>%
  group_by(datatype,sample)%>%
  mutate(name=if_else(percent<0.5,"<0.5% abundances in the sample",name))%>%
  group_by(datatype,sample,name)%>%
  summarise(percent=sum(percent))%>%
  ggplot(aes(x=sample,y=percent,fill=reorder(name,-percent)))+
  geom_bar(position="fill",stat="identity",color="black",width=1)+
  scale_y_continuous(labels = scales::percent_format())+
  labs(y="Relative Abundance",x=NULL,fill="Phylum",title="Phylum Relative abundances: Ankan's Bracken Data")+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(values=P12)+theme_classic()+
  ggh4x::facet_nested(~datatype)+
  #scale_x_discrete(expand = c(-2, -2),labels = c("Raw" = "Raw", "Trimmed-filtered" = "Trimmed\nfiltered", "Trimmed-filtered-contigs" = "Trimmed\nfiltered\ncontigs"))+
  #coord_equal(9,expand = FALSE,xlim=c(0,4))+
  theme(
    panel.background = element_rect(fill="white"),
    text = element_text(size=22,color="black",face="bold"),
    axis.text.x = element_blank(),
    panel.border = element_rect(size=1,fill="NA"),
    strip.text = element_text(size=22,color="black",face="bold"),
    legend.text = element_text(size=15,color="black",face="italic"),
    legend.margin = unit(0,"cm"),
    legend.title=element_text(size=20,face="bold",color="black"),
    legend.key.size= unit(1, 'cm'),
    plot.title = element_text(hjust=0, size=20),
    plot.margin = margin(0,0,0,0),
    panel.spacing = unit(0,"mm"))

PhylumData_Mahmud_summary%>%
  ggplot(aes(x=datatype,y=mean_percent,fill=name))+
  geom_bar(position="fill",stat="identity",color="black",width=0.5)+
  scale_y_continuous(labels = scales::percent_format())+
  labs(y="Relative Abundance",x="Data Type",fill="Phylum",title="Relative abundances: Mahmud's Data")+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(values=P12)+theme_classic()+
  theme(
    panel.background = element_rect(fill="white"),
    text = element_text(size=35,color="black",face="bold"),
    panel.border = element_rect(size=1,fill="NA"),
    strip.text = element_text(size=22,color="black",face="bold"),
    legend.text = element_text(size=20,color="black",face="italic"),
    legend.margin = unit(0,"cm"),
    legend.title=element_text(size=25,face="bold",color="black"),
    legend.key.size = unit(1, 'cm'),
    panel.spacing=unit(0,"cm"),
    plot.title = element_text(hjust=0.5, size=40))  

#################################Boxplots faceted by phylum name, with datatype as x-axis########################################

PhylumData_Ankan%>%
  filter(name %in% PhylumData_Ankan_summary$name)%>%
  ggplot(aes(x=datatype,y=percent,fill=datatype))+
  geom_boxplot()+
  scale_x_discrete(labels = c("Raw" = "Raw", "Trimmed-filtered" = "Trimmed\nfiltered", "Trimmed-filtered-contigs" = "Trimmed\nfiltered\ncontigs"))+
  facet_wrap(~name,scales="free_y")+
  guides(fill="none")+
  labs(y="Relative Abundance %",x="Data Type",title="Relative abundances: Ankan's Data")+
  theme_classic()+
  theme(
    panel.background = element_rect(fill="white"),
    text = element_text(size=25,color="black",face="bold",hjust=0.5),
    axis.text = element_text(size=20,color="black",face="bold"),
    panel.border = element_rect(size=1,fill="NA"),
    strip.text = element_text(size=20,color="black",face="bold"),
    panel.spacing=unit(0,"cm"),
    plot.title = element_text(hjust=0.5, size=25)) 

PhylumData_Mahmud%>%
  filter(name %in% PhylumData_Mahmud_summary$name)%>%
  ggplot(aes(x=datatype,y=percent,fill=datatype))+
  geom_boxplot()+
  scale_x_discrete(labels = c("Raw" = "Raw", "Trimmed-filtered" = "Trimmed\nfiltered", "Trimmed-filtered-contigs" = "Trimmed\nfiltered\ncontigs"))+
  facet_wrap(~name,scales="free_y")+
  guides(fill="none")+
  labs(y="Relative Abundance %",x="Data Type",title="Relative abundances: Mahmud's Data")+
  theme_classic()+
  theme(
    panel.background = element_rect(fill="white"),
    text = element_text(size=25,color="black",face="bold",hjust=0.5),
    axis.text = element_text(size=20,color="black",face="bold"),
    panel.border = element_rect(size=1,fill="NA"),
    strip.text = element_text(size=20,color="black",face="bold"),
    panel.spacing=unit(0,"cm"),
    plot.title = element_text(hjust=0.5, size=25)) 


#################################Boxplots faceted by phylum, with data source as x-axis########################################
#switch the 'datatype ==' statement to get figures for Raw, Trimmed-filtered and Trimmed-filtered-contigs data

PhylumData_Ankan %>%
  mutate(source = "Ankan") %>%
  rbind(PhylumData_Mahmud %>% mutate(source = "Mahmud")) %>%
  filter(name %in% PhylumData_Mahmud_summary$name,
         datatype == "Trimmed-filtered-contigs") %>%
  ggplot(aes(x = source, y = percent, fill = datatype)) +
  geom_boxplot(aes(x = source, y = percent, fill = source)) +
  facet_wrap(~name, scales = "free_y") +
  guides(fill = "none") +
  labs(y = "Relative Abundance %", x = "Source", title = paste("Relative abundances: Trimmed, filtered & assembled")) +
  ggpubr::stat_compare_means(method = "aov",size = 7, label.y.npc = 0.9, hjust = 0) +
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "white"),
    text = element_text(size = 25, color = "black", face = "bold", hjust = 0.5),
    axis.text = element_text(size = 20, color = "black", face = "bold"),
    panel.border = element_rect(size = 1, fill = "NA"),
    strip.text = element_text(size = 20, color = "black", face = "bold"),
    panel.spacing = unit(0, "cm"),
    plot.title = element_text(hjust = 0.5, size = 25)
  )


#######################################Creating Phylum level abundance data from kraken reports (Ankan)##########################
#Folder names and data-types for Ankan, will switch the path and folder variable accordingly: kraken_reports-demux-paired=Raw,kraken_reports-trimmed-filtered=Trimmed-filtered,kraken_reports-trimmed-filtered-contigs=Trimmed-filtered-contigs

kraken<-read.table("Kraken Outputs/New/kraken_reports_trimmed-filtered/1-MDcol002-5-2.report.txt",sep='\t',header = FALSE)
PhylumDataKraken_Ankan<-NULL
path<-"Kraken Outputs/New/kraken_reports-trimmed-filtered-contigs"
filenames<-list.files(path)
folder<-"Trimmed-filtered-contigs"
for(each in filenames){
  file<-substring(each,1,(nchar(each)-11))
  file1<-read.table(paste(path,each,sep="/"),sep='\t',header = FALSE)%>%
    as.data.frame()%>%
    setNames(c("percent","clade_reads", "direct_reads", "rank_code", "tax_id", "name"))%>%
    mutate(name=trimws(name),
           sample=file,
           datatype=folder)%>%
    filter(rank_code %in% c("U","P"))%>%
    select(-c(3,4,5))%>%
    relocate(datatype,sample,name)
  PhylumDataKraken_Ankan<-PhylumDataKraken_Ankan%>%
    rbind(file1)
}
PhylumDataKraken_Ankan_summary<-PhylumDataKraken_Ankan%>%
  group_by(datatype,name)%>%
  summarize(mean_percent=round(mean(percent),2),
            SE=sd(percent)/sqrt(n()))%>%
  filter(mean_percent>0.01)  

#######################################Creating Phylum level abundance data from kraken reports (Mahmud)##########################
#Folder names and data-types for Mahmud, will switch the path and folder variable accordingly: bc_kraken_results=Raw,kraken_results=Trimmed-filtered,A_kraken_results=Trimmed-filtered-contigs


PhylumDataKraken_Mahmud<-NULL
path<-"C:/Users/Ankan_Choudhury/Box/M_COLON/Stool Sample Sequencing/MRDNA/Mahmud/A_kraken_results"
filenames<-list.files(path,pattern="*kraken_report.txt")
folder<-"Trimmed-filtered-contigs"
file2<-NULL
for(each in filenames){
  print(paste("Processing:",each))
  file<-substring(each,1,(nchar(each)-26))
  file1<-tryCatch({read.table(paste(path,each,sep="/"),sep='\t',header = FALSE)%>%
      as.data.frame()%>%
      setNames(c("percent","clade_reads", "direct_reads", "rank_code", "tax_id", "name"))%>%
      mutate(name=trimws(name),
             sample=file,
             datatype=folder)%>%
      filter(rank_code %in% c("U","P"))%>%
      select(-c(3,4,5))%>%
      relocate(datatype,sample,name)
  }, warning = function(w) {
    message("Warning in file: ", each, "\n", w)  # Print file name if warning occurs
    return(NULL)  # Return NULL so the loop continues
  }, error = function(e) {
    message("Error in file: ", each, "\n", e)  # Print file name if error occurs
    return(NULL)
  })
  file2<-file2%>%
    rbind(file1)
}
PhylumDataKraken_Mahmud<-PhylumDataKraken_Mahmud%>%
  rbind(file2)
PhylumDataKraken_Mahmud_summary<-PhylumDataKraken_Mahmud%>%
  group_by(datatype,name)%>%
  summarize(mean_percent=round(mean(percent),2),
            SE=sd(percent)/sqrt(n()))%>%
  filter(mean_percent>0.01)  

#####################################Kracken Stacked taxa bar plots##############################################################  

PhylumDataKraken_Mahmud_summary%>%
  ggplot(aes(x=datatype,y=mean_percent,fill=name))+
  geom_bar(position="fill",stat="identity",color="black",width=0.5)+
  scale_y_continuous(labels = scales::percent_format())+
  labs(y="Relative Abundance",x="Data Type",fill="Phylum",title="Relative abundances: Mahmud's Kraken Data")+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(values=P12)+theme_classic()+
  theme(
    panel.background = element_rect(fill="white"),
    text = element_text(size=35,color="black",face="bold"),
    panel.border = element_rect(size=1,fill="NA"),
    strip.text = element_text(size=22,color="black",face="bold"),
    legend.text = element_text(size=20,color="black",face="italic"),
    legend.margin = unit(0,"cm"),
    legend.title=element_text(size=25,face="bold",color="black"),
    legend.key.size = unit(1, 'cm'),
    panel.spacing=unit(0,"cm"),
    plot.title = element_text(hjust=0.5, size=40))

PhylumDataKraken_Ankan_summary%>%
  ggplot(aes(x=datatype,y=mean_percent,fill=name))+
  geom_bar(position="fill",stat="identity",color="black",width=0.5)+
  scale_y_continuous(labels = scales::percent_format())+
  labs(y="Relative Abundance",x="Data Type",fill="Phylum",title="Relative abundances: Ankan's Kraken Data")+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(values=P12)+theme_classic()+
  theme(
    panel.background = element_rect(fill="white"),
    text = element_text(size=35,color="black",face="bold"),
    panel.border = element_rect(size=1,fill="NA"),
    strip.text = element_text(size=22,color="black",face="bold"),
    legend.text = element_text(size=20,color="black",face="italic"),
    legend.margin = unit(0,"cm"),
    legend.title=element_text(size=25,face="bold",color="black"),
    legend.key.size = unit(1, 'cm'),
    panel.spacing=unit(0,"cm"),
    plot.title = element_text(hjust=0.5, size=40))

#################################Boxplots faceted by phylum name, with datatype as x-axis########################################


PhylumDataKraken_Ankan%>%
  filter(name %in% PhylumDataKraken_Ankan_summary$name)%>%
  ggplot(aes(x=datatype,y=percent,fill=datatype))+
  geom_boxplot()+
  scale_x_discrete(labels = c("Raw" = "Raw", "Trimmed-filtered" = "Trimmed\nfiltered", "Trimmed-filtered-contigs" = "Trimmed\nfiltered\ncontigs"))+
  facet_wrap(~name,scales="free_y")+
  guides(fill="none")+
  labs(y="Relative Abundance %",x="Data Type",title="Relative abundances: Ankan's Kraken Data")+
  theme_classic()+
  theme(
    panel.background = element_rect(fill="white"),
    text = element_text(size=25,color="black",face="bold",hjust=0.5),
    axis.text = element_text(size=20,color="black",face="bold"),
    panel.border = element_rect(size=1,fill="NA"),
    strip.text = element_text(size=20,color="black",face="bold"),
    panel.spacing=unit(0,"cm"),
    plot.title = element_text(hjust=0.5, size=25)) 

PhylumDataKraken_Mahmud%>%
  filter(name %in% PhylumDataKraken_Mahmud_summary$name)%>%
  ggplot(aes(x=datatype,y=percent,fill=datatype))+
  geom_boxplot()+
  scale_x_discrete(labels = c("Raw" = "Raw", "Trimmed-filtered" = "Trimmed\nfiltered", "Trimmed-filtered-contigs" = "Trimmed\nfiltered\ncontigs"))+
  facet_wrap(~name,scales="free_y")+
  guides(fill="none")+
  labs(y="Relative Abundance %",x="Data Type",title="Relative abundances: Mahmud's Kraken Data")+
  theme_classic()+
  theme(
    panel.background = element_rect(fill="white"),
    text = element_text(size=25,color="black",face="bold",hjust=0.5),
    axis.text = element_text(size=20,color="black",face="bold"),
    panel.border = element_rect(size=1,fill="NA"),
    strip.text = element_text(size=20,color="black",face="bold"),
    panel.spacing=unit(0,"cm"),
    plot.title = element_text(hjust=0.5, size=25)) 


#################################Boxplots faceted by phylum, with data source as x-axis########################################
#switch the 'datatype ==' statement to get figures for Raw, Trimmed-filtered and Trimmed-filtered-contigs data

PhylumDataKraken_Ankan %>%
  mutate(source = "Ankan") %>%
  rbind(PhylumDataKraken_Mahmud %>% mutate(source = "Mahmud")) %>%
  filter(name %in% PhylumDataKraken_Mahmud_summary$name,
         datatype == "Raw") %>%
  ggplot(aes(x = source, y = percent, fill = datatype)) +
  geom_boxplot(aes(x = source, y = percent, fill = source)) +
  facet_wrap(~name, scales = "free_y") +
  guides(fill = "none") +
  labs(y = "Relative Abundance %", x = "Source", title = paste("Relative abundances: Raw")) +
  ggpubr::stat_compare_means(method = "aov",size = 7, label.y.npc = 0.9, hjust = 0) +
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "white"),
    text = element_text(size = 25, color = "black", face = "bold", hjust = 0.5),
    axis.text = element_text(size = 20, color = "black", face = "bold"),
    panel.border = element_rect(size = 1, fill = "NA"),
    strip.text = element_text(size = 20, color = "black", face = "bold"),
    panel.spacing = unit(0, "cm"),
    plot.title = element_text(hjust = 0.5, size = 25)
  )

combined_PhylumDataKraken_summary<-PhylumDataKraken_Mahmud_summary%>%
  inner_join(PhylumDataKraken_Ankan_summary,by=c("datatype","name"))

############################################################Kracken-Bracken Species level plots#####################################################
SpeciesData_Ankan<-NULL
path<-"C:/Users/Ankan_Choudhury/Box/M_COLON/Stool Sample Sequencing/MRDNA/Ankan Output/BrackenOutputs/bracken-reports-trimmed-filtered-contigs"
filenames<-list.files(path,pattern="*report.txt")
folder<-"Trimmed-filtered-contigs"
for(each in filenames){
  print(paste("Processing:",each))
  file<-substring(each,1,(nchar(each)-11))
  file1<-tryCatch({read.table(paste(path,each,sep="/"),sep='\t',header = FALSE)%>%
      as.data.frame()%>%
      setNames(c("percent","clade_reads", "direct_reads", "rank_code", "tax_id", "name"))%>%
      mutate(name=trimws(name),
             sample=file,
             datatype=folder)%>%
      filter(rank_code=="S")%>%
      select(-c(2,4,5))%>%
      relocate(datatype,sample,name)
  }, warning = function(w) {
    message("Warning in file: ", each, "\n", w)  # Print file name if warning occurs
    return(NULL)  # Return NULL so the loop continues
  }, error = function(e) {
    message("Error in file: ", each, "\n", e)  # Print file name if error occurs
    return(NULL)
  })
  SpeciesData_Ankan<-SpeciesData_Ankan%>%
    rbind(file1)
}



SpeciesData_Ankan_summary<-SpeciesData_Ankan%>%
  group_by(datatype,sample)%>%
  mutate(percent=direct_reads/sum(direct_reads)*100,
         sum=sum(direct_reads))%>%
  group_by(name)%>%
  mutate(mean_percent=mean(percent))%>%
  group_by(datatype,sample)%>%
  mutate(name=if_else(mean_percent<1,"Below 1% mean abundance",name))%>%
  group_by(datatype,sample,name)%>%
  summarise(reads=sum(direct_reads))%>%
  mutate(percent=reads/sum(reads)*100)%>%
  group_by(datatype,name)%>%
  summarize(mean_percent=mean(percent))

SpeciesData_Ankan_summary%>%
  ggplot(aes(x=datatype,y=mean_percent,fill=reorder(name,-mean_percent)))+
  geom_bar(position="fill",stat="identity",color="black",width=1)+
  scale_y_continuous(labels = scales::percent_format())+
  labs(y="Relative Abundance",x="Data Type",fill="Species",title="Species Relative abundances: Ankan's Bracken Data")+
  guides(fill=guide_legend(ncol=2))+
  scale_x_discrete(expand = c(-2, -2),labels = c("Raw" = "Raw", "Trimmed-filtered" = "Trimmed\nfiltered", "Trimmed-filtered-contigs" = "Trimmed\nfiltered\ncontigs"))+
  scale_fill_manual(values=P12)+theme_classic()+
  coord_equal(7.5,expand = FALSE,xlim=c(0,4))+
  theme(
    panel.background = element_rect(fill="white"),
    text = element_text(size=22,color="black",face="bold"),
    panel.border = element_rect(size=1,fill="NA"),
    strip.text = element_text(size=22,color="black",face="bold"),
    legend.text = element_text(size=15,color="black",face="italic"),
    legend.margin = unit(0,"cm"),
    legend.title=element_text(size=20,face="bold",color="black"),
    legend.key.size= unit(1, 'cm'),
    plot.title = element_text(hjust=-1, size=20))

SpeciesData_Ankan%>%
  group_by(datatype,sample)%>%
  mutate(name=if_else(percent<5,"<5% abundances in the sample",name))%>%
  group_by(datatype,sample,name)%>%
  summarise(percent=sum(percent))%>%
  group_by(name)%>%
  mutate(mean_percent=mean(percent))%>%
  ggplot(aes(x=sample,y=percent,fill=reorder(name,-mean_percent)))+
  geom_bar(position="fill",stat="identity",color="black",width=1)+
  scale_y_continuous(labels = scales::percent_format())+
  labs(y="Relative Abundance",x=NULL,fill="Species",title="Species Relative abundances: Ankan's Bracken Data")+
  guides(fill=guide_legend(ncol=3))+
  scale_fill_manual(values=P12)+theme_classic()+
  ggh4x::facet_nested(~datatype)+
  #scale_x_discrete(expand = c(-2, -2),labels = c("Raw" = "Raw", "Trimmed-filtered" = "Trimmed\nfiltered", "Trimmed-filtered-contigs" = "Trimmed\nfiltered\ncontigs"))+
  #coord_equal(9,expand = FALSE,xlim=c(0,4))+
  theme(
    panel.background = element_rect(fill="white"),
    text = element_text(size=22,color="black",face="bold"),
    axis.text.x = element_blank(),
    panel.border = element_rect(size=1,fill="NA"),
    strip.text = element_text(size=15,color="black",face="bold"),
    legend.text = element_text(size=12,color="black",face="italic"),
    legend.margin = unit(0,"cm"),
    legend.title=element_text(size=20,face="bold",color="black"),
    legend.key.size= unit(0.75, 'cm'),
    plot.title = element_text(hjust=0, size=20),
    plot.margin = margin(0,0,0,0),
    panel.spacing = unit(0,"mm"))

#########################################################################Kaiju###################################################################################################

library(readxl)
kaijufiles<-list.files("Kaiju Taxonomy Abundance",pattern = "*.xlsx")
PhylumKaiju<-NULL
SpeciesKaiju<-NULL
for(each in kaijufiles){
  table<-read_excel(paste("Kaiju Taxonomy Abundance/",each,sep=""))%>%
    filter(Domain!="unclassified")%>%
    filter(Domain!="cannot be assigned to a (non-viral) species")
  Phylumtable<-table%>%
    select(-c(1,3:7))%>%
    pivot_longer(-Phylum)%>%
    group_by(name,Phylum)%>%
    summarize(reads=sum(value))%>%
    mutate(datatype=substr(each,12,nchar(each)-5))
  PhylumKaiju<-PhylumKaiju%>%
    rbind(Phylumtable)
  Speciestable<-table%>%
    select(-c(1:6))%>%
    pivot_longer(-Species)%>%
    group_by(name,Species)%>%
    summarize(reads=sum(value))%>%
    mutate(datatype=substr(each,12,nchar(each)-5))
  SpeciesKaiju<-SpeciesKaiju%>%
    rbind(Speciestable)
  
}

PhylumKaijuSummary<-PhylumKaiju%>%
  group_by(datatype,name)%>%
  mutate(percent=reads/sum(reads)*100,
         sum=sum(reads))%>%
  group_by(datatype,Phylum)%>%
  summarize(mean_percent=mean(percent))

PhylumKaijuSummary%>%
  ggplot(aes(x=datatype,y=mean_percent,fill=reorder(Phylum,-mean_percent)))+
  geom_bar(position="fill",stat="identity",color="black",width=1)+
  scale_y_continuous(labels = scales::percent_format())+
  labs(y="Relative Abundance",x="Data Type",fill="Species",title="Phylum Relative abundances: Ankan's Kaiju Data")+
  guides(fill=guide_legend(ncol=1))+
  scale_x_discrete(limits=c("reads","trimmed","filtered","trimmed-filtered"),
                   expand = c(-2, -2),
                   labels = c("reads" = "Raw", "trimmed-filtered" = "Trimmed\nfiltered", "trimmed" = "Trimmed","filtered"="Filtered"))+
  scale_fill_manual(values=P12)+theme_classic()+
  coord_equal(10,expand = FALSE,xlim=c(-0.5,5.5))+
  theme(
    panel.background = element_rect(fill="white"),
    text = element_text(size=22,color="black",face="bold"),
    panel.border = element_rect(size=1,fill="NA"),
    strip.text = element_text(size=22,color="black",face="bold"),
    legend.text = element_text(size=15,color="black",face="italic"),
    legend.margin = unit(0,"cm"),
    legend.title=element_text(size=20,face="bold",color="black"),
    legend.key.size= unit(1, 'cm'),
    plot.title = element_text(hjust=-1, size=20))  

PhylumKaijuSummary<-PhylumKaiju%>%
  mutate(Phylum=if_else(Phylum=="belong to a (non-viral) species with less than 0.01% of all reads","<0.01% of all reads",Phylum))%>%
  group_by(datatype,name)%>%
  mutate(percent=reads/sum(reads)*100,
         sum=sum(reads))%>%
  group_by(datatype,Phylum)%>%
  summarize(mean_percent=mean(percent))

PhylumKaijuSummary%>%
  ggplot(aes(x=datatype,y=mean_percent,fill=reorder(Phylum,-mean_percent)))+
  geom_bar(position="fill",stat="identity",color="black",width=1)+
  scale_y_continuous(labels = scales::percent_format())+
  labs(y="Relative Abundance",x="Data Type",fill="Species",title="Phylum Relative abundances: Ankan's Kaiju Data")+
  guides(fill=guide_legend(ncol=1))+
  scale_x_discrete(limits=c("reads","trimmed","filtered","trimmed-filtered"),
                   expand = c(-1, -1),
                   labels = c("reads" = "Raw", "trimmed-filtered" = "Trimmed\nfiltered", "trimmed" = "Trimmed","filtered"="Filtered"))+
  scale_fill_manual(values=P12)+theme_classic()+
  coord_equal(10,expand = FALSE,xlim=c(-0.5,5.5))+
  theme(
    panel.background = element_rect(fill="white"),
    text = element_text(size=22,color="black",face="bold"),
    panel.border = element_rect(size=1,fill="NA"),
    strip.text = element_text(size=22,color="black",face="bold"),
    legend.text = element_text(size=15,color="black",face="italic"),
    legend.margin = unit(0,"cm"),
    legend.title=element_text(size=20,face="bold",color="black"),
    legend.key.size= unit(1, 'cm'),
    plot.title = element_text(hjust=-1, size=20)) 

PhylumKaiju%>%
  mutate(Phylum=if_else(Phylum=="belong to a (non-viral) species with less than 0.01% of all reads","<0.01% of all reads",Phylum))%>%
  group_by(datatype,name)%>%
  mutate(percent=reads/sum(reads)*100,
         sum=sum(reads))%>%
  group_by(datatype,Phylum)%>%
  mutate(mean_percent=mean(percent))%>%  
  # group_by(datatype,sample)%>%
  # # mutate(name=if_else(percent<0.5,"<0.5% abundances in the sample",name))%>%
  # # group_by(datatype,sample,name)%>%
  # # summarise(percent=sum(percent))%>%
  ggplot(aes(x=name,y=percent,fill=reorder(Phylum,-mean_percent)))+
  geom_bar(position="fill",stat="identity",color="black",width=1)+
  scale_y_continuous(labels = scales::percent_format())+
  labs(y="Relative Abundance",x=NULL,fill="Phylum",title="Phylum Relative abundances: Ankan's Kaiju Data")+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(values=P12)+theme_classic()+
  ggh4x::facet_nested(~factor(datatype, levels = c("reads", "trimmed", "filtered", "trimmed-filtered"),
                              labels = c("Raw", "Trimmed", "Filtered", "Trimmed\nfiltered")))+
  theme(
    panel.background = element_rect(fill="white"),
    text = element_text(size=22,color="black",face="bold"),
    axis.text.x = element_blank(),
    panel.border = element_rect(size=1,fill="NA"),
    strip.text = element_text(size=22,color="black",face="bold"),
    legend.text = element_text(size=15,color="black",face="italic"),
    legend.margin = unit(0,"cm"),
    legend.title=element_text(size=20,face="bold",color="black"),
    legend.key.size= unit(1, 'cm'),
    plot.title = element_text(hjust=0, size=20),
    plot.margin = margin(0,0,0,0),
    panel.spacing = unit(0,"mm"))

SpeciesKaijuSummary<-SpeciesKaiju%>%
  mutate(Species=if_else(Species=="belong to a (non-viral) species with less than 0.01% of all reads","<0.01% of all reads",Species))%>%
  group_by(datatype,name)%>%
  mutate(percent=reads/sum(reads)*100,
         sum=sum(reads))%>%
  group_by(datatype,Species)%>%
  summarize(mean_percent=mean(percent))%>%
  mutate(Species=if_else(mean_percent<0.75,"Below 0.75% mean abundance",Species))%>%
  group_by(datatype,Species)%>%
  summarise(mean_percent=sum(mean_percent))%>%
  mutate(Species=str_replace_all(Species,"_"," "))

SpeciesKaijuSummary%>%
  ggplot(aes(x=datatype,y=mean_percent,fill=reorder(Species,-mean_percent)))+
  geom_bar(position="fill",stat="identity",color="black",width=1)+
  scale_y_continuous(labels = scales::percent_format())+
  labs(y="Relative Abundance",x="Data Type",fill="Species",title="Species Relative abundances: Ankan's Kaiju Data")+
  guides(fill=guide_legend(ncol=2))+
  scale_x_discrete(limits=c("reads","trimmed","filtered","trimmed-filtered"),
                   expand = c(-1, -1),
                   labels = c("reads" = "Raw", "trimmed-filtered" = "Trimmed\nfiltered", "trimmed" = "Trimmed","filtered"="Filtered"))+
  scale_fill_manual(values=P12)+theme_classic()+
  coord_equal(10,expand = FALSE,xlim=c(-0.5,5.5))+
  theme(
    panel.background = element_rect(fill="white"),
    text = element_text(size=22,color="black",face="bold"),
    axis.text = element_text(size=15,color="black",face="bold"),
    panel.border = element_rect(size=1,fill="NA"),
    strip.text = element_text(size=22,color="black",face="bold"),
    legend.text = element_text(size=15,color="black",face="italic"),
    legend.margin = unit(0,"cm"),
    legend.title=element_text(size=20,face="bold",color="black"),
    legend.key.size= unit(1, 'cm'),
    plot.title = element_text(hjust=0, size=20)) 


SpeciesKaiju%>%
  mutate(Species=if_else(Species=="belong to a (non-viral) species with less than 0.01% of all reads","<0.01% of all reads",Species))%>%
  group_by(datatype,name)%>%
  mutate(percent=reads/sum(reads)*100,
         sum=sum(reads))%>%
  # group_by(datatype,Species)%>%
  # mutate(mean_percent=mean(percent))%>%  
  group_by(datatype,name)%>%
  mutate(Species=if_else(percent<5,"<5% abundances in the sample",Species))%>%
  group_by(datatype,name,Species)%>%
  summarise(percent=sum(percent))%>%
  ggplot(aes(x=name,y=percent,fill=reorder(Species,-percent)))+
  geom_bar(position="fill",stat="identity",color="black",width=1)+
  scale_y_continuous(labels = scales::percent_format())+
  labs(y="Relative Abundance",x=NULL,fill="Species",title="Species Relative abundances: Ankan's Kaiju Data")+
  guides(fill=guide_legend(ncol=2))+
  scale_fill_manual(values=P12)+theme_classic()+
  ggh4x::facet_nested(~factor(datatype, levels = c("reads", "trimmed", "filtered", "trimmed-filtered"),
                              labels = c("Raw", "Trimmed", "Filtered", "Trimmed\nfiltered")))+
  theme(
    panel.background = element_rect(fill="white"),
    text = element_text(size=22,color="black",face="bold"),
    axis.text.x = element_blank(),
    panel.border = element_rect(size=1,fill="NA"),
    strip.text = element_text(size=22,color="black",face="bold"),
    legend.text = element_text(size=12,color="black",face="italic"),
    legend.margin = unit(0,"cm"),
    legend.title=element_text(size=20,face="bold",color="black"),
    legend.key.size= unit(0.75, 'cm'),
    plot.title = element_text(hjust=0, size=20),
    plot.margin = margin(0,0,0,0),
    panel.spacing = unit(0,"mm"))


##############################################################################Metaphlan######################################################################################
MetaphlanOutputs<-dir("MetaphlanOutput")
PhylumMetaphlan<-NULL
SpeciesMetaphlan<-NULL
for(folder in MetaphlanOutputs){
  for(file in list.files(paste("MetaphlanOutput",folder,sep="/"))){
    print(paste("Processing: MetaphlanOutput",folder,file,sep="/"))
    content<-readLines(paste("MetaphlanOutput",folder,file,sep="/"))[-c(1:4)]
    table<-read.table(text=content,header=FALSE,sep="\t")%>%
      select(-c(2,4))%>%
      mutate(V1=str_replace_all(V1,".__",""),
             V1 = str_replace_all(V1, "\\|", ";")) %>%
      separate(V1,c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),sep=";")
    Phylum<-table%>%
      filter(!(is.na(Phylum))&is.na(Class))%>%
      select(2,8)%>%
      rename(percent=V3)%>%
      mutate(Sample=substr(file,1,nchar(file)-11),
             datatype=folder)
    PhylumMetaphlan<-PhylumMetaphlan%>%
      rbind(Phylum)
    Species<-table%>%
      filter(!(is.na(Species)))%>%
      select(7,8)%>%
      unique()%>%
      group_by(Species)%>%
      summarise(sum=sum(V3))%>%
      rename(percent=sum)%>%
      mutate(Sample=substr(file,1,nchar(file)-11),
             datatype=folder)
    SpeciesMetaphlan<-SpeciesMetaphlan%>%
      rbind(Species)
    
  }
}

PhylumMetaphlan_summary<-PhylumMetaphlan%>%
  group_by(datatype,Phylum)%>%
  summarize(mean_percent=mean(percent))%>%
  mutate(Phylum=if_else(Phylum=="Bacteria_unclassified","Unclassified Bacteria",Phylum))
PhylumMetaphlan_summary%>%
  ggplot(aes(x=datatype,y=mean_percent,fill=reorder(Phylum,-mean_percent)))+
  geom_bar(position="fill",stat="identity",color="black",width=1)+
  scale_y_continuous(labels = scales::percent_format())+
  labs(y="Relative Abundance",x="Data Type",fill="Species",title="Phylum Relative abundances: Ankan's Metaphlan Data")+
  guides(fill=guide_legend(ncol=1))+
  scale_x_discrete(limits=c("raw-paired","trimmed-paired","filtered-paired","trimmed-filtered-paired"),
                   expand = c(-1, -1),
                   labels = c("raw-paired" = "Raw", "trimmed-filtered-paired" = "Trimmed\nfiltered", "trimmed-paired" = "Trimmed","filtered-paired"="Filtered"))+
  scale_fill_manual(values=P12)+theme_classic()+
  coord_equal(10,expand = FALSE,xlim=c(-0.5,5.5))+
  theme(
    panel.background = element_rect(fill="white"),
    text = element_text(size=22,color="black",face="bold"),
    axis.text = element_text(size=15,color="black",face="bold"),
    panel.border = element_rect(size=1,fill="NA"),
    strip.text = element_text(size=22,color="black",face="bold"),
    legend.text = element_text(size=15,color="black",face="italic"),
    legend.margin = unit(0,"cm"),
    legend.title=element_text(size=20,face="bold",color="black"),
    legend.key.size= unit(1, 'cm'),
    plot.title = element_text(hjust=0, size=20)) 

PhylumMetaphlan%>%
  group_by(datatype,Phylum)%>%
  mutate(mean_percent=mean(percent))%>%  
  # group_by(datatype,sample)%>%
  # # mutate(name=if_else(percent<0.5,"<0.5% abundances in the sample",name))%>%
  # # group_by(datatype,sample,name)%>%
  # # summarise(percent=sum(percent))%>%
  ggplot(aes(x=Sample,y=percent,fill=reorder(Phylum,-mean_percent)))+
  geom_bar(position="fill",stat="identity",color="black",width=1)+
  scale_y_continuous(labels = scales::percent_format())+
  labs(y="Relative Abundance",x=NULL,fill="Phylum",title="Phylum Relative abundances: Ankan's Metaphlan Data")+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(values=P12)+theme_classic()+
  ggh4x::facet_nested(~factor(datatype, levels = c("raw-paired", "trimmed-paired", "filtered-paired", "trimmed-filtered-paired"),
                              labels = c("Raw", "Trimmed", "Filtered", "Trimmed\nfiltered")))+
  theme(
    panel.background = element_rect(fill="white"),
    text = element_text(size=22,color="black",face="bold"),
    axis.text.x = element_blank(),
    panel.border = element_rect(size=1,fill="NA"),
    strip.text = element_text(size=22,color="black",face="bold"),
    legend.text = element_text(size=15,color="black",face="italic"),
    legend.margin = unit(0,"cm"),
    legend.title=element_text(size=20,face="bold",color="black"),
    legend.key.size= unit(1, 'cm'),
    plot.title = element_text(hjust=0, size=20),
    plot.margin = margin(0,0,0,0),
    panel.spacing = unit(0,"mm"))

SpeciesMetaphlan_summary<-SpeciesMetaphlan%>%
  # group_by(datatype,Sample)%>%
  # mutate(Species=if_else(percent<10,"Below 10% abundance",Species))%>%
  # group_by(datatype,Sample,Species)%>%
  # summarise(percent=sum(percent))%>%
  group_by(datatype,Species)%>%
  summarize(mean_percent=mean(percent))%>%
  mutate(Species=if_else(mean_percent<2.5,"Below 2.5% mean abundance",Species))%>%
  group_by(datatype,Species)%>%
  summarise(mean_percent=sum(mean_percent))%>%
  mutate(Species=str_replace_all(Species,"_"," "))


SpeciesMetaphlan_summary%>%
  ggplot(aes(x=datatype,y=mean_percent,fill=reorder(Species,-mean_percent)))+
  geom_bar(position="fill",stat="identity",color="black",width=1)+
  scale_y_continuous(labels = scales::percent_format())+
  labs(y="Relative Abundance",x="Data Type",fill="Species",title="Species Relative abundances: Ankan's Metaphlan Data")+
  guides(fill=guide_legend(ncol=2))+
  scale_x_discrete(limits=c("raw-paired","trimmed-paired","filtered-paired","trimmed-filtered-paired"),
                   expand = c(-1, -1),
                   labels = c("raw-paired" = "Raw", "trimmed-filtered-paired" = "Trimmed\nfiltered", "trimmed-paired" = "Trimmed","filtered-paired"="Filtered"))+
  scale_fill_manual(values=P12)+theme_classic()+
  coord_equal(10,expand = FALSE,xlim=c(-0.5,5.5))+
  theme(
    panel.background = element_rect(fill="white"),
    text = element_text(size=22,color="black",face="bold"),
    axis.text = element_text(size=15,color="black",face="bold"),
    panel.border = element_rect(size=1,fill="NA"),
    strip.text = element_text(size=22,color="black",face="bold"),
    legend.text = element_text(size=15,color="black",face="italic"),
    legend.margin = unit(0,"cm"),
    legend.title=element_text(size=20,face="bold",color="black"),
    legend.key.size= unit(1, 'cm'),
    plot.title = element_text(hjust=0, size=20)) 

SpeciesMetaphlan%>%
  mutate(Species=str_replace_all(Species,"_"," "))%>%
  group_by(datatype,Species)%>%
  mutate(mean_percent=mean(percent))%>%
  group_by(datatype,Sample)%>%
  mutate(Species=if_else(percent<5,"<5% abundances in the sample",Species))%>%
  group_by(datatype,Sample,Species)%>%
  summarise(percent=sum(percent))%>%
  ggplot(aes(x=Sample,y=percent,fill=reorder(Species,-percent)))+
  geom_bar(position="fill",stat="identity",color="black",width=1)+
  scale_y_continuous(labels = scales::percent_format())+
  labs(y="Relative Abundance",x=NULL,fill="Species",title="Species Relative abundances: Ankan's Metaphlan Data")+
  guides(fill=guide_legend(ncol=3))+
  scale_fill_manual(values=P12)+theme_classic()+
  ggh4x::facet_nested(~factor(datatype, levels = c("raw-paired", "trimmed-paired", "filtered-paired", "trimmed-filtered-paired"),
                              labels = c("Raw", "Trimmed", "Filtered", "Trimmed\nfiltered")))+
  theme(
    panel.background = element_rect(fill="white"),
    text = element_text(size=22,color="black",face="bold"),
    axis.text.x = element_blank(),
    panel.border = element_rect(size=1,fill="NA"),
    strip.text = element_text(size=22,color="black",face="bold"),
    legend.text = element_text(size=12,color="black",face="italic"),
    legend.margin = unit(0,"cm"),
    legend.title=element_text(size=20,face="bold",color="black"),
    legend.key.size= unit(0.75, 'cm'),
    plot.title = element_text(hjust=0, size=20),
    plot.margin = margin(0,0,0,0),
    panel.spacing = unit(0,"mm"))

SpeciesMetaphlan%>%
  group_by(datatype,Species)%>%
  summarize(mean_percent=mean(percent))%>%
  ggplot(aes(x=datatype,y=reorder(Species,mean_percent),fill=log2(mean_percent)))+
  geom_tile()+
  # facet_grid(~datatype)+
  coord_equal(0.01,expand = FALSE)+
  theme_classic()+
  theme(axis.text.y=element_blank())


##############################################################################Metaphlan Phyloseq######################################################################################
MetaphlanOutputs<-dir("MetaphlanOutput")
dir.create("MetaphlanPhyloseq")
MetaphlanTree<-read_tree("20201013_mpav30_speciesMod.nwk")
new_names <- MetaphlanTree$tip.label
new_names<-new_names%>%
  str_replace_all("s__","")
MetaphlanTree$tip.label<-new_names
for(folder in MetaphlanOutputs){
  Alltaxonomy<-NULL
  AllOTU<-NULL
  for(file in list.files(paste("MetaphlanOutput",folder,sep="/"))){
    content<-readLines(paste("MetaphlanOutput",folder,file,sep="/"))[-c(1:4)]
    id<-substring(basename(file),1,nchar(basename(file))-11)
    totalreads<-readLines(paste("MetaphlanOutput",folder,file,sep="/"))[3]
    totalreads<-as.numeric(gsub("^#(\\d+).*", "\\1", totalreads))
    table<-read.table(text=content,header=FALSE,sep="\t")%>%
      select(-4)%>%
      mutate(V1=str_replace_all(V1,".__",""),
             V1 = str_replace_all(V1, "\\|", ";"),
             reads=round(V3/100*totalreads)) %>%
      separate(V1,c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),sep=";")%>%
      separate(V2,c("KID","PID","CID","OID","FID","GID","SID"),sep="\\|")%>%
      unique()%>%
      filter(!(is.na(Species)))%>%
      mutate(PID=if_else(PID==""|is.na(PID),KID,PID),
             CID=if_else(CID==""|is.na(CID),PID,CID),
             OID=if_else(OID==""|is.na(OID),CID,OID),
             FID=if_else(FID==""|is.na(FID),OID,FID),
             GID=if_else(GID==""|is.na(GID),FID,GID),
             SID=if_else(SID==""|is.na(SID),GID,SID))
    taxonomy<-table[,c(1:7)]
    OTU<-table[,c(7,16)]%>%
      unique()%>%
      mutate(Sample=id)
    Alltaxonomy<-rbind(Alltaxonomy,taxonomy)
    AllOTU<-rbind(AllOTU,OTU)
    print(paste("Processed: MetaphlanOutput",folder,file,sep="/"))
  }
  
  taxa_table<-Alltaxonomy%>%
    arrange("Species")%>%
    unique()
  rownames(taxa_table)<-taxa_table$Species
  otu_table<-AllOTU%>%
    group_by(Sample,Species)%>%
    summarize(reads=sum(reads))%>%
    pivot_wider(id_cols = "Species",values_from = "reads",names_from = "Sample",values_fill = 0)%>%
    column_to_rownames("Species")
  physeqMetaphlan<-phyloseq(otu_table(as.matrix(otu_table),taxa_are_rows = TRUE),tax_table(as.matrix(taxa_table)))
  physeqMetaphlan<-merge_phyloseq(physeqMetaphlan,phy_tree(MetaphlanTree))
  saveRDS(physeqMetaphlan,paste("MetaphlanPhyloseq/",folder,".rds",sep=""))
}

#######################################Humann Phyloseq###############################################################################
HumannOutput<-dir("HumannOutput")
for (folder in HumannOutput){
  AllCoverage<-data.frame()
  AllAbundance<-data.frame()
  counter=1
  for(file in list.files(paste("HumannOutput",folder,sep="/"),pattern="*_pathabundance.tsv")){
    id<-substring(basename(file),1,nchar(basename(file))-18)
    abundance<-read.delim(paste("HumannOutput",folder,file,sep="/"),sep='\t',header = TRUE,comment.char = "")
    df<-colnames(abundance)%>%
      str_replace_all("X\\.\\.","")%>%
      str_replace_all("X","")%>%
      str_replace_all("\\.","_")%>%
      str_replace_all("_Abundance","")
    colnames(abundance)<-df
    coverage<-read.delim(paste("HumannOutput/",folder,"/",id,"_pathcoverage.tsv",sep=""),sep='\t',header = TRUE,comment.char = "")
    df<-colnames(coverage)%>%
  str_replace_all("X\\.\\.","")%>%
  str_replace_all("X","")%>%
  str_replace_all("\\.","-")%>%
  str_replace_all("_Coverage","")
    colnames(coverage)<-df
    if (counter==1){
      AllAbundance<-abundance}
    else{AllAbundance<-left_join(AllAbundance,abundance)}
    if(counter==1){
    AllCoverage<-coverage}
    else{AllCoverage<-left_join(AllCoverage,coverage)}
    counter=counter+1
      print(paste("Processed",folder,":",id))
  }
  
  AllAbundance<-AllAbundance%>%
    column_to_rownames("Pathway")%>%
    as.matrix()
  AllCoverage<-AllCoverage%>%
    column_to_rownames("Pathway")%>%
    as.matrix()
  physeqcoverage<-phyloseq(otu_table(AllCoverage,taxa_are_rows=TRUE))
  physeqabundance<-phyloseq(otu_table(AllAbundance,taxa_are_rows=TRUE))
  saveRDS(physeqabundance,paste("HumannPhyloseq/",folder,"_pathwayabundance.rds",sep=""))
  saveRDS(physeqcoverage,paste("HumannPhyloseq/",folder,"_pathwaycoverage.rds",sep=""))
}

for (folder in HumannOutput){
  AllRPK<-data.frame()
  counter=1
  for(file in list.files(paste("HumannOutput",folder,sep="/"),pattern="*_genefamilies.tsv")){
    id<-substring(basename(file),1,nchar(basename(file))-17)
    RPK<-read.delim(paste("HumannOutput",folder,file,sep="/"),sep='\t',header = TRUE,comment.char = "")
    df<-colnames(RPK)%>%
      str_replace_all("X\\.\\.","")%>%
      str_replace_all("X","")%>%
      str_replace_all("\\.","-")%>%
      str_replace_all("_Abundance-RPKs","")
    colnames(RPK)<-df
    if (counter==1){
      AllRPK<-RPK}
    else{AllRPK<-left_join(AllRPK,RPK)}
    counter=counter+1
    print(paste("Processed",folder,":",id))
  }
  
  AllRPK<-AllRPK%>%
    column_to_rownames("Gene-Family")%>%
    as.matrix()
  physeqRPK<-phyloseq(otu_table(AllRPK,taxa_are_rows=TRUE))
  saveRDS(physeqRPK,paste("HumannPhyloseq/",folder,"_genefamilies.rds",sep=""))
  rm(AllRPK)
}


