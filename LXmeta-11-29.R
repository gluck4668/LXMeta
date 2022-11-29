
LXmeta <- function(file_data,species){

# install the dependency packages using a function [metanr_packages()]
  metanr_packages <- function(){
    metr_pkgs <- c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz",
                   "preprocessCore", "genefilter", "sva", "limma", "KEGGgraph",
                   "siggenes","BiocParallel", "MSnbase", "multtest", "RBGL", "edgeR",
                   "fgsea", "crmn")
    list_installed <- installed.packages()
    new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
    if(length(new_pkgs)!=0){if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
      library(BiocManager)
      BiocManager::install(new_pkgs,force=TRUE,quietly = TRUE)
      print(c(new_pkgs, " packages added..."))
    }

    if((length(new_pkgs)<1)){
      print("No new BiocManager packages added...")
    }
  }


  metanr_packages()

# list all the packages that have been installed

  R_packs_install <- function(){

    metr_pkgs <- c("devtools", "httr", "openxlsx", "magrittr", "dplyr",
                 "psych", "ggplot2", "ggrepel", "RColorBrewer", "ggthemes",
                 "tidyverse","roxygen2", "XML", "RCurl", "curl", "stringr",
                 "plyr", "momr","patchwork","ggpubr","scales","ggplotify")

   list_installed <- installed.packages()

   new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))

   if(length(new_pkgs)!=0){install.packages(new_pkgs,force=TRUE,quietly = TRUE)
                         print(c(new_pkgs, " packages added..."))}

  if((length(new_pkgs)<1)){print("No new dependency packages added...")}

 }

  R_packs_install()

  library(devtools)
  if (!requireNamespace("MetaboAnalystR", quietly = TRUE))
  install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = TRUE, build_manual =T)


  packages <- c("MetaboAnalystR", "httr", "openxlsx", "magrittr", "dplyr",
                "psych", "ggplot2", "ggrepel", "RColorBrewer", "ggthemes",
                "tidyverse","roxygen2", "XML", "RCurl", "curl", "stringr",
                "plyr", "momr","patchwork","ggpubr","scales","ggplotify","conflicted")

  for(i in packages){
    library(i, character.only = T)
  }

  rm(i)


  #-------------------------------------------------------------------


# To select a library.
  #(1) Chemcial structures:
      # "main_class"：464 main chemical class metabolite sets or lipid sets.
      # "super_class"：35 super chemical class metabolite sets or lipid sets.
      # "sub_class": 1072 sub chemical class metabolite sets or lipid sets.
  #(2) Pathway based:
      # "smpdb_pathway":99 metabolite sets based on normal human metabolic pathways.
      # "kegg_pathway": 84 metabolite sets based on KEGG human metabolic pathways (Oct. 2019).
      # "drug": 461 metabolite sets based on drug pathways from SMPDB.

  libraries="main_class"

# reading the data file
  meta_set <- read.xlsx(file_data)

# To judge the type of the data, "HMDB" or "name"
  if(substr(meta_set[1,1],1,4)=="HMDB")
    inputType = "hmdb" else
      inputType = "name"

# translate the first column as a character list
  m_set <- t(meta_set[,1]) %>% as.character()

# perform the packages below to get the mapping metabolites from the KEGG database.
  m_obj<-InitDataObjects("conc", "msetora", FALSE)
  m_obj<-Setup.MapData(m_obj, m_set)
  m_obj<-CrossReferencing(m_obj, inputType)
  m_obj<-CreateMappingResultTable(m_obj)
  m_Filter<-SetMetabolomeFilter(m_obj, F)

# screen out the mapping result from the m_Filter list.
  mapping_result_all <- m_Filter$dataSet$map.table %>% data.frame()

# screen out the mapping result without the KEGG ID "NA".
  mapping_result_all$KEGG[which(mapping_result_all$KEGG=="NA")]=NA
  mapping_result <- dplyr::filter(mapping_result_all,!is.na(KEGG))

# create a folder to save the files of the analysis results.
  if(dir.exists("analysis results")==FALSE)
    dir.create("analysis results")

  write.xlsx(mapping_result_all,"analysis results/mapping_result_all.xlsx")
  write.xlsx(mapping_result,"analysis results/mapping_result.xlsx")

# name the libraries as a file with a type of ".qs" . In this case, it is "main_class.qs".
  database <- paste0(libraries,".qs")

# download the library file from the website
  url <- paste0("https://www.metaboanalyst.ca/resources/libs/msets/",database)
  curl_download(url, database)

#mSet0<-SetCurrentMsetLib(m_Filter, libraries, 2)

# Set Current Metabolite set library
#---------------------------- begin
  on_public_web <- FALSE  # only TRUE when on metaboanalyst web server
                           # note, this is usually used at the end of a function
                           # for local, return itself; for web, push to global environment

  set_mSet <- function(mSetObj=m_Filter){
    if(on_public_web){mSet <<- mSetObj
                      return (1)}
    return(mSetObj)
    }


  get_mSet <- function(mSetObj=m_Filter){
    if(on_public_web){
      return(mSet)
    }else{
      return(mSetObj)
    }
   }


    mSetObj=m_Filter
    libname=libraries
    excludeNum=2

    mSetObj <- get_mSet(mSetObj)
#---------------
  if(!on_public_web & grepl("kegg", libname)){ # api only for KEGG msets
        mSetObj$api$libname <- libname
        mSetObj$api$excludeNum = excludeNum
        mSetObj$analSet$msetlibname <- libname
        return(set_mSet(mSetObj))
       }

  current.msetlib <- qs::qread(database)
  mSetObj$analSet$msetlibname <- libname

  # create a named list, use the ids for list names
  ms_list <- strsplit(current.msetlib[,3],"; ", fixed=TRUE)
  names(ms_list) <- current.msetlib[,2]

#--------------

  #cmpd.count <- lapply(ms_list, length)
  #sel.inx <- cmpd.count >= excludeNum
  #ms_list <- ms_list[sel.inx]

#--------------

# total uniq cmpds in the mset lib
mSetObj$dataSet$uniq.count <- length(unique(unlist(ms_list, use.names = FALSE)))

# update current.mset and push to global env
current.msetlib$member <- ms_list
current.msetlib <<- current.msetlib

mSet0 <- set_mSet(mSetObj)

#---------------------------- end

  mSet<-CalculateHyperScore(mSet0)

  file.remove("name_map.csv")
  file.rename("msea_ora_result.csv","analysis results/enrich_result.csv")

  enrich_df <- read.csv( "analysis results/enrich_result.csv",header=T)

  enrich_df$Enrichment_Ratio <- enrich_df$hits/enrich_df$expected

  colnames(enrich_df) <- c("Metabolite_Set","Total","Expect","Hits","Pvalue","Holm p","FDR","Enrichment_Ratio")

  enrich_df <- data.frame(enrich_df)

#-------------geom_pointr-------------------#
  title_text_size <- case_when(nrow(enrich_df)>30 ~12,
                               nrow(enrich_df)>20 ~11,
                               TRUE ~10)

  xy_text_size <- case_when(   nrow(enrich_df)>30 ~9,
                               nrow(enrich_df)>20 ~10,
                               TRUE ~10)

  mytheme<-theme_bw()+
    theme(text=element_text(family = "sans",colour ="black",face="bold",size =title_text_size),
          panel.grid = element_blank(),
          panel.border = element_rect(linewidth = 0.8,colour = "gray30"),
          axis.line = element_blank(),
          axis.ticks = element_line(linewidth = 0.6,colour = "gray30"),
          axis.ticks.length = unit(1.5,units = "mm"))+
    theme(plot.title = element_text(hjust = 0.5))


  xytheme <-theme(axis.text.x = element_text(face="bold",color="black",size=xy_text_size,angle =0,hjust=1))+
    theme(axis.text.y = element_text(face="bold",color="black",size=xy_text_size))+
    theme(legend.text=element_text(face="bold",color="black",size=10))

  grid_theme <- theme(panel.grid =element_line(colour="#dcdcdc",linewidth=0.2,linetype = "dashed"))

  rev_set <- rev(enrich_df$Metabolite_Set)

  enrich_df$term <- factor(enrich_df$Metabolite_Set,levels=rev_set)

  x_max <- max(-log10(enrich_df$Pvalue))
  y_max <- nrow(enrich_df)+1

if(x_max<1)
  my_color <- scale_color_gradient2(midpoint = 0,low = "#33ffff",mid = "#ffff00",high = "#ff0000") else
    my_color <- scale_color_gradient2(midpoint = 0,low = "#33ffff",mid = "#f08080",high = "#ff0000")

 enrich_point <- ggplot(enrich_df,aes(x=-log10(Pvalue),y=term))+
    geom_point(aes(color=-log10(Pvalue),size=Enrichment_Ratio))+
    my_color+
    #scale_color_gradient2(midpoint = mid_point,low = "#33ffff",mid = "#ffff00",high = "#ff0000")+
    labs(x = '-log10(Pvalue)', y = '',title="Overview of Enriched Metabolite Sets")+
    coord_cartesian(xlim=c(0, x_max),ylim=c(0,y_max),expand = T)+
    mytheme+xytheme+grid_theme+
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

enrich_point

#-------------geom_bar-------------------#
  legend_theme <- theme(
    legend.title = element_blank(),
    legend.position = "none"
    #legend.text = element_text(size = 10, face = "bold"),
    #legend.direction = "horizontal",
    #legend.position = c(0.5,0.9),
    #legend.background = element_blank()
  )

  text_theme<-theme_bw()+
    theme(panel.grid=element_blank())+
    theme(text=element_text(family = "sans",colour ="black",face="bold",size =title_text_size))+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text.x = element_text(face="bold",color="black",size=xy_text_size,angle =0,hjust=1))+
    theme(axis.text.y = element_text(face="bold",color="black",size=xy_text_size))


  cols<-brewer.pal(3,"YlOrRd")
  pal<-colorRampPalette(cols)
  mycolors<-pal(nrow(enrich_df))

  x_m <- max(enrich_df$Enrichment_Ratio)
  x_expand <- case_when(x_m>100 ~x_m*0.005,
                        x_m>50 ~x_m*0.001,
                        TRUE ~0.0001)
  x_lim <- x_m*1.05

  #enrich_data <- filter(enrich_df,Enrichment_Ratio>1)

  enrich_bar <- ggplot(enrich_df, aes(x =Enrichment_Ratio, y =term,fill=term))+
    geom_bar(position = "dodge",stat = "identity",width = 0.8)+
    scale_fill_manual(values=mycolors)+
    scale_x_continuous(expand = c(0, x_expand))+
    labs(x = 'Enrichment Ratio', y = '',title="Overview of Enriched Metabolite Sets")+
    text_theme+legend_theme+grid_theme+
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

  enrich_bar

#-----------------pie-------------------#

enrich_df_pie <- enrich_df[order(enrich_df$Hits,decreasing = T),]
sum_Hits <- sum(enrich_df_pie$Hits)
enrich_df_pie$percent <- round(100*enrich_df_pie$Hits/sum_Hits,2)
enrich_df_pie$percent100 <- paste0(enrich_df_pie$percent,"%")

enrich_df_pie <- dplyr::filter(enrich_df_pie,enrich_df_pie$percent>1)

enrich_df_pie$term <- paste0(enrich_df_pie$term," (",enrich_df_pie$percent100,")")

enrich_df_pie$term <- factor(enrich_df_pie$term,levels=enrich_df_pie$term,ordered = T)

col_table <- c("#ff6347","#6b8e23","#daa520","#87cefa","#a0522d","#98fb98","#ffff00","#556b2f","#adff2f","#cc0000",
               "#7fffd4","#008b8b","#ba55d3","#0000cd","#ffb6c1","#7b68ee","#00bfff","#5f9ea0","#e0ffff","#32cd32",
               "#ffd700","#330033","#ee82ee","#9370db","#7cfc00","#afeeee","#bc8f8f","#00ff7f","#2e8b57","#b8860b")

my_col <- col_table[c(1:nrow(enrich_df_pie))]

leng_text <- case_when(nrow(enrich_df_pie)>20 ~9,
                       nrow(enrich_df_pie)>10 ~10,
                       TRUE ~11)
leng_theme <- theme(legend.position="right",legend.title = element_blank(),
                    legend.text = element_text(colour="black",size=leng_text))


title_pie <- case_when(nrow(enrich_df_pie)>30 ~14,
                             nrow(enrich_df)>20 ~13,
                             TRUE ~12)

xy_pie <- case_when(nrow(enrich_df_pie)>30 ~9,
                         nrow(enrich_df)>20 ~10,
                         TRUE ~10)

pie_theme<-theme_bw()+
     theme(panel.grid=element_blank())+
     theme(text=element_text(family = "sans",colour ="black",face="bold",size =title_pie))+
     theme(plot.title = element_text(hjust = 0.9))+
     theme(axis.text.x = element_text(face="bold",color="black",size=xy_pie,angle =0,hjust=1))+
     theme(axis.text.y = element_text(face="bold",color="black",size=xy_pie))

blank_theme_pie <-theme(
    #axis.title.x = element_blank(),
    #axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank()
  )


enrich_pie <- ggplot(enrich_df_pie, aes(x ="", y =Hits,fill=term))+
              geom_bar(stat = "identity",width = 1,color="white")+
              coord_polar(theta="y",start = 0,direction = -1)+ # direction1, clockwise; -1, anticlockwise
              labs(x = '', y = '',title="The main enriched metabolite sets")+
              scale_fill_manual(values=my_col)+
              pie_theme +
              blank_theme_pie +
              leng_theme +
              theme(plot.title = element_text(hjust = 0.5))
              theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

enrich_pie

pic_point <- paste0("analysis results/Metabolites Enrichment_point (",libraries,").png")
pic_bar <- paste0("analysis results/Metabolites Enrichment_bar (",libraries,").png")
pic_pie <- paste0("analysis results/Metabolites Enrichment_pie (",libraries,").png")


ggsave(pic_point, enrich_point, width=1400, height =1200, dpi=150,units = "px")
ggsave(pic_bar, enrich_bar, width=1400, height =1200, dpi=150,units = "px")
ggsave(pic_pie, enrich_pie, width=1400, height =1200, dpi=150,units = "px")

enrich_all <- enrich_bar+enrich_point+enrich_pie+plot_layout(ncol=3,nrow=1,widths = c(1,1,1.4))

#--------------------------------------------------------
#pic_pie02 <- paste0("analysis results/Metabolites Enrichment_pie02 (",libraries,").png")
#png(file=pic_pie02, width=1000, height =800)

pie(enrich_df_pie$Hits, labels=enrich_df_pie$term,init.angle = 90,
            edges=200,radius=1,clockwise = T,border = "white",col = my_col,cex=0.8,
            main = "The main enriched metabolite sets")


#-------------pathways analsysis-----------------------------------------------

spe <- c("human","mouse","rat")
species <- tolower(species)

if(species %in% spe == FALSE)
  stop ("The species is not found, and the pathways analysis can not be performed! Note: The species should be human, mouse, or rat.") else
  
  {if(species=="human")
    {m_lib<-SetKEGG.PathLib(m_obj, "hsa", "current")
    web="https://rest.kegg.jp/list/pathway/hsa"
    path_file <- c("analysis results/human_metabolic_pathways.xlsx")
    path_pic <- c("analysis results/human_metabolic_pathways.png")
    meta_enrichment <- c("analysis results/human_metabolites_enrichment_analysis.png")
    t <- url(web,encoding ="UTF-8")  # library(XML)
    d<-scan(t, what=character()) %>% paste0(collapse = " ")
    d_c <- gsub("- Homo sapiens (human)",";",d,fixed = TRUE)} else
  
  {if(species=="mouse")
    {m_lib<-SetKEGG.PathLib(m_obj, "mmu", "current")
    web="https://rest.kegg.jp/list/pathway/mmu"
    path_file <- c("analysis results/mouse_metabolic_pathways.xlsx")
    path_pic <- c("analysis results/mouse_metabolic_pathways.png")
    meta_enrichment <- c("analysis results/mouse_metabolites_enrichment_analysis.png")
    t <- url(web,encoding ="UTF-8")  # library(XML)
    d<-scan(t, what=character()) %>% paste0(collapse = " ")
    d_c <- gsub("- Mus musculus (house mouse)",";",d,fixed = TRUE)} else
   
   {m_lib<-SetKEGG.PathLib(m_obj, "rno", "current")
    web="https://rest.kegg.jp/list/pathway/rno"
    path_file <- c("analysis results/rat_metabolic_pathways.xlsx")
    path_pic <- c("analysis results/rat_metabolic_pathways.png")
    meta_enrichment <- c("analysis results/rat_metabolites_enrichment_analysis.png")
    t <- url(web,encoding ="UTF-8")  # library(XML)
    d<-scan(t, what=character()) %>% paste0(collapse = " ")
    d_c <- gsub("- Rattus norvegicus (rat)",";",d,fixed = TRUE)
    }
    }
  }


#------------------------------------------------------------------------------
 
  m_filt<-SetMetabolomeFilter(m_lib, F);
  
  m_scor<-CalculateOraScore(m_filt, "rbc", "hyperg")

  path_result <- read.csv("pathway_results.csv",header=T)

  file.remove("pathway_results.csv")

  colnames(path_result) <- c("path_id","Total","Expected","Hits","p value","-log10(p)","Holm p","FDR","Impact")


    dw <- strsplit(d_c, "[;]") %>% data.frame()

    dw[,1] <- trimws(dw[,1])

    dw$path_id <- substring(dw[,1],6,14)

    dw$path_name <- substring(dw[,1], 15, str_length(dw[,1]))  %>% trimws()

    hsa_paths <- dw[-1]

    write.xlsx(hsa_paths,path_file)


  path_ls <- read.xlsx(path_file,colNames=T) %>% data.frame()

  path_ls$path_id <- trimws(as.character(path_ls$path_id))

  path_result$path_id <- trimws(as.character(path_result$path_id))

  path_df <- left_join(path_result,path_ls,by="path_id")

  head(path_df)

  path_data <-cbind(path_df[1],path_df[10],path_df[,2:9])
  colnames(path_data) <- c("path_id","path_name","Totla","Expected","Hits","Pvalue","-log10(p)","Holm p","FDR","Impact")

  write.xlsx(path_data,path_file,rownames=T,colnames=T)


  title_size <- case_when(nrow(path_data)>30 ~13,
                          nrow(path_data)>20 ~12,
                          TRUE ~11)

  xy_size <- case_when(nrow(path_data)>30 ~8,
                       nrow(path_data)>20 ~9,
                       TRUE ~10)

  my_theme<-theme_bw()+
    theme(text=element_text(family = "sans",colour ="black",face="bold",size =title_size),
          panel.grid = element_blank(),
          panel.border = element_rect(linewidth = 0.8,colour = "gray30"),
          axis.line = element_blank(),
          axis.ticks = element_line(linewidth = 0.6,colour = "gray30"),
          axis.ticks.length = unit(1.5,units = "mm"))+
    theme(plot.title = element_text(hjust = 0.5))


  xy_theme <-theme(axis.text.x = element_text(face="bold",color="black",size=xy_size,angle =0,hjust=1))+
    theme(axis.text.y = element_text(face="bold",color="black",size=xy_size))+
    theme(legend.text=element_text(face="bold",color="black",size=10))


  x_ran <- max(path_data$Impact)*1.05
  y_ran <- nrow(path_data)

  meta_path <- ggplot(path_data,
                      aes(x=Impact,
                          y=fct_reorder(path_name,-log10(Pvalue))))+
    geom_point(aes(color=-log10(Pvalue),size=Impact))+
    scale_color_gradient2(midpoint = 0,low = "#33ffff",mid = "#ffff00",high = "#ff0000")+
    labs(x = 'Pathway Impact', y = '',title="Overview of Metabolism Pathway Anlysis")+
    my_theme+xy_theme+
    theme(panel.grid =element_line(colour="#dcdcdc",linewidth=0.2,linetype = "dashed"))+
    coord_cartesian(xlim=c(0, x_ran),expand = T)+
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

  meta_path

  ggsave(path_pic, meta_path, width=1200, height =1000, dpi=150,units = "px")

  #file.remove("compound_db.qs","smpdb_pathway.qs","tosend.rds","syn_nms.qs")

  print("---------------------------------------------------------------")
  print("The analysis results can be found in the folder of <analysis result>")

  allfiles <- dir()
  qs <- grep("*.qs",allfiles)
  rds <- grep("*.rds",allfiles)
  file.remove(allfiles[qs], allfiles[rds])


  pp <- enrich_point+enrich_pie+meta_path+plot_layout(ncol=3,nrow=1,widths = c(1,1.2,1))

  ggsave(meta_enrichment, pp, width=3400, height =1000, dpi=150,units = "px")

  dev.new(width=18,height=7, noRStudioGD = T)

  pp

  }

