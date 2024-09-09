
FEAST_Paral_two_table<-function(otu_source,otu_target, Env=NA , cores=4 ,COVERAGE=10000,EM_iterations = 500){
  ##otu_source and otu_source are dataframe with sample as column and OTU/ASV as rows
  if (!requireNamespace("FEAST", quietly = TRUE)) {
    devtools::install_github("cozygene/FEAST")
  }

  library(dplyr)
  library(foreach)
  library(doParallel)
  library(FEAST)

  ###combine table
  otu_all<-merge(otu_target,otu_source, by="row.names", all=T); rownames(otu_all)<-otu_all[,1];otu_all<-otu_all[,-1] ; otu_all[is.na(otu_all)]<-0

  rownames(otu_all)<-1:dim(otu_all)[1]  ##
  name_source<-colnames(otu_source)
  name_target<-colnames(otu_target)

  ######
  metadata_feast<-matrix(,nrow = length((name_source))+1,ncol=3) %>% as.data.frame()
  rownames(metadata_feast)<-c("target",name_source)
  colnames(metadata_feast)<-c("Env","SourceSink","id")
  metadata_feast[,3]<-1
  metadata_feast[1,2]<-"Sink"
  metadata_feast[-1,2]<-"Source"
  metadata_my<-metadata_feast

	print("Start calculation!")
  ##
  cl <- makeCluster(cores)
  registerDoParallel(cl, cores=cores)
  FEAST_output <- data.frame(stringsAsFactors = F)


  #@######################################
  FEAST_output <-foreach(i=1:length(name_target), .combine='cbind',.verbose = F)  %dopar% {
    library(FEAST)
    rownames(metadata_my)[1]<-name_target[i]

    FEAST_output1 <- FEAST(
      C =t(otu_all),
      metadata = metadata_my,
      different_sources_flag = 1,dir_path=getwd(),
      outfile = "demo", COVERAGE=COVERAGE,EM_iterations = EM_iterations
    )
    print(FEAST_output1)

    return(FEAST_output1$data_prop[,1])
  }
  if(length(name_target)==1) FEAST_output<-as.data.frame(FEAST_output)
  stopImplicitCluster()
  stopCluster(cl)


  colnames(FEAST_output)<-name_target
  rownames(FEAST_output)<-c(colnames(otu_source), "Unknown")
  print("Analysis complete!")
  print("Please cite: Shenhav L, Thompson M, Joseph TA, Briscoe L, Furman O, Bogumil D, et al., 2019. FEAST: fast expectation-maximization for microbial source tracking. Nat. Methods 16, 627-632. https://doi.org/10.1038/s41592-019-0431-x.")
  print("He J, Zhang N, Muhammad A, Shen X, Sun C, Li Q, et al., 2022. From surviving to thriving, the assembly processes of microbial communities in stone biodeterioration: A case study of the West Lake UNESCO World Heritage area in China. Sci. Total Environ. 805, 150395. https://doi.org/10.1016/j.scitotenv.2021.150395.")
  print("He J, Zhang N, Shen X, Muhammad A, Shao Y, 2022. Deciphering environmental resistome and mobilome risks on the stone monument: A reservoir of antimicrobial resistance genes. Sci. Total Environ. 838, 156443. https://doi.org/10.1016/j.scitotenv.2022.156443.")
  return(as.data.frame(FEAST_output))

}
