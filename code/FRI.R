##otu is dataframe column as sample row as ASV/OTU, phy is a tree object.
##ref_kegg is dataframe with column of KO, PathwayL1 PathwayL2 PathwayL3 KoDescription  ##Here filter only those element-related KO
# res_fri<-calculateFunctionalRedundancy_parallel_metagenome(otu,phy,cores=12,
                                                           # ref_kegg =ref_kegg %>% filter(KO%in%tmp_ref$KO), 
                                                           # path_to_temp_folder="tmp1")
														   
calculateFunctionalRedundancy_parallel_metagenome<-function (otu,phy,ref_kegg_file ,path_to_temp_folder="tmp",  cores=12)
{
  # dir.create(path_to_temp_folder)
  otu[otu>0]<-1
  otu_table_reduced_aggregated=otu %>% as.data.frame %>% rownames_to_column()
  distance_matrix = cophenetic(phy)
  distance_matrix_reduced = distance_matrix[rownames(otu),rownames(otu)]
  distance_matrix_mean = mean(as.dist(distance_matrix))
  distance_matrix_reduced_mean = mean(as.dist(distance_matrix_reduced))
  rm(distance_matrix)

  reference_profile = as.data.frame(t(otu_table(ps_mag_kegg)))
  dim(reference_profile)
  ko_list = ref_kegg_file %>% select(ko=KO, description=KoDescription) %>% distinct(ko,.keep_all = T)
  reference_profile<-reference_profile[,colnames(reference_profile)%in%ko_list$ko]
  ko_list = ko_list[match(colnames(reference_profile),ko_list$ko),]


  abs_functional_redundancy_list = list()
  rel_functional_redundancy_list = list()
  library(doParallel)
  library(foreach)
  registerDoParallel(cores = cores)
  results = foreach(sample = 2:ncol(otu_table_reduced_aggregated), .combine = 'c', .packages = c("ape")) %dopar% {
    print(paste("Calculate functional redundancy for sample", colnames(otu_table_reduced_aggregated)[sample]))
    write.table(sample, file = "tmp/tmp.tsv",append = T)
    functional_prediction_sample = reference_profile * as.numeric(otu_table_reduced_aggregated[, sample])
    functional_prediction_sample_mod = ifelse(functional_prediction_sample >= 1, 1, 0)

    abs_functional_redundancy_sample = vector("list", length = nrow(ko_list))
    rel_functional_redundancy_sample = vector("list", length = nrow(ko_list))


    for (i in 1:nrow(ko_list)) {
      ko_count = functional_prediction_sample_mod[, i]
      sum(ko_count)
      FRI<-(mean(as.dist(distance_matrix_reduced * ko_count)) * (sum(ko_count)/length(ko_count)))
      abs_functional_redundancy_sample[[i]] = FRI/distance_matrix_mean
      rel_functional_redundancy_sample[[i]] = FRI/distance_matrix_reduced_mean
      print(i)
    }


    list(abs_functional_redundancy_sample = do.call(cbind, abs_functional_redundancy_sample),
         rel_functional_redundancy_sample = do.call(cbind, rel_functional_redundancy_sample))
  }
  stopImplicitCluster()
  # Post-processing to format results correctly
  # abs_functional_redundancy_tab = do.call(cbind, lapply(results, `[[`, "abs_functional_redundancy_sample"))
  # rel_functional_redundancy_tab = do.call(cbind, lapply(results, `[[`, "rel_functional_redundancy_sample"))

  abs_functional_redundancy_tab = do.call(rbind, results[2*(1:length(names(results)))-1]) %>% t()
  rel_functional_redundancy_tab = do.call(rbind, results[2*(1:length(names(results)))])%>% t()

  abs_functional_redundancy_tab = data.frame(abs_functional_redundancy_tab)
  rel_functional_redundancy_tab = data.frame(rel_functional_redundancy_tab)

  colnames(abs_functional_redundancy_tab) = names(otu_table_reduced_aggregated)[2:ncol(otu_table_reduced_aggregated)]
  colnames(rel_functional_redundancy_tab) = names(otu_table_reduced_aggregated)[2:ncol(otu_table_reduced_aggregated)]

  # Stop the parallel backend

  abs_functional_redundancy_final = data.frame(KO = ko_list$ko, abs_functional_redundancy_tab)
  rel_functional_redundancy_final = data.frame(KO = ko_list$ko,  rel_functional_redundancy_tab)

  # write.table(x = abs_functional_redundancy_final, file = file.path(path_to_temp_folder,
  #                                                                   "absolute_functional_redundancy.txt"), append = F, quote = F,  sep = "\t", row.names = F, col.names = T)
  # write.table(x = rel_functional_redundancy_final, file = file.path(path_to_temp_folder,
  #                                                                   "relative_functional_redundancy.txt"), append = F, quote = F,  sep = "\t", row.names = F, col.names = T)
  return(list(abs_functional_redundancy_final,rel_functional_redundancy_final))
  }

##Script is modified from function of Tax4Fun2 package

