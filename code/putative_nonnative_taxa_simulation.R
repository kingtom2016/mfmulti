library(tidyverse)
library(furrr)
plan(multisession, workers = 10)
for(pri in c("Bacteria","Fungi")){
  ##choose samples
  if(pri=="Bacteria") { 
    ps<-ps1
  }else{
    ps<-ps2
  }
  ##Nich width
  df_niche_breadth<-foreach(i=group_name_all,.combine = cbind) %do% {
  otu<-ps %>% subset_samples(Habitat==i) %>% otu_table()
  tmp<-t(spaa::niche.width(as.data.frame(t(otu)), "levins"))
  };colnames(df_niche_breadth)<-group_name_all;df_niche_breadth[is.infinite(df_niche_breadth)]<-0;df_niche_breadth<-as.data.frame(df_niche_breadth)
  df_niche_breadth<- apply(df_niche_breadth,2,function(x){x/max(x)}) %>% as.data.frame()
  df_niche_breadth$high_niche_habitat<-apply(df_niche_breadth, 1, function(x) {return( names(x)[which.max(x)])})
  ##
  for(i in group_name_all) { 
    otu<-data.frame(otu_table(ps %>% transform("compositional") %>% merge_samples2("Habitat",mean)))
    df<-otu
    df$occurring_habitat_number=rowSums(binary_otu(otu,F))
    df$occurring_this_habitat<-ifelse(otu[,i]>0,1,0)
    df$enriched_habitat<-apply(otu, 1, function(x) {return( names(x)[which.max(x)])})
    
    df$high_niche_habitat<-df_niche_breadth[rownames(df),]$high_niche_habitat
    df$potential_transient<-ifelse(df$enriched_habitat!=i&df$occurring_this_habitat==1&df$high_niche_habitat!=i, "potential transient","not transient")
    tmp_selected=rownames(df[df$potential_transient=="potential transient",])
    selected<-tmp_selected
    
      otu<-ps %>% subset_samples(Habitat==i) %>% prune_taxa_nonzero() %>% otu_table()
      tmp<-t(spaa::niche.width(as.data.frame(t(otu)), "levins"))
      tmp1<-tmp[tmp_selected,1]
      tmp2<-tmp
      df<-future_map_dfr(names(tmp1), ~{
      test_result <- wilcox.test(x = tmp1, mu = tmp1[.x]) 
      tibble(
        name = .x,
        vec1_value = tmp1[.x],
        statistic = test_result$statistic,
        p.value = test_result$p.value)
    }, .progress = TRUE)
      selected<-df %>% mutate(p.adj=p.adjust(p.value,"BH")) %>% filter(p.adj<0.001,vec1_value<mean(tmp1)) %>% pull("name")
    
    list_tmp[[i]]<-selected
    
  }
  list_incidant[[pri]]<-list_tmp
}
plan(sequential)

binary_otu <- function (data, physeq_data = T, boundary = 0) 
{
    if (physeq_data) {
        otu_table(data)[otu_table(data) <= boundary] <- 0
        otu_table(data)[otu_table(data) > boundary] <- 1
    }
    else {
        data[data <= boundary] <- 0
        data[data > boundary] <- 1
    }
    return(data)
}

alpha_diversity <- function (physeq, apped_group = T, phy_div = F) 
{
    df = cbind(Shannon = microbiome::diversity(physeq, "Shannon")[, 
        1], Richness = sample_sums(binary_otu(physeq)))
    if (phy_div) {
        df = cbind(df, PD_whole_tree = pd(t(data.frame(otu_table(physeq))), 
            phy_tree(physeq), include.root = FALSE)[1])
    }
    if (apped_group) {
        df = cbind(df, microbiome::meta(physeq))
    }
    return(as.data.frame(df))
}

res_incidant <- foreach(pri=c("Bacteria","Fungi"),.combine = rbind) %:% foreach(i=group_name_all,.combine = rbind) %do% {
  if(pri=="Bacteria"){
    selected<-list_incidant[["Bacteria"]][[i]]
    ps<-ps1
  } else {
    selected<-list_incidant[["Fungi"]][[i]]
    ps<-ps2
  }
  ps<-ps %>% subset_samples(Habitat==i)
  ps_tmp<-ps %>% prune_taxa(!taxa_names(.)%in%selected,.)
  ##calculate richness of observed and simulated community
  df_res<-rbind(
    data.frame(alpha_diversity(ps),type="observed",sample_sums=sample_sums(ps)),
    data.frame(alpha_diversity(ps_tmp),type="simulated",sample_sums=sample_sums(ps_tmp)))#df_res
  return(df_res %>% mutate(Habitat=i) %>% mutate(primer=pri))
}