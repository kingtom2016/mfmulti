wd="H:/SJYT_metagenomics/unique_MAG/"
df_mag_contig_snp<-fread(paste0(wd,"MetaPop/all_sample_01_snv/07.Cleaned_SNPs/genic_snps.tsv"), data.table=F) %>%
  mutate(contig_bin=contig,
         # contig_name=str_sub(contig,1,str_locate_my_first(contig,"_SGB_")-1),
         MAG=str_sub(contig,str_locate_my_first(contig,"SGB_"))) %>% 
  select(-contig,-contig_gene,-contig_bin) %>% mutate(sample_name=str_remove(source,".sort")) %>% mutate(Habitat=short_name[sample_name]) %>% select(-source)
# df_mag_contig_snp<-fread(paste0(wd,"MetaPop/07.Cleaned_SNPs/non_genic_snps.tsv"), data.table=F)
# df_mag_contig<-fread(paste0(wd,"stat/bins_contigs.tsv"),data.table = F,header = F) %>% select(contig=V1,MAG=V2) %>% mutate(MAG=str_sub(MAG,1,-4)) 


df_tmp <- rbind(df_mag_contig_snp %>% filter(str_length(snps)==1),
                df_mag_contig_snp %>% filter(str_length(snps)>1) %>% separate_rows(snps,sep = '(?<=.)(?=.)')) ;gc()
df_tmp<-df_tmp %>% mutate(snps_num=ifelse(snps=="A",a_ct,
                                          ifelse(snps=="T",t_ct,
                                                 ifelse(snps=="C",c_ct,
                                                        ifelse(snps=="G",g_ct,NA)))))
otu_snps<- df_tmp %>% mutate(contig_pos_snps=paste0(contig_pos,"_",snps)) %>% select(contig_pos_snps,Habitat,snps_num) %>% group_by(contig_pos_snps, Habitat) %>% summarise(snps_num = sum(snps_num)) 
otu_snps<- otu_snps %>% pivot_wider(names_from =Habitat,values_from = snps_num,values_fill = 0) ;gc()


otu_snps_info<-otu_snps %>% ungroup() %>% #head() %>% 
  mutate(contig_name=str_sub(contig_pos_snps,1,str_locate_my_first(contig_pos_snps,"_",F)-1),
         contig_name=str_sub(contig_name,1,str_locate_my_first(contig_name,"_",F)-1)) %>% 
  mutate(MAG=str_sub(contig_name,str_locate_my_first(contig_name,"SGB_"))) %>% select(-1)


list_tmp<-list()
for(i in c( "Water" ,   "Sediment", "Soil"   ,  "Leaf" ,    "Fish"   )){
  otu<-otu_snps_info %>% filter(MAG %in% list_incidant[[i]])
  df<- otu %>% select(-contig_name) %>% group_by(MAG) %>% summarise_all(sum) %>% data.frame()
  df<- df %>% mutate(snp_called=!!sym(i)>0) %>% mutate(snp_called_also_other=rowSums(binary_otu(df %>% select(-!!sym(i),-MAG),F))>0)

  print("SNP-called mags number");print(sum(df$snp_called)/nrow(df));print(sum(df$snp_called_also_other)/nrow(df))
  if(T) otu<-otu %>% filter(MAG%in%df[df$snp_called==T&df$snp_called_also_other==T,]$MAG)

  otu<-otu[,1:5] 

  tmp_res<-FEAST_Paral_two_table(otu %>% select(-!!sym(i)),otu %>% select(!!sym(i)),COVERAGE = 1e8)
  list_tmp[[i]] <-tmp_res %>% rownames_to_column("source") %>% melt(variable.name = "sink")

}
df_snp_source<-bind_rows(list_tmp)

df_snp_source %>% mutate(source=fct_relevel(source,names(col_habitat))) %>% 
  mutate(source=ifelse(source=="Unknown","Unknown source","Known source from other habtiats"),source=as.character(source)) %>%
  ggplot(aes(x=sink,y=value))+
  geom_col(aes(fill=source))+
  scale_fill_manual(values=c(col_habitat,"Unknown"="gray","#7f7f7f"))+
  scale_fill_manual(values=c("#7f7f7f","gray"))+
  labs(y="Relative contribution",fill="")
ggsave("p_snp_source.pdf",width=7,height=4);rstudioapi::viewer("p_snp_source.pdf")
