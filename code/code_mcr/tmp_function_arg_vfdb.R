tmp_function_arg<-function(path="H:/SJYT_metagenomics/unique_gene/resistance/",
                           path_silva="H:/SJYT_metagenomics/unique_gene/silva/",
                           no_normalize=F,group_meta=NULL,sample_name=NULL, self_index = F,
                           path_db="G:\\stone_meta\\database\\resistome\\ARG_CARD/ARG_SARG/Long_subdatabase/"){
  library(tidyverse)
  library(phyloseq)
  tmp<-data.table::fread(paste0(path,"abundance_table_reads.tsv"),data.table=F)
  rownames(tmp)<-data.table::fread(paste0(path,"gene_name.txt"),header = F,sep=NULL,data.table=F)[,1] ; tmp_count<-tmp#;colnames(tmp_count)<-sample_name##ARG
  tmp2<-data.table::fread(paste0(path,"gene_length.txt"),data.table=F)[,1]# ;colnames(tmp)<-sample_name;colSums(tmp)
  mat<-tmp/tmp2
  tmp<-data.table::fread(paste0(path_silva,"abundance_table_reads.tsv"),data.table=F) ##16S rRNA gene
  tmp2<-data.table::fread(paste0(path_silva,"gene_length.txt"),data.table=F)[,1]# ;colnames(tmp)<-sample_name;colSums(tmp)
  # tmp2<-data.table::fread("Y:/database/silva_index/silva_index_length.txt",data.table=F)[,2]# ;colnames(tmp)<-sample_name;colSums(tmp)
  copies<-colSums(tmp/tmp2)
  for(i in 1:ncol(mat)){mat[,i]<-mat[,i]/copies[i]}
  print(paste0("Removed samples without silva 16S rRNA gene: ",names(which(is.na(colSums(mat))))))
  print(paste0("Removed samples without ARG: ",names(which(colSums(mat)==0))))

  if(no_normalize) mat<-tmp_count
  mat<-mat[,colSums(mat)!=0&!is.na(colSums(mat))] 
  # return(which(is.infinite(colSums(mat))))

  if(self_index){
    df<-read.blast(paste0(path,"unique_gene.faa.diamond.blast")) 
    df<-cbind(df,mat)
    df<-df %>% filter(pident>90,evalue<1e-10,length>=25) %>%  select(-3:-12);dim(df) # %>% melt()
    rownames(df)<- df$qseqid
    
    db_arg<-read.delim(paste0(path_db,"3.SARG_v3.2_20220917_Long_subdatabase_structure.txt"))[,-1]
    rownames(db_arg)<-db_arg$SARG.Seq.ID
    db_arg$Type[db_arg$Type=="macrolide-lincosamide-streptogramin"]<-"MLS"
    db_arg<-db_arg %>% mutate(Name=str_sub(Subtype,str_locate_my_first(Subtype,"__")+2)) 
    tax<-db_arg[df$sseqid,];rownames(tax)<-df$qseqid
    
    df<-df %>% select(-1:-2)
    db_arg<-db_arg %>% select(Type:Mechanism.subgroup2)

  } else {
    db_arg<-read.delim(paste0(path_db,"Nucleotide_SARG_v3.2_Long_structure.txt"))
    db_arg$Type[db_arg$Type=="macrolide-lincosamide-streptogramin"]<-"MLS"
    db_arg<-db_arg %>% mutate(Name=str_sub(Subtype,str_locate_my_first(Subtype,"__")+2)) 
    rownames(db_arg)<-db_arg$NucleotideID
    db_arg<-db_arg %>% select(Type:Mechanism.subgroup2,Name,NucleotideID)
    df <-mat; tax <- db_arg
  }
  
  
  
  if(!is.null(group_meta)){
    ps_arg<-phyloseq(otu_table(df,T),
                     sample_data(group_meta),
                     tax_table(as.matrix(tax))) %>% prune_taxa(taxa_sums(.)>0,.)
  } else {
    ps_arg<-phyloseq(otu_table(df,T),
                     tax_table(as.matrix(tax))) %>% prune_taxa(taxa_sums(.)>0,.)
  }
 
  
  return(ps_arg)
}

# path="H:/analysis/sjyt_yh/unique_gene/resistance/"
# path_silva="H:/SJYT_metagenomics/unique_gene/silva/"
# no_normalize=F
# group_meta=df
# sample_name=df$sample_name


##20240215
tmp_function_vfdb<-function(path="H:/SJYT_metagenomics/unique_gene/vfdb/",
                            path_silva="H:/SJYT_metagenomics/unique_gene/silva/",
                            no_normalize=F,group_meta=NULL,sample_name=NULL,self_index=F){
  tmp<-data.table::fread(paste0(path,"abundance_table_reads.tsv"),data.table=F) ##ARG
  
  if(no_normalize) {
    mat<-tmp#;colnames(mat)<-sample_name
  } else {
    tmp2<-data.table::fread(paste0(path,"gene_length.txt"),data.table=F)[,1]# ;colnames(tmp)<-sample_name;colSums(tmp)
    mat<-tmp/tmp2
    tmp<-data.table::fread(paste0(path_silva,"abundance_table_reads.tsv"),data.table=F)##16S rRNA gene
    tmp2<-data.table::fread(paste0(path_silva,"gene_length.txt"),data.table=F)[,1]# ;colnames(tmp)<-sample_name;colSums(tmp)
    copies<-colSums(tmp/tmp2)
    for(i in 1:ncol(mat)){mat[,i]<-mat[,i]/copies[i]}
    #colnames(mat)<-sample_name
  }
  
  rownames(mat)<-data.table::fread(paste0(path,"gene_name.txt"),data.table=F,header = F)[,1]
  
  ##database
  db<-fread("G:\\stone_meta\\database\\resistome\\VF/VFDB_setB_pro.txt",data.table = F,header = F,sep = ";") %>% 
    mutate(info=str_sub(V1,str_locate_my_first(V1,"\\[")+1,str_locate_my_first(V1,"\\]")-1)) %>% 
    mutate(VFC=str_sub(info,str_locate_my_first(info,"-")+2),
           VF=str_sub(info,1,str_locate_my_first(info,"\\)")),
           name=str_sub(V1,1,str_locate_my_first(V1,"\\(")-1)) %>% 
    column_to_rownames("name") %>% select(-V1,-info) %>% 
    mutate(VFID=str_sub(VF,str_locate_my_first(VF,"[(]")+1,str_length(VF)-1)) ;db$entry=rownames(db)
  tmp<-read.xlsx("G:\\stone_meta\\database\\resistome\\VF/VFs.xlsx",1) %>% select(VFID,Type=VFcategory)
  db<-db %>% left_join(tmp,by=c("VFID"));rownames(db)<-db$entry
  
  
  if(self_index){
    df<-read.blast(paste0(path,"unique_gene.faa.diamond.blast")) %>% mutate(sseqid=str_sub(sseqid,1,str_locate_my_first(sseqid,"\\(")-1))
    df<-cbind(df,mat)
    df<-df %>% filter(pident>90,evalue<1e-10,length>=25) %>%  select(-3:-12);dim(df) # %>% melt()
    
    tax<-db[df$sseqid,];rownames(tax)<-df$qseqid
    df<-df %>% select(-1,-2)
  } else {
    df <-mat; tax <- db
  }
  
  
  
  
  if(!is.null(group_meta)){
    ps_vfdb<-phyloseq(otu_table(df,T),
                      sample_data(group_meta),
                      tax_table(as.matrix(tax))) %>% prune_taxa(taxa_sums(.)>0,.)
  } else {
    ps_vfdb<-phyloseq(otu_table(df,T),
                      tax_table(as.matrix(tax))) %>% prune_taxa(taxa_sums(.)>0,.)
  }
  
  return(ps_vfdb)
}

##
# tmp_function_vfdb<-function(path="H:/SJYT_metagenomics/unique_gene/vfdb/",
#          path_silva="H:/SJYT_metagenomics/unique_gene/silva/",
#          no_normalize=T,group_meta=group_meta,sample_name=sample_name){
#   tmp<-data.table::fread(paste0(path,"abundance_table_reads.tsv"),data.table=F);tmp_count<-tmp;colnames(tmp_count)<-sample_name ##ARG
#   tmp2<-data.table::fread(paste0(path,"gene_length.tsv"),data.table=F)[,1]# ;colnames(tmp)<-sample_name;colSums(tmp)
#   mat<-tmp/tmp2
#   tmp<-data.table::fread(paste0(path_silva,"abundance_table_reads.tsv"),data.table=F)##16S rRNA gene
#   tmp2<-data.table::fread("Y:/database/silva_index/length.txt",data.table=F)[,2]# ;colnames(tmp)<-sample_name;colSums(tmp)
#   copies<-colSums(tmp/tmp2)
#   for(i in 1:ncol(mat)){mat[,i]<-mat[,i]/copies[i]}
#   colnames(mat)<-sample_name
#   
#   if(no_normalize) mat<-tmp_count
#   
#   df<-read.blast(paste0(path,"unique_gene.faa.diamond.blast")) %>% mutate(sseqid=str_sub(sseqid,1,str_locate_my_first(sseqid,"\\(")-1))
#   df<-cbind(df,mat)
#   df<-df %>% filter(pident>90,evalue<1e-10,length>=25) %>%  select(-3:-12);dim(df) # %>% melt()
#   db<-fread("G:\\stone_meta\\database\\resistome\\VF/VFDB_setB_pro.txt",data.table = F,header = F,sep = ";") %>% 
#     mutate(info=str_sub(V1,str_locate_my_first(V1,"\\[")+1,str_locate_my_first(V1,"\\]")-1)) %>% 
#     mutate(VFC=str_sub(info,str_locate_my_first(info,"-")+2),
#            VF=str_sub(info,1,str_locate_my_first(info,"\\)")),
#            name=str_sub(V1,1,str_locate_my_first(V1,"\\(")-1)) %>% 
#     column_to_rownames("name") %>% select(-V1,-info);db$entry=rownames(db)
#   tmp<-read.xlsx("G:\\stone_meta\\database\\resistome\\VF/VFs.xlsx",1) %>% select(VF=VFID,Type=VFcategory)
#   db %>% left_join(tmp,by=c("VF"))
#   tax<-db[df$sseqid,];rownames(tax)<-df$qseqid
#   df<-df %>% column_to_rownames("qseqid") %>% select(-sseqid)
#   ps_vfdb<-phyloseq(otu_table(df,T),
#                     sample_data(group_meta),
#                     tax_table(as.matrix(tax))) %>% prune_taxa(taxa_sums(.)>0,.)
#   return(ps_vfdb)
# }
