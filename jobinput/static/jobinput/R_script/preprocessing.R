
#' read the data
#' @param file_name
#'
read_file<-function(file_name){
  library(data.table)
  df_use<-fread(file_name)
  return(df_use)

}


#' calculate the cpm
#' @param gene_counts count df with rows as genes and columns as samples, the first column should be genes,colnames is samples ID

calculate_cpm<-function(gene_counts){
  library(edgeR)
  gene_counts<-gene_counts[order(gene_counts[,1]),]
  expression_count<-DGEList(as.matrix(gene_counts[,2:ncol(gene_counts)]),genes = gene_counts[,1])
  isexpr <- rowSums(cpm(expression_count)>=1) >=2
  expression_count<- expression_count[isexpr,]
  expression_count$samples$lib.size <- colSums(expression_count$counts)
  expression_count <- calcNormFactors(expression_count, method="upperquartile")
  final_cpm<-cpm(expression_count, normalized.lib.sizes=TRUE, log=T)
  row.names(final_cpm)<-expression_count$genes$genes
  colnames(final_cpm)<-colnames(gene_counts)[2:ncol(gene_counts)]
  return(final_cpm)


}


#' calculate the raw count
#' @param cpm_matrix a number matrix with the cpm values.
change_from_cpm_to_raw_count<-function(cpm_matrix){
  #we assume library size for each sample is 2000000
  not_log_cpm<-2 ^ as.matrix(cpm_matrix)
  pseudocount<-not_log_cpm*2000000/1000000
  mode(pseudocount) <- "integer"
  pseudocount[is.na(pseudocount)]<-0
  return(as.data.frame(pseudocount))
}



#' transfer IDs to gene symbol
#' @param df1 character as IDs need to be transfer
#'
transfer_IDs<-function(df1){
  library(org.Hs.eg.db)
  genes_ID<-df1
  if (grepl('ENSG*.',genes_ID[1])) {
    #delete the gene ID with version information. Like ENSG00000186092.4 remove the .4
    genes_ID<-gsub("\\.[0-9]+","",genes_ID)
    gene_id_return<-mapIds(org.Hs.eg.db,keys=genes_ID,keytype='ENSEMBL',column="SYMBOL",multiVals = 'first')

  }
  else if (!grepl("\\D",genes_ID[1])){
    gene_id_return<-mapIds(org.Hs.eg.db,keys=genes_ID,keytype='ENTREZID',column="SYMBOL",multiVals = 'first')

  }
  else {
    gene_id_return<-genes_ID
  }
  return(gene_id_return)

}

#' update gene symbol used
#' @param df  df with rows as genes and columns as samples, colnames is samples ID
#' @param column_num which column has gene symbol
update_gene_symbol<-function(df,column_num){
  library(HGNChelper)
  library(dplyr)
  colnames(df)[column_num]<-'symbol'
  mapping_df<-checkGeneSymbols(df[,column_num])
  df_new<-merge(mapping_df,df,by.x='x',by.y='symbol')
  #remove genes without new symbol
  df_new<-df_new[!is.na(df_new$Suggested.Symbol),]
  #selected the first one in Suggested symbol
  df_new$Suggested.Symbol<-unlist(lapply(strsplit(df_new$Suggested.Symbol,split=" /// "),function(x){return(x[1])}))
  df_new<-df_new[,3:ncol(df_new)]
  colnames(df_new)[1]<-'symbol'
  df_new<-df_new %>% group_by(symbol) %>% summarise_all(sum)
  return(df_new)
}




#' batch effect removal
#' @param count_df count df with rows as genes and columns as samples, the first column should be genes,colnames is samples ID
#'
batch_effect_removal<-function(count_df,dir_store_files,methods='combat',batch_effect=NULL){
  if (methods=='ruvg'){
  #ruvg
  library(RUVSeq)
  IMU_KRT_not_sig_genes<-read.csv(paste0(dir_store_files,'/IMU_KRT_not_sig_genes.csv'),row.names = 1)
  row.names(count_df)<-count_df[,1]
  set <- newSeqExpressionSet(as.matrix(count_df[,2:ncol(count_df)]))
  idx  <- rowSums(counts(set) > 5) >= 2
  set  <- set[idx, ]
  temp_genes<-IMU_KRT_not_sig_genes[,1][IMU_KRT_not_sig_genes[,1] %in% row.names(set)]
  set2<-RUVg(set,temp_genes,k=1)
  Ruvg_imporved<-as.data.frame(set2@assayData$normalizedCounts)
  Ruvg_imporved<-cbind(row.names(Ruvg_imporved),Ruvg_imporved)
  colnames(Ruvg_imporved)[1]<-'symbol'
  return(Ruvg_imporved)
  }
  else if (methods=='combat'){
    #combat
    library(sva)
    colnames(count_df)[1]<-'symbol'
    if (is.null(batch_effect)){
    FF_df<-read.csv(paste0(dir_store_files,'/FF_df_symbol_corrected.csv'),row.names = NULL,check.names = F)
    batch_effect<-c(rep('other',ncol(count_df)-1),rep('FF',18))
    final_df<-merge(count_df,FF_df,by='symbol')
    adjusted_combat <- ComBat_seq(as.matrix(final_df[,2:ncol(final_df)]), batch=batch_effect, group=NULL)
    adjusted_combat<-as.data.frame(adjusted_combat)
    adjusted_combat<-cbind(final_df$symbol,adjusted_combat)
    colnames(adjusted_combat)[1]<-'symbol'
}
    else{
      batch_effect=batch_effect
      final_df<-count_df
      adjusted_combat <- ComBat_seq(as.matrix(final_df[,2:ncol(final_df)]), batch=batch_effect, group=NULL)
      adjusted_combat<-as.data.frame(adjusted_combat)
      adjusted_combat<-cbind(final_df$symbol,adjusted_combat)
      colnames(adjusted_combat)[1]<-'symbol'
    }

    return(adjusted_combat)
  }
}

#' examine genes available in gene set
#' @param count_df count cpm
#' @param geneset_list list store all gene sets and
#'

whether_perform_PCA_for_next_step<-function(count_df,json_file){
  library(rjson)
  gene_set_list <- fromJSON(file=json_file)
  all_gene_set_num_of_genes<-names(gene_set_list)
  each_gene_set_array<-lapply(gene_set_list,function(x){sum(x %in% rownames(count_df))})
  judge_array_all<-sum(sapply(gene_set_list,function(x){all(TRUE,x %in% rownames(count_df))}))
  #judge the result
  if (judge_array_all<length(each_gene_set_array)){
    PCA=1
  }else{PCA=0}
  return_list=list()
  return_list['PCA']<-PCA
  return_list<-append(return_list,each_gene_set_array)
  return(return_list)
}



message('usage: fileinput[absolute path],cpmornot[0,1],batchremovalornot(has to provide raw count)[combat,ruvg,not],batch_effect_file[either input a NULL or a directory to files store batcheffect],includeFF18ornot[0,1],output_directory[absolute path],folder_store_required_information[abosulte path]')

options_a <- commandArgs(trailingOnly = TRUE)

#read data and deal with gene id
df_value<-as.data.frame(read_file(options_a[1]))
gene_symbol_final<-as.character(transfer_IDs(df_value[,1]))
df_value[,1]<-gene_symbol_final
df_value<-df_value[!is.na(df_value[,1]),]
colnames(df_value)[1]<-'symbol'
#sum duplicated
df_value<-aggregate(df_value[,2:ncol(df_value)],by=list(df_value[,1]),sum)
colnames(df_value)[1]<-'symbol'
#update gene symbol to the same verion use same hgnchelper
df_value<-update_gene_symbol(df_value,1)




#judge whether user provide raw count or log2cpm value
if (options_a[2]=='0'){
  df_value_cpm<-calculate_cpm(df_value)

} else if (options_a[2]=='1'){
  print('user has provided log2cpm')
  df_value_cpm<-df_value
  df_value<-change_from_cpm_to_raw_count(df_value_cpm[,2:ncol(df_value_cpm)])
  colnames(df_value)<-colnames(df_value_cpm)[2:ncol(df_value_cpm)]
  df_value<-cbind(df_value_cpm[,1],df_value)
}

# different situations
# 1: user use combat and provide his own meta but not include FF samples
if (options_a[5]=='0'){
  print('not include FF samples')
  
if (options_a[3]=='combat' & options_a[4]=='NULL'){
  df_value_fixed<-batch_effect_removal(df_value,options_a[7],methods='combat')
  df_value_cpm<-as.data.frame(calculate_cpm(df_value_fixed))
  df_value_cpm<-df_value_cpm[,1:(ncol(df_value_cpm)-18)]
  
}else if (options_a[3]=='combat' & !options_a[4]=='NULL'){
  #the batch file should not have anyheader
  batch_effect<-fread(options_a[4],header=F)$V2
  df_value_fixed<-batch_effect_removal(df_value,options_a[7],methods='combat',batch_effect=batch_effect)
  df_value_cpm<-calculate_cpm(df_value_fixed)
  df_value_cpm<-df_value_cpm[,1:ncol(df_value_cpm)-18]
  
}else if (options_a[3]=='ruvg'){
  df_value_fixed<-batch_effect_removal(df_value,options_a[7],methods='ruvg')
  df_value_cpm<-calculate_cpm(df_value_fixed)
  df_value_cpm<-df_value_cpm[,1:ncol(df_value_cpm)-18]
  
}}else if (options_a[5]=='1' & options_a[4]=='NULL' & options_a[3]!='ruvg'){
    df_value_fixed<-batch_effect_removal(df_value,options_a[7],methods='combat')
    df_value_cpm<-calculate_cpm(df_value_fixed)
    
}else if (options_a[5]=='1' & !options_a[4]=='NULL'){
    #make the FF samples meta and add to meta files
    FF_meta<-fread(paste0(options_a[7],'/FF_meta.csv'),header=F)
    batch_effect<-fread(options_a[4],header=F)
    new_meta_final<-as.data.frame(rbind(FF_meta,batch_effect))$V2
    FF_df<-read.csv(paste0(options_a[7],'/FF_df_symbol_corrected.csv'),row.names = NULL,check.names = F)
    colnames(df_value)[1]<-'symbol'
    final_df<-merge(FF_df,df_value,by='symbol')
    df_value_fixed<-batch_effect_removal(final_df,options_a[7],methods='combat',batch_effect=new_meta_final)
    df_value_cpm<-calculate_cpm(df_value_fixed)
    
}else if (options_a[5]=='1' & options_a[4]=='NULL' & options_a[3]=='ruvg'){
    FF_df<-read.csv(paste0(options_a[7],'/FF_df_symbol_corrected.csv'),row.names = 1,check.names = F)
    colnames(df_value)[1]<-'symbol'
    final_df<-merge(FF_df,df_value,by='symbol')
    df_value_fixed<-batch_effect_removal(final_df,options_a[7],methods='ruvg')
    df_value_cpm<-calculate_cpm(df_value_fixed)
    
}
df_value_cpm<-as.data.frame(df_value_cpm)

list_judgement<-whether_perform_PCA_for_next_step(df_value_cpm,paste0(options_a[7],'/gene_sets.json'))
jsonData <- toJSON(list_judgement)
write(jsonData, paste0(options_a[6],"/PCA_related.json"))
#print(df_value_cpm)
write.csv(df_value_cpm,paste0(options_a[6],'/cpm_input_classifier.csv'))
