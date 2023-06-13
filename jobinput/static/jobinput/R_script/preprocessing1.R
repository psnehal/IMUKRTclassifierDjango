#' read the data
#' @param file_name
#'
read_file<-function(file_name){
  library(data.table)
  # There are some trick part here
  # 1. the first column is the gene name, so we need to skip it
  # 2. the first columnnames is the gene names. We need to judge this.
  df_use<-fread(file_name)
  return(df_use)

}
calculate_cpm<-function(gene_counts){
  library(edgeR)
  gene_counts<-gene_counts[order(gene_counts[,1]),]
  gene_counts<-as.data.frame(gene_counts)
  expression_count<-DGEList(as.matrix(gene_counts[,2:ncol(gene_counts)]),genes = as.character(gene_counts[,1]))
  isexpr <- rowSums(cpm(expression_count)>=1) >=2
  expression_count<- expression_count[isexpr,]
  expression_count$samples$lib.size <- colSums(expression_count$counts)
  expression_count <- calcNormFactors(expression_count, method="upperquartile")
  final_cpm<-cpm(expression_count, normalized.lib.sizes=TRUE, log=T)
  rownames(final_cpm)<-expression_count$genes$genes
  colnames(final_cpm)<-colnames(gene_counts)[2:ncol(gene_counts)]
  #add gene symbol as the first column
  final_cpm<-cbind(rownames(final_cpm),final_cpm)
  colnames(final_cpm)[1]<-'symbol'
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



whether_perform_PCA_for_next_step<-function(count_df,json_file){
  library(rjson)
  gene_set_list <- fromJSON(file=json_file)
  all_gene_set_num_of_genes<-names(gene_set_list)
  count_df<-as.data.frame(count_df)
  each_gene_set_array<-lapply(gene_set_list,function(x){sum(x %in% count_df$symbol)})
  #all(TRUE,x %in% count_df$symbol)
  judge_array_all<-sum(sapply(gene_set_list,function(x){
    ifelse(sum(x %in% count_df$symbol)>length(x)*0.9,1,0)}))
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
if (options_a[2]=='1'){
  print('user has provided log2cpm')
  df_value_cpm<-df_value
  df_value<-change_from_cpm_to_raw_count(df_value_cpm[,2:ncol(df_value_cpm)])
  colnames(df_value)<-colnames(df_value_cpm)[2:ncol(df_value_cpm)]
  df_value<-cbind(df_value_cpm[,1],df_value)
}


list_judgement<-whether_perform_PCA_for_next_step(df_value,paste0(options_a[7],'/gene_sets.json'))

jsonData <- toJSON(list_judgement)
write(jsonData, paste0(options_a[6],"/PCA_related.json"))

df_value_cpm<-calculate_cpm(df_value)

write.csv(df_value,paste0(options_a[6],'/raw_count_input.csv'),row.names=F)
write.csv(df_value_cpm,paste0(options_a[6],'/cpm_input_classifier_no_batchre.csv'),row.names=F)
