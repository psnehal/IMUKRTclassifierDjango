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

#' calculate the cpm
#' @param gene_counts count df with rows as genes and columns as samples, the first column should be genes,colnames is samples ID

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

#' batch effect removal
#' @param count_df count df with rows as genes and columns as samples, the first column should be genes,colnames is samples ID
#'
batch_effect_removal<-function(count_df,dir_store_files,methods='combat',batch_effect=NULL){
  if (methods=='ruvg'){
  #ruvg
  library(RUVSeq)
  IMU_KRT_not_sig_genes<-read.csv(paste0(dir_store_files,'/IMU_KRT_not_sig_genes.csv'),row.names = 1)
  rownames(count_df)<-count_df[,1]
  set <- newSeqExpressionSet(as.matrix(count_df[,2:ncol(count_df)]))
  idx  <- rowSums(counts(set) > 5) >= 2
  set  <- set[idx, ]
  temp_genes<-IMU_KRT_not_sig_genes[,1][IMU_KRT_not_sig_genes[,1] %in% rownames(set)]
  set2<-RUVg(set,temp_genes,k=1)
  Ruvg_imporved<-as.data.frame(set2@assayData$normalizedCounts)
  Ruvg_imporved<-cbind(rownames(Ruvg_imporved),Ruvg_imporved)
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



message('usage: fileinput[absolute path],cpmornot[0,1],batchremovalornot(has to provide raw count)[combat,ruvg,not],batch_effect_file[either input a NULL or a directory to files store batcheffect],includeFF18ornot[0,1],output_directory[absolute path],folder_store_required_information[abosulte path]')
options_a <- commandArgs(trailingOnly = TRUE)

df_value<-read_file(options_a[1])

df_value<-as.data.frame(df_value)

df_value_cpm<-calculate_cpm(df_value)



# different situations
# 1: user use combat and provide his own meta but not include FF samples
if (options_a[5]=='0'){
  print('not include FF samples')
  
if (options_a[3]=='combat' & options_a[4]=='NULL'){
  df_value_fixed<-batch_effect_removal(df_value,options_a[7],methods='combat')
  df_value_cpm<-calculate_cpm(df_value_fixed)
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
write.csv(df_value_cpm,paste0(options_a[6],'/cpm_input_classifier.csv'),row.names=F)
