'''
the current training model is from sklearn 0.24.2


'''
import sklearn
import pickle
import pandas as pd
from sklearn.decomposition import PCA
import seaborn as sns
import numpy as np
import collections
from sklearn.preprocessing import normalize
import matplotlib.pyplot as plt
import argparse
from matplotlib.colors import ListedColormap
import sys
def data_loading(dir_store_data):
  with open(f'{dir_store_data}/list_of_model_rf_knn_gnb_svm_elastic.pkl','rb') as f:
    all_models_list=pickle.load(f)
  with open(f'{dir_store_data}/list_of_gene_sets_used.pkl','rb') as f:
    all_gene_set_list=pickle.load(f)
  return all_models_list, all_gene_set_list

def read_data_and_extract_matrix(gene_set_list,input_df_file,PCA_version=False):
  '''
  input df should be gene symbols*samples,csv file,with the first row as samples and the first column is gene symbol log2cpm value
  '''

  df_for_classifying_list=[]
  if input_df_file.endswith('.csv'):
      df_input=pd.read_csv(input_df_file,index_col=0).T.sort_index()
  elif input_df_file.endswith('.tsv'):
      df_input=pd.read_csv(input_df_file,index_col=0,sep='\t').T.sort_index()
  all_samples=df_input.index
  # go though all the list and extract the df
  for i in gene_set_list:
    print(f'gene set list right now has {len(i)} genes')
    if PCA_version:
      print('since not all genes in gene sets being included,pick genes as much as possible and perform PCA and classify by PCA results, but the result may not be as accuract as original one')

      overlap_part=set(i) &set(df_input.columns)
      df_for_classifying=df_input[overlap_part]
      print(f'gene set with {len(i)} genes has {df_for_classifying.shape[1]} of genes being used to perform PCA')
      df_for_classifying_PCA=PCA(n_components=2).fit_transform(np.array(df_for_classifying))
      df_for_classifying_list.append(df_for_classifying_PCA)
    else:
      try:
        print('use not PCA perform classify')
        df_for_classifying=df_input[i]
        df_for_classifying_list.append(df_for_classifying)
      except KeyError as err:
        print('cant get all genes in this gene sets')
        print(err)
        exit()
  return df_input,df_for_classifying_list,all_samples


def pick_most_frequent_result_show_count(array1):
  Counter_part=collections.Counter(array1)
  most_frequent_subtype=Counter_part.most_common()[0][0]
  count_string='|'.join([str(j) for i,j in dict(Counter_part).items()])
  return count_string,most_frequent_subtype



def predict_by_all_models_and_summary_result(model_list,a,sample_array,PCA_df=False,model_interest=['rf','knn','gnb','svm','elasticnet'],gene_set_interest=['7','168','960']):
  '''
  the model list is composed of multiple lists, the first is for not PCA and the second is for PCA-based training.
  Each list represents a ML method, and the gene set list is mapped to the model list
  some model return 1,0.
  1 is KRT and 0 is IMU
  '''
  print('start classifying')
  #all models and gene set
  model_used=['rf','knn','gnb','svm','elasticnet']
  gene_set_used=['7','168','960']
  #get the index of those we want to use
  model_used_index=[i for i, item in enumerate(model_used) if item in set(model_interest)]
  gene_set_used_index=[i for i, item in enumerate(gene_set_used) if item in set(gene_set_interest)]
  #obtain the new gene set
  gene_set_used_new=[a[i] for i in gene_set_used_index]
  predict_result_list=[]

  if PCA_df:
    model_list_new=[model_list[i] for i in range(1,len(model_list),2)]
    model_list_used_new=[model_list_new[i] for i in model_used_index]

    for i in range(len(model_list_used_new)):
      for j in gene_set_used_index:
        predict_results=model_list_used_new[i][j].predict(a[j])
        predict_result_list.append(predict_results)


  else:
    model_list_new=[model_list[i] for i in range(0,len(model_list),2)]
    model_list_used_new=[model_list_new[i] for i in model_used_index]
    for i in range(len(model_list_used_new)):
      for j in gene_set_used_index:
        predict_results=model_list_used_new[i][j].predict(a[j])
        predict_result_list.append(predict_results)

  predict_array=np.array(predict_result_list)
  predict_array[predict_array==0]='IMU'
  predict_array[predict_array==1]='KRT'
  count_string,result_subtypes=np.apply_along_axis(pick_most_frequent_result_show_count,0,predict_array)
  sample_subtype_array=np.stack((sample_array,result_subtypes,count_string))
  return sample_subtype_array

def visulization_part_for_provement(input_df_file,sample_subtype_array,df_for_classifying_list_PCA_version,output_dir,another_color=None):
  '''
  heatmaps and PCA plot will be provided, genes used are coming from yanxiao's paper
  '''
  genes_for_heatmap=['CD40','CD4','IL7R','HLA-DQB2','KIT','FOXP3','NFKB2','BCL2','VIM','SNAI1','ZEB1','SPRR3','TGM1','CDH3','SFN','TP63','MAPK14','AKT1','MAOA','PPARD','CDH1','TJP1','DSG3','KRT16']
  df_input=pd.read_csv(input_df_file,index_col=0).T.sort_index()
  overlap_genes_temp=df_input.columns.intersection(genes_for_heatmap)
  final_genes_list=sorted(overlap_genes_temp,key=genes_for_heatmap.index)
  heatmap_part=df_input[final_genes_list].T
  result_subtypes_sort=sample_subtype_array[:,sample_subtype_array[1, :].argsort()]
  heatmap_for_plot_final=heatmap_part[result_subtypes_sort[0]]
  #heatmap_for_plot_final.columns=result_subtypes_sort[1]
  heatmap_for_plot_final_nor=heatmap_for_plot_final.sub(heatmap_for_plot_final.mean(axis=1), axis=0)
  #print(heatmap_for_plot_final)
  #heatmap_for_plot_final_nor=normalize(heatmap_for_plot_final,axis=1,norm='l1')
  #print(np.mean(heatmap_for_plot_final_nor,axis=1))
  #print(np.mean(heatmap_for_plot_final_nor,axis=0))
  heatmap_for_plot_final_nor=pd.DataFrame(heatmap_for_plot_final_nor,index=heatmap_for_plot_final.index,columns=heatmap_for_plot_final.columns)
  fig = plt.figure(figsize=(15,15))
  sns.heatmap(heatmap_for_plot_final_nor,cmap='coolwarm',center=0,vmax=5,vmin=-5)
  #axr.tick_params(right=True, top=True,bottom=True,labeltop=True,labelright=True,labelbottom=True,rotation=0)
  #axr
  plt.savefig(f'{output_dir}/heatmap for clustering result.png',dpi=300)

  fig1=plt.figure(figsize=(10,10))
  for i in range(len(df_for_classifying_list_PCA_version)):
    plt.subplot(2,2,i+1)
    sns.scatterplot(x=df_for_classifying_list_PCA_version[i][:,0], y=df_for_classifying_list_PCA_version[i][:,1],hue=sample_subtype_array[1,:])
    plt.title(f'PCA plot of different num of genes model {i}')
    plt.xlabel('PC1')
    # Set y-axis label
    plt.ylabel('PC2')


  plt.savefig(f'{output_dir}/PCA classification result of IMU KRT.png',dpi=300)

  if another_color is not None:
        fig2=plt.figure(figsize=(10,10))
        for i in range(len(df_for_classifying_list_PCA_version)):
          plt.subplot(2,2,i+1)
          sns.scatterplot(x=df_for_classifying_list_PCA_version[i][:,0], y=df_for_classifying_list_PCA_version[i][:,1],hue=another_color)
          plt.title(f'PCA plot of different num of genes model {i}')
          plt.xlabel('PC1')
          # Set y-axis label
          plt.ylabel('PC2')


        plt.savefig(f'{output_dir}/PCA classification result of other.png',dpi=300)


  return heatmap_for_plot_final_nor

if __name__ == '__main__':




    parser = argparse.ArgumentParser(description="This is a script use well trained ML methods to classifiy HPV positive HNSCC subtypes, especially for oropharynx cancers, \
                                      current settings are using rf,knn,gnb,svm,elasticnet logistic regression, and three different models from different ML methods designed for different features (7,168,960 genes) \
                                      two pickle files storing model and gene set lists are required. The recommendation startup is to use non-PCA model, but if genes information is not comprehensive. PCA based models should be considered \
                                      ")
    parser.add_argument('-dir',type=str,required=True,help='directory store models and gene lists in pickle dump format',dest='dir')
    parser.add_argument('-PCA',type=int,required=True,help='whether use PCA based model',dest='PCA')
    parser.add_argument('-log2cpmmatrix',type=str,required=True,help='input file, has to be csv file and log2cpm value. the first row is sample ID and the first column is gene symbols',dest='input')
    parser.add_argument('-output_dir',type=str,required=True,help='output directory',dest='output_dir')
    parser.add_argument('-models',type=str,help='specific model applied(rf,knn,gnb,svm,elasticnet), multiple ones use , to separate.',dest='model')
    parser.add_argument('-genes',type=str,help='specific gene set interest(7,168,960), multiple ones use , to separate.',dest='gene_set')
    result=parser.parse_args()
    models,gene_sets=data_loading(result.dir)
    if result.model==None:
        final_result_model=''
    else:
        final_result_model=result.model
    if result.gene_set==None:
        final_result_gene_set=''
    else:
        final_result_gene_set=result.gene_set

    if result.dir and result.PCA and result.input:
        print('PCA start')
        model_used_final_list=final_result_model.split(',')
        gene_set_used_final_list=final_result_gene_set.split(',')
        original_df,a,samples=read_data_and_extract_matrix(gene_sets,result.input,PCA_version=True)
        if (final_result_model !='') & (final_result_gene_set!=''):
            result_subtypes=predict_by_all_models_and_summary_result(models,a,samples,PCA_df=True,model_interest=model_used_final_list,gene_set_interest=gene_set_used_final_list)
        elif (final_result_model!='') & (final_result_gene_set==''):
            result_subtypes=predict_by_all_models_and_summary_result(models,a,samples,PCA_df=True,model_interest=model_used_final_list)
        elif (final_result_model=='') & (final_result_gene_set!=''):
            result_subtypes=predict_by_all_models_and_summary_result(models,a,samples,PCA_df=True,gene_set_interest=gene_set_used_final_list)

        elif (final_result_model=='') & (final_result_gene_set==''):
            result_subtypes=predict_by_all_models_and_summary_result(models,a,samples,PCA_df=True)

        #since the PCA only capture the variance,it can separate the subtype of the HNSCC but the direction is not clear.
        #here we manually correct the predict result.
        temp_df=pd.DataFrame(result_subtypes).T
        IMU_predicted_samples=temp_df[temp_df[1]=='IMU'][0]
        KRT_predicted_samples=temp_df[temp_df[1]=='KRT'][0]
        another_color=None
        #result_subtypes[1:,]=np.where((result_subtypes[1:,]=='IMU'), 'KRT', 'IMU')
        heatmap_data=visulization_part_for_provement(result.input,result_subtypes,a,result.output_dir,another_color)
        IMU_heatmap_gene_data=heatmap_data[IMU_predicted_samples]
        KRT_heatmap_gene_data=heatmap_data[KRT_predicted_samples]
        high_IMU_genes=['SPRR3','TGM1','CDH3','SFN','TP63','MAPK14','AKT1','MAOA','PPARD','CDH1','TJP1','DSG3','KRT16']
        IMU_heatmap_gene_data_high=IMU_heatmap_gene_data[IMU_heatmap_gene_data.index.isin(high_IMU_genes)]
        KRT_heatmap_gene_data_low=KRT_heatmap_gene_data[KRT_heatmap_gene_data.index.isin(high_IMU_genes)]
        IMU_mean=IMU_heatmap_gene_data.mean().mean()
        KRT_mean=KRT_heatmap_gene_data_low.mean().mean()
        if IMU_mean>KRT_mean:
            result_subtypes[1,]=np.where((result_subtypes[1,]=='IMU'), 'KRT', 'IMU')
            heatmap_data=visulization_part_for_provement(result.input,result_subtypes,a,result.output_dir,another_color)
            result_subtypes=pd.DataFrame(result_subtypes)
            result_subtypes.to_csv(f'{result.output_dir}/meta_heatmap.csv')
            result_subtypes.T.to_csv(f'{result.output_dir}/predict_result.txt',index=None,sep='\t')
            heatmap_data.to_csv(f'{result.output_dir}/heatmap_data.csv')
        else:
            result_subtypes=pd.DataFrame(result_subtypes)
            result_subtypes.to_csv(f'{result.output_dir}/meta_heatmap.csv')
            result_subtypes.T.to_csv(f'{result.output_dir}/predict_result.txt',index=None,sep='\t')
            heatmap_data.to_csv(f'{result.output_dir}/heatmap_data.csv')


    else:
        print('non PCA start')
        original_df,a,samples=read_data_and_extract_matrix(gene_sets,result.input,PCA_version=False)
        model_used_final_list=final_result_model.split(',')
        gene_set_used_final_list=final_result_gene_set.split(',')
        if (final_result_model !='') & (final_result_gene_set!=''):
            result_subtypes=predict_by_all_models_and_summary_result(models,a,samples,PCA_df=False,model_interest=model_used_final_list,gene_set_interest=gene_set_used_final_list)
        elif (final_result_model!='') & (final_result_gene_set==''):
            result_subtypes=predict_by_all_models_and_summary_result(models,a,samples,PCA_df=False,model_interest=model_used_final_list)
        elif (final_result_model=='') & (final_result_gene_set!=''):
            result_subtypes=predict_by_all_models_and_summary_result(models,a,samples,PCA_df=False,gene_set_interest=gene_set_used_final_list)
        elif (final_result_model=='') & (final_result_gene_set==''):
            result_subtypes=predict_by_all_models_and_summary_result(models,a,samples,PCA_df=False)
        heatmap_data=visulization_part_for_provement(result.input,result_subtypes,a,result.output_dir,another_color)

        result_subtypes.to_csv(f'{result.output_dir}/meta_heatmap.csv')
        result_subtypes.T.to_csv(f'{result.output_dir}/predict_result.txt',index=None,sep='\t')
        heatmap_data.to_csv(f'{result.output_dir}/heatmap_data.csv')
        result_subtypes.T.to_csv(f'{result.output_dir}/predict_result.txt',index=None,sep='\t')
