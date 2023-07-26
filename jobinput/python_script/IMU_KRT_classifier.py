'''
the current training model is from sklearn 0.24.2


'''
import sklearn
import pickle
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import normalize
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
import seaborn as sns
import numpy as np
import collections
import matplotlib.pyplot as plt
import argparse
from matplotlib.colors import ListedColormap
import sys
import warnings


def data_loading(dir_store_data):
    with open(f'{dir_store_data}/list_of_model_rf_knn_gnb_svm_elastic.pkl','rb') as f:
        all_models_list=pickle.load(f)
    with open(f'{dir_store_data}/list_of_gene_sets_used.pkl','rb') as f:
        all_gene_set_list=pickle.load(f)
    with open(f'{dir_store_data}/list_of_pca.pkl','rb') as f:
        all_PCA_model_list=pickle.load(f)
    return all_models_list, all_gene_set_list,all_PCA_model_list

def read_data_and_extract_matrix(gene_set_list,input_df_file):
    '''
    input df should be gene symbols*samples,csv file,with the first row as samples and the first column is gene symbol log2cpm value
    '''
    df_for_classifying_list=[]
    df_input=pd.read_csv(input_df_file,index_col=0).T.sort_index()
    all_samples=df_input.index
    # go though all the list and extract the df
    for i in gene_set_list:
        print(f'gene set list right now has {len(i)} genes')
        print('since not all genes in gene sets being included,pick genes as much as possible then perform simple imputation')
        #build a df with na values
        df_output = pd.DataFrame(np.nan, index=all_samples, columns=i)
        overlap_part=set(i) &set(df_input.columns)
        print(f'user provided matirx include {len(overlap_part)} genes in gene set')
        #judge whether user provided gene is less than 95% percent of genes we used for model trainig
        #if less, print a warning.
        if len(overlap_part)<len(i):
            warnings.warn("User provided gene set is less than 95% of genes we used for model training, please check your gene set and make sure it is correct.")

        df_for_classifying=df_input[overlap_part]
        df_for_classifying_new=pd.DataFrame(df_for_classifying,index=df_for_classifying.index,columns=df_for_classifying.columns)
        df_output=df_output.combine_first(df_for_classifying_new)
        df_output=df_output[i]
        #imputation
        imp = SimpleImputer(missing_values=np.nan, strategy='mean')
        df_output_new=imp.fit_transform(df_output.T)
        print(df_output_new.T.shape)
        sc = StandardScaler()
        features_scale = sc.fit_transform(df_output_new.T)
        df_for_classifying_list.append(features_scale)
    return df_input,df_for_classifying_list,all_samples



def pick_most_frequent_result_show_count(array1):
    Counter_part=collections.Counter(array1)
    most_frequent_subtype=Counter_part.most_common()[0][0]

    count_string=';'.join([str(j) for i,j in dict(Counter_part).items()])
    #to avoid numpy delete information in later steps
    count_string=count_string+"  "
    return count_string,most_frequent_subtype



def predict_by_all_models_and_summary_result(model_list,a,sample_array,PCA_model_list,gene_set_interest,PCA_df=False,model_interest=['rf','knn','gnb','svm','elasticnet']):
    '''
    the model list is composed of multiple lists, the first is for not PCA and the second is for PCA-based training.
    Each list represents a ML method, and the gene set list is mapped to the model list
    some model return 1,0.
    1 is KRT and 0 is IMU
    '''
    print('start classifying')
    #all models and gene set
    model_used=['rf','knn','gnb','svm','elasticnet']
    gene_set_used=gene_set_interest
    #get the index of those we want to use
    model_used_index=[i for i, item in enumerate(model_used) if item in set(model_interest)]
    gene_set_used_index=[i for i, item in enumerate(gene_set_used) if item in set(gene_set_interest)]
    #obtain the new gene set
    gene_set_used_new=[a[i] for i in gene_set_used_index]
    predict_result_list=[]

    if PCA_df:

        model_list_new=[model_list[i] for i in range(1,len(model_list),2)]
        model_list_used_new=[model_list_new[i] for i in model_used_index]
        PCA_model_list_new=[PCA_model_list[i] for i in model_used_index]
        for i in range(len(model_list_used_new)):
            for j in gene_set_used_index:
                # use previous PCA model to perform PCA
                PCA_input=PCA_model_list_new[i][j].transform(a[j])
                predict_results=model_list_used_new[i][j].predict(PCA_input)
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
    pathways_for_repeat=['immune response','Mesenchymal differentiation','Keratinization','Oxidation-reduction process','cell adhension']
    pathways_for_genes=list(np.repeat(pathways_for_repeat, [8,4,5,4,4], axis=0))
    df_input=pd.read_csv(input_df_file,index_col=0).T.sort_index()
    overlap_genes_temp=df_input.columns.intersection(genes_for_heatmap)
    final_genes_list=sorted(overlap_genes_temp,key=genes_for_heatmap.index)
    heatmap_part=df_input[final_genes_list].T
    result_subtypes_sort=sample_subtype_array[:,sample_subtype_array[1, :].argsort()]
    heatmap_for_plot_final=heatmap_part[result_subtypes_sort[0]]
    place_add_line=np.where(result_subtypes_sort[1]== 'KRT')[0][0]
    #heatmap_for_plot_final.columns=result_subtypes_sort[1]
    heatmap_for_plot_final_nor=heatmap_for_plot_final.sub(heatmap_for_plot_final.mean(axis=1), axis=0)
    #print(heatmap_for_plot_final)
    #heatmap_for_plot_final_nor=normalize(heatmap_for_plot_final,axis=1,norm='l1')
    #print(np.mean(heatmap_for_plot_final_nor,axis=1))
    #print(np.mean(heatmap_for_plot_final_nor,axis=0))
    columns_1=[[f'Sample ID: {i}' for i in heatmap_for_plot_final.columns],[f'Subtype: {j}' for j in result_subtypes_sort[1]]]
    columns_tuple=list(zip(*columns_1))
    columns_multiple=pd.MultiIndex.from_tuples(columns_tuple)
    indexs_1=[[f'Genes: {i}' for i in heatmap_for_plot_final.index],[f'Pathways: {j}' for j in pathways_for_genes]]
    indexs_tuple=list(zip(*indexs_1))
    indexs_multiple=pd.MultiIndex.from_tuples(indexs_tuple)

    heatmap_for_plot_final_nor_output=pd.DataFrame(np.array(heatmap_for_plot_final_nor),index=indexs_multiple,columns=columns_multiple)
    #make the heatmap_for_plot_final_nor to fit the https://github.com/MaayanLab/clustergrammer-py




    #fig = plt.figure(figsize=(25,25))

    heatmap_temp=sns.heatmap(heatmap_for_plot_final_nor,cmap='coolwarm',center=0,vmax=5,vmin=-5)
    heatmap_temp.vlines([place_add_line],*heatmap_temp.get_ylim())
    #axr.tick_params(right=True, top=True,bottom=True,labeltop=True,labelright=True,labelbottom=True,rotation=0)
    #axr
    plt.savefig(f'{output_dir}/heatmap_for_clustering_result.png',dpi=300)

    #fig1=plt.figure(figsize=(10,10))
    for i in range(len(df_for_classifying_list_PCA_version)):
        plt.subplot(2,2,i+1)
        sns.scatterplot(x=df_for_classifying_list_PCA_version[i][:,0], y=df_for_classifying_list_PCA_version[i][:,1],hue=sample_subtype_array[1,:])
        plt.title(f'PCA plot of different num of genes model {i}')
        plt.xlabel('PC1')
        # Set y-axis label
        plt.ylabel('PC2')


    plt.savefig(f'{output_dir}/PCA_classification_result_of_IMU_KRT.png',dpi=300)

    if another_color is not None:
        #fig2=plt.figure(figsize=(10,10))
        for i in range(len(df_for_classifying_list_PCA_version)):
            plt.subplot(2,2,i+1)
            sns.scatterplot(x=df_for_classifying_list_PCA_version[i][:,0], y=df_for_classifying_list_PCA_version[i][:,1],hue=another_color)
            plt.title(f'PCA plot of different num of genes model {i}')
            plt.xlabel('PC1')
            # Set y-axis label
            plt.ylabel('PC2')
        plt.savefig(f'{output_dir}/PCA_classification_result_of_other.png',dpi=300)


    return heatmap_for_plot_final_nor_output

def runfile(workdir, pca, input,output_dir,imumodel,gene_set):
    #dir  workdir
    #pCA PCA pca
    #log2cpmmatrix:input:input
    #output_dir:output_dir:output_dir
    #models:model:models
    #genes :gene_set:gene_set
    models,gene_sets,PCA_models=data_loading(workdir)
    if imumodel==None:
        final_result_model=''
    else:
        final_result_model=imumodel
    if gene_set==None:
        final_result_gene_set=''
    else:
        final_result_gene_set=gene_set
    if workdir and pca and input:
        print('PCA start')
        model_used_final_list=final_result_model.split(',')
        gene_set_used_final_list=final_result_gene_set.split(',')
        original_df,a,samples=read_data_and_extract_matrix(gene_sets,input)
        if (final_result_model !='') & (final_result_gene_set!=''):
            result_subtypes=predict_by_all_models_and_summary_result(models,a,samples,PCA_models,PCA_df=True,model_interest=model_used_final_list,gene_set_interest=gene_set_used_final_list)
        elif (final_result_model!='') & (final_result_gene_set==''):
            result_subtypes=predict_by_all_models_and_summary_result(models,a,samples,PCA_models,PCA_df=True,model_interest=model_used_final_list)
        elif (final_result_model=='') & (final_result_gene_set!=''):
            result_subtypes=predict_by_all_models_and_summary_result(models,a,samples,PCA_models,PCA_df=True,gene_set_interest=gene_set_used_final_list)
        elif (final_result_model=='') & (final_result_gene_set==''):
            result_subtypes=predict_by_all_models_and_summary_result(models,a,samples,PCA_models,PCA_df=True,gene_set_interest=[len(i) for i in gene_sets])
        #use the first PCA to train rf
        PCA_a=[]
        for j in range(len(a)):
            pca_temp=PCA_models[0][j]
            df_for_classifying_PCA=pca_temp.transform(np.array(a[j]))
            PCA_a.append(df_for_classifying_PCA)
        # temp_df=pd.DataFrame(result_subtypes).T
        # IMU_predicted_samples=temp_df[temp_df[1]=='IMU'][0]
        # KRT_predicted_samples=temp_df[temp_df[1]=='KRT'][0]
        another_color=None
        #result_subtypes[1:,]=np.where((result_subtypes[1:,]=='IMU'), 'KRT', 'IMU')
        # heatmap_data=visulization_part_for_provement(result.input,result_subtypes,a,result.output_dir,another_color)
        # IMU_heatmap_gene_data=heatmap_data[IMU_predicted_samples]
        # KRT_heatmap_gene_data=heatmap_data[KRT_predicted_samples]
        # high_IMU_genes=['SPRR3','TGM1','CDH3','SFN','TP63','MAPK14','AKT1','MAOA','PPARD','CDH1','TJP1','DSG3','KRT16']
        # IMU_heatmap_gene_data_high=IMU_heatmap_gene_data[IMU_heatmap_gene_data.index.isin(high_IMU_genes)]
        # KRT_heatmap_gene_data_low=KRT_heatmap_gene_data[KRT_heatmap_gene_data.index.isin(high_IMU_genes)]
        # IMU_mean=IMU_heatmap_gene_data.mean().mean()
        # KRT_mean=KRT_heatmap_gene_data_low.mean().mean()
        # if IMU_mean>KRT_mean:
        #     result_subtypes[1,]=np.where((result_subtypes[1,]=='IMU'), 'KRT', 'IMU')
        heatmap_data=visulization_part_for_provement(input,result_subtypes,PCA_a,output_dir,another_color)
        result_subtypes=pd.DataFrame(result_subtypes)
        result_subtypes.to_csv(f'{output_dir}/meta_heatmap.csv')
        result_subtypes.T.to_csv(f'{output_dir}/predict_result.txt',index=None,sep='\t')
        heatmap_data.to_csv(f'{output_dir}/heatmap_data.txt',sep='\t')
        # else:
        #     result_subtypes=pd.DataFrame(result_subtypes)
        #     result_subtypes.to_csv(f'{result.output_dir}/meta_heatmap.csv')
        #     result_subtypes.T.to_csv(f'{result.output_dir}/predict_result.txt',index=None,sep='\t')
        #     heatmap_data.to_csv(f'{result.output_dir}/heatmap_data.csv')
    else:
        print('non PCA start')
        model_used_final_list=final_result_model.split(',')
        gene_set_used_final_list=final_result_gene_set.split(',')
        original_df,a,samples=read_data_and_extract_matrix(gene_sets,input)
        if (final_result_model !='') & (final_result_gene_set!=''):
            result_subtypes=predict_by_all_models_and_summary_result(models,a,samples,PCA_models,PCA_df=False,model_interest=model_used_final_list,gene_set_interest=gene_set_used_final_list)
        elif (final_result_model!='') & (final_result_gene_set==''):
            result_subtypes=predict_by_all_models_and_summary_result(models,a,samples,PCA_models,PCA_df=False,model_interest=model_used_final_list)
        elif (final_result_model=='') & (final_result_gene_set!=''):
            result_subtypes=predict_by_all_models_and_summary_result(models,a,samples,PCA_models,PCA_df=False,gene_set_interest=gene_set_used_final_list)
        elif (final_result_model=='') & (final_result_gene_set==''):
            result_subtypes=predict_by_all_models_and_summary_result(models,a,samples,PCA_models,PCA_df=False,gene_set_interest=[len(i) for i in gene_sets])
        PCA_a=[]
        for j in range(len(a)):
            pca_temp=PCA_models[0][j]
            df_for_classifying_PCA=pca_temp.transform(np.array(a[j]))
            PCA_a.append(df_for_classifying_PCA)
        another_color=None
        #result_subtypes[1:,]=np.where((result_subtypes[1:,]=='IMU'), 'KRT', 'IMU')
        heatmap_data=visulization_part_for_provement(input,result_subtypes,PCA_a,output_dir,another_color)
        result_subtypes=pd.DataFrame(result_subtypes)
        result_subtypes.to_csv(f'{output_dir}/meta_heatmap.csv')
        result_subtypes.T.to_csv(f'{output_dir}/predict_result.txt',index=None,sep='\t')
        heatmap_data.to_csv(f'{output_dir}/heatmap_data.txt',sep='\t')





if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="This is a script use well trained ML methods to classifiy HPV positive HNSCC subtypes, especially for oropharynx cancers, \
                                      current settings are using rf,knn,gnb,svm,elasticnet logistic regression, and three different gene sets from different ML methods designed for different features \
                                      two pickle files storing model and gene set lists are required. The recommendation startup is to use non-PCA model, but if genes information is not comprehensive. PCA based models should be considered \
                                      ")
    parser.add_argument('-dir',type=str,required=True,help='directory store models and gene lists in pickle dump format',dest='dir')
    parser.add_argument('-PCA',type=int,required=True,help='whether use PCA based model',dest='PCA')
    parser.add_argument('-log2cpmmatrix',type=str,required=True,help='input file, has to be csv file and log2cpm value. the first row is sample ID and the first column is gene symbols',dest='input')
    parser.add_argument('-output_dir',type=str,required=True,help='output directory',dest='output_dir')
    parser.add_argument('-models',type=str,help='specific model applied(rf,knn,gnb,svm,elasticnet), multiple ones use , to separate.',dest='model')
    parser.add_argument('-genes',type=str,help='specific gene set interest(10,29,165), multiple ones use , to separate.',dest='gene_set')
    result=parser.parse_args()
    runfile(result.dir, result.PCA, result.input,result.output_dir,result.model, result.gene_set)


    # models,gene_sets,PCA_models=data_loading(result.dir)
    # if result.model==None:
    #     final_result_model=''
    # else:
    #     final_result_model=result.model
    # if result.gene_set==None:
    #     final_result_gene_set=''
    # else:
    #     final_result_gene_set=result.gene_set

    # if result.dir and result.PCA and result.input:
    #     print('PCA start')
    #     model_used_final_list=final_result_model.split(',')
    #     gene_set_used_final_list=final_result_gene_set.split(',')
    #     original_df,a,samples=read_data_and_extract_matrix(gene_sets,result.input)
    #     if (final_result_model !='') & (final_result_gene_set!=''):
    #         result_subtypes=predict_by_all_models_and_summary_result(models,a,samples,PCA_models,PCA_df=True,model_interest=model_used_final_list,gene_set_interest=gene_set_used_final_list)
    #     elif (final_result_model!='') & (final_result_gene_set==''):
    #         result_subtypes=predict_by_all_models_and_summary_result(models,a,samples,PCA_models,PCA_df=True,model_interest=model_used_final_list)
    #     elif (final_result_model=='') & (final_result_gene_set!=''):
    #         result_subtypes=predict_by_all_models_and_summary_result(models,a,samples,PCA_models,PCA_df=True,gene_set_interest=gene_set_used_final_list)

    #     elif (final_result_model=='') & (final_result_gene_set==''):
    #         result_subtypes=predict_by_all_models_and_summary_result(models,a,samples,PCA_models,PCA_df=True,gene_set_interest=[len(i) for i in gene_sets])

    #     #use the first PCA to train rf
    #     PCA_a=[]
    #     for j in range(len(a)):
    #         pca_temp=PCA_models[0][j]
    #         df_for_classifying_PCA=pca_temp.transform(np.array(a[j]))
    #         PCA_a.append(df_for_classifying_PCA)

    #     # temp_df=pd.DataFrame(result_subtypes).T
    #     # IMU_predicted_samples=temp_df[temp_df[1]=='IMU'][0]
    #     # KRT_predicted_samples=temp_df[temp_df[1]=='KRT'][0]
    #     another_color=None
    #     #result_subtypes[1:,]=np.where((result_subtypes[1:,]=='IMU'), 'KRT', 'IMU')
    #     # heatmap_data=visulization_part_for_provement(result.input,result_subtypes,a,result.output_dir,another_color)
    #     # IMU_heatmap_gene_data=heatmap_data[IMU_predicted_samples]
    #     # KRT_heatmap_gene_data=heatmap_data[KRT_predicted_samples]
    #     # high_IMU_genes=['SPRR3','TGM1','CDH3','SFN','TP63','MAPK14','AKT1','MAOA','PPARD','CDH1','TJP1','DSG3','KRT16']
    #     # IMU_heatmap_gene_data_high=IMU_heatmap_gene_data[IMU_heatmap_gene_data.index.isin(high_IMU_genes)]
    #     # KRT_heatmap_gene_data_low=KRT_heatmap_gene_data[KRT_heatmap_gene_data.index.isin(high_IMU_genes)]
    #     # IMU_mean=IMU_heatmap_gene_data.mean().mean()
    #     # KRT_mean=KRT_heatmap_gene_data_low.mean().mean()
    #     # if IMU_mean>KRT_mean:
    #     #     result_subtypes[1,]=np.where((result_subtypes[1,]=='IMU'), 'KRT', 'IMU')
    #     heatmap_data=visulization_part_for_provement(result.input,result_subtypes,PCA_a,result.output_dir,another_color)
    #     result_subtypes=pd.DataFrame(result_subtypes)
    #     result_subtypes.to_csv(f'{result.output_dir}/meta_heatmap.csv')
    #     result_subtypes.T.to_csv(f'{result.output_dir}/predict_result.txt',index=None,sep='\t')
    #     heatmap_data.to_csv(f'{result.output_dir}/heatmap_data.txt',sep='\t')
    #     # else:
    #     #     result_subtypes=pd.DataFrame(result_subtypes)
    #     #     result_subtypes.to_csv(f'{result.output_dir}/meta_heatmap.csv')
    #     #     result_subtypes.T.to_csv(f'{result.output_dir}/predict_result.txt',index=None,sep='\t')
    #     #     heatmap_data.to_csv(f'{result.output_dir}/heatmap_data.csv')


    # else:
    #     print('non PCA start')
    #     model_used_final_list=final_result_model.split(',')
    #     gene_set_used_final_list=final_result_gene_set.split(',')
    #     original_df,a,samples=read_data_and_extract_matrix(gene_sets,result.input)
    #     if (final_result_model !='') & (final_result_gene_set!=''):
    #         result_subtypes=predict_by_all_models_and_summary_result(models,a,samples,PCA_models,PCA_df=False,model_interest=model_used_final_list,gene_set_interest=gene_set_used_final_list)
    #     elif (final_result_model!='') & (final_result_gene_set==''):
    #         result_subtypes=predict_by_all_models_and_summary_result(models,a,samples,PCA_models,PCA_df=False,model_interest=model_used_final_list)
    #     elif (final_result_model=='') & (final_result_gene_set!=''):
    #         result_subtypes=predict_by_all_models_and_summary_result(models,a,samples,PCA_models,PCA_df=False,gene_set_interest=gene_set_used_final_list)
    #     elif (final_result_model=='') & (final_result_gene_set==''):
    #         result_subtypes=predict_by_all_models_and_summary_result(models,a,samples,PCA_models,PCA_df=False,gene_set_interest=[len(i) for i in gene_sets])
    #     PCA_a=[]
    #     for j in range(len(a)):
    #         pca_temp=PCA_models[0][j]
    #         df_for_classifying_PCA=pca_temp.transform(np.array(a[j]))
    #         PCA_a.append(df_for_classifying_PCA)
    #     another_color=None
    #     #result_subtypes[1:,]=np.where((result_subtypes[1:,]=='IMU'), 'KRT', 'IMU')
    #     heatmap_data=visulization_part_for_provement(result.input,result_subtypes,PCA_a,result.output_dir,another_color)
    #     result_subtypes=pd.DataFrame(result_subtypes)
    #     result_subtypes.to_csv(f'{result.output_dir}/meta_heatmap.csv')
    #     result_subtypes.T.to_csv(f'{result.output_dir}/predict_result.txt',index=None,sep='\t')
    #     heatmap_data.to_csv(f'{result.output_dir}/heatmap_data.txt')



