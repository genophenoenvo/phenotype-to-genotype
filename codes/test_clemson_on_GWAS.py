# Test clemson data
import argparse
import os
import datetime
from datetime import date
import matplotlib.pyplot as plt
import pandas as pd
import tarfile
import itertools
import numpy as np
import torch

# Return intersection of two lists
def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

def correlation_coefficient(X, Y):
    denominator = X.std() * Y.std()
    numerator = torch.mean((X - X.mean()) * (Y - Y.mean()))
    result = numerator / denominator
    return result

def main(args):
    # For faster computation
    if args.cuda:
        import torch.backends.cudnn as cudnn
        os.environ['CUDA_VISIBLE_DEVICES'] = args.gpu_device
    device = torch.device('cuda' if args.cuda else 'cpu')
    
    device = torch.device("cuda" if (args.cuda) and torch.cuda.is_available() else "cpu")
    
    # Retrieve current working directory
    curr_wd = os.path.dirname(__file__)

    # Read the phenotype file into a dataframe
    phenotype_file_df=pd.read_csv(args.phenotype_file)
    phenotype_file_df = phenotype_file_df.sort_values(by=['cultivar']).rename(columns={"cultivar": "genotype"})
    phenotype_file_df = phenotype_file_df.rename(columns={"plant_height_cm": args.phenotype})
    MAC_target_corr = pd.read_csv(args.MAC_target_corr)
     
    row_matrices_df_list = []
    tar = tarfile.open(os.path.join(curr_wd, args.row_matrices_path, args.row_matrices))
    row_matrices_index_list = tar.getnames()
    row_matrices_index_list = [item.split('.')[5] for item in row_matrices_index_list]
    print("Number of row matrices:"+str(len(row_matrices_index_list)))
    for row_matrix in tar.getmembers():
        f = tar.extractfile(row_matrix)
        row_matrices_df_list.append(pd.read_csv(f, sep ='\t').set_index("Unnamed: 0"))
    tar.close()
    
    unique_cultivar_list = list(phenotype_file_df["genotype"].unique())
    row_matrices_cultivar_list = list(row_matrices_df_list[0].index)
    print("Number of unique cultivars/genptypes in row_matrices:" + str(len(row_matrices_cultivar_list)))
    unique_cultivar_filtered = intersection(row_matrices_cultivar_list, unique_cultivar_list)
    print("Number of cultivars under consideration:" + str(len(unique_cultivar_filtered)))
    
    pairwise_cultivar = list(itertools.permutations(unique_cultivar_filtered, 2))
    print("No. of pairwise cultivars:" + str(len(pairwise_cultivar)))

    cultivar_estimated_corr_df = pd.DataFrame(pairwise_cultivar, columns =['CULTIVAR_1', 'CULTIVAR_2'])
    cultivar_estimated_corr_df = pd.merge(left = cultivar_estimated_corr_df, right = phenotype_file_df[['genotype', args.phenotype]], 
                            how='inner', left_on=['CULTIVAR_1'], right_on=['genotype'])
    cultivar_estimated_corr_df = cultivar_estimated_corr_df.rename(columns={args.phenotype:args.phenotype+"_1"})
    cultivar_estimated_corr_df = pd.merge(left = cultivar_estimated_corr_df, right = phenotype_file_df[['genotype', args.phenotype]], 
                               how='inner', left_on=['CULTIVAR_2'], right_on=['genotype'])
    cultivar_estimated_corr_df = cultivar_estimated_corr_df.rename(columns={args.phenotype:args.phenotype+"_2"})
    cultivar_estimated_corr_df["difference"] = abs(cultivar_estimated_corr_df[args.phenotype+"_1"] - cultivar_estimated_corr_df[args.phenotype+"_2"])
    cultivar_estimated_corr_df = cultivar_estimated_corr_df[["CULTIVAR_1", "CULTIVAR_2", "difference"]]
    max_diff = cultivar_estimated_corr_df["difference"].max()
    min_diff = cultivar_estimated_corr_df["difference"].min()
    cultivar_estimated_corr_df["estimated_corr"] = round((cultivar_estimated_corr_df["difference"]-min_diff)/(max_diff - min_diff),2)
    cultivar_estimated_corr_df["estimated_corr"] = 1 - cultivar_estimated_corr_df["estimated_corr"]
    cultivar_estimated_corr_df.to_csv(os.path.join(args.target_path,"clemson_target_correlation_df_"+args.phenotype+".csv"))
    cultivar_estimated_corr_df_pivot = cultivar_estimated_corr_df.pivot_table(index="CULTIVAR_1", columns="CULTIVAR_2", values="estimated_corr", fill_value=1)
    cultivar_estimated_corr_df_pivot.to_csv(os.path.join(args.target_path,"clemson_target_correlation_matrix_"+args.phenotype+".csv"))
    cultivar_estimated_corr = torch.FloatTensor(np.array(cultivar_estimated_corr_df_pivot)).to(device)
    
    for i in range(0, len(row_matrices_index_list)):
        row_matrices_df_list[i] = row_matrices_df_list[i].filter(items = unique_cultivar_filtered, axis=0)
        row_matrices_df_list[i] = row_matrices_df_list[i].filter(items = unique_cultivar_filtered, axis=1)
        if list(cultivar_estimated_corr_df_pivot.index) != list(row_matrices_df_list[i].index):
            print("Index Not Matched -" + row_matrices_index_list[i])
    
    eff_snp_list = list(pd.read_csv(os.path.join(args.snp_list_path, args.effective_snp_list))["SNP_ID"])
    eff_snp_list_ind = [row_matrices_index_list.index(i) for i in eff_snp_list]
    
    row_matrices_processed = torch.FloatTensor(np.array([row_matrices_df_list[i] for i in eff_snp_list_ind])).to(device)
    row_matrices_processed = row_matrices_processed.sum(dim=0)
    corr = correlation_coefficient(cultivar_estimated_corr, row_matrices_processed)
                
    print("Correlation Achieved:" + str(corr))
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="phenotype-to-genotype")
    parser.add_argument('--cuda', default= True)
    parser.add_argument('--gpu_device', default="7", help="GPU devices to be used")
    
    # input details
    parser.add_argument('--row_matrices_path', type=str, default= 'row_matrices_v0.0.4',help = 'full path of row_matrices given in tar.gz format')
    parser.add_argument('--row_matrices', type=str, default= 'sorghum.filtered.season4.season6.max_height_cm_p0001.row_matrices.tar.gz',help = 'full path of row_matrices given in tar.gz format')
    parser.add_argument('--phenotype_file', type=str, default='/research/iprobe-paldebas/Research Work/phenotype-to-genotype/phenotype_data/clemson_plant_height.csv')
    parser.add_argument('--MAC_target_corr', type=str, default='/research/iprobe-paldebas/Research Work/phenotype-to-genotype/target_correlation_matrix/target_correlation_matrix_reducedmax_height_cm.csv')
    parser.add_argument('--phenotype', type=str, default='max_height_cm', choices=['max_growth_cm_gdd','max_height_cm'])
    parser.add_argument('--snp_list_path', type=str, default='/research/iprobe-paldebas/Research Work/phenotype-to-genotype/SNP_outputs_v0.0.4/significant_SNPs')
    parser.add_argument('--effective_snp_list', type=str, default='max_height_cm_p0001_significant_SNP_list.csv')
    # output path details
    parser.add_argument('--target_path', type=str, default='target_correlation_matrix/Clemson')
    args = parser.parse_args()
    main(args)