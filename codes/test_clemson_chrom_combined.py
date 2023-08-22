# Test significant SNPs extracted combining all chromosomes based on MAC season 6 dat on clemson
import argparse
import os
import datetime
from datetime import date
import time
#import seaborn as sns
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
    unique_cultivar_list = list(phenotype_file_df["genotype"].unique())
    

    # Read the list of accessions from vcf files
    vcf_df_head = pd.read_csv("https://github.com/genophenoenvo/sorghum_data/releases/download/v0.0.5/sorghum.filtered.season4.season6.vcf.gz", sep="\t", skiprows=82, dtype='string', nrows=100)
    row_matrices_cultivar_list = list(vcf_df_head.columns)[9:]
    # Consider the cultivars common in both dataset
    unique_cultivar_filtered = intersection(row_matrices_cultivar_list, unique_cultivar_list)
    print("Number of cultivars under consideration:" + str(len(unique_cultivar_filtered)))
    
    # Create "target" similarity matrix for test data
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
    #cultivar_estimated_corr_df.to_csv(os.path.join(args.target_path,"clemson_target_correlation_df_"+args.phenotype+".csv"))
    cultivar_estimated_corr_df_pivot = cultivar_estimated_corr_df.pivot_table(index="CULTIVAR_1", columns="CULTIVAR_2", values="estimated_corr", fill_value=1)
    cultivar_estimated_corr_df_pivot.to_csv(os.path.join(args.target_path,"clemson_target_correlation_matrix_"+args.phenotype+".csv"))
    cultivar_estimated_corr = torch.FloatTensor(np.array(cultivar_estimated_corr_df_pivot)).to(device)

    # Read significant SNP list to be tested    
    significant_snp_list_df = pd.read_csv(os.path.join(args.snp_list_path, args.significant_snp_list))
    significant_row_matrices = pd.DataFrame()
    for chromosome_number in range(1,11):
        chromosome = str(chromosome_number).zfill(2)
        matrices_df = pd.read_pickle(f'https://storage.googleapis.com/gpe-sorghum/whole-vcf-snp-arrays/{chromosome}-matrix.pickle.gz')
        matrices_df = pd.merge(significant_snp_list_df, matrices_df, how='inner', left_on=['SNP_ID'], right_on=['ID'])
        significant_row_matrices = pd.concat([significant_row_matrices, matrices_df])
    
    significant_row_matrices = significant_row_matrices.reset_index(drop=True)
    print("Number of significant row matrices:" + str(significant_row_matrices["SNP_ID"].count()))

    # Keep only the accessions of test data in "significant" similarity matrices to compute correlation coeffcient
    cultivars_removed = list(set(row_matrices_cultivar_list) - set(unique_cultivar_list))
    ind_removed = [row_matrices_cultivar_list.index(i) for i in cultivars_removed]
    significant_row_matrices["matrix_reduced"] = significant_row_matrices.apply(lambda x: np.delete(x["matrix"], ind_removed, axis=0), axis=1)
    significant_row_matrices["matrix_reduced"] = significant_row_matrices.apply(lambda x: np.delete(x["matrix_reduced"], ind_removed, axis=1), axis=1)
    for i in ind_removed:
        row_matrices_cultivar_list.pop(i)
    if list(cultivar_estimated_corr_df_pivot.index) != row_matrices_cultivar_list:
        print("Index Not Matched.")
        exit()

    # Compute the correlation coeffcient 
    significant_row_matrices = torch.FloatTensor(np.array(list(significant_row_matrices["matrix_reduced"]))).to(device)
    print("Size of final similarity matrix: " + str(significant_row_matrices.size()))
    significant_row_matrices = significant_row_matrices.sum(dim=0)
    corr = correlation_coefficient(cultivar_estimated_corr, significant_row_matrices)
                
    print("Correlation Achieved:" + str(corr.item()))
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="phenotype-to-genotype")
    parser.add_argument('--cuda', default= True)
    parser.add_argument('--gpu_device', default="7", help="GPU devices to be used")
    
    # input details
    parser.add_argument('--phenotype_file', type=str, default='/research/iprobe-paldebas/Research Work/phenotype-to-genotype/phenotype_data/clemson_plant_height.csv')
    parser.add_argument('--phenotype', type=str, default='max_height_cm', choices=['max_growth_cm_gdd','max_height_cm'], help='PhenoCam Site ID')
    parser.add_argument('--snp_list_path', type=str, default='/research/iprobe-paldebas/Research Work/phenotype-to-genotype/SNP_outputs_v0.0.5/significant_SNPs')
    parser.add_argument('--significant_snp_list', type=str, default='Chr_combined_significant_SNP_list.csv')
    # output path details
    parser.add_argument('--target_path', type=str, default='target_correlation_matrix/Clemson')
    args = parser.parse_args()
    main(args)