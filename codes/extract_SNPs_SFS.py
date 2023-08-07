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

def correlation_coefficient(target_corr, row_matrices):
    target_corr = target_corr.repeat(row_matrices.size(0),1,1)
    denominator = target_corr.std(dim=[1,2]) * row_matrices.std(dim=[1,2])
    target_corr = target_corr - target_corr.mean(dim=[1,2]).unsqueeze(1).unsqueeze(2)
    row_matrices = row_matrices - row_matrices.mean(dim=[1,2]).unsqueeze(1).unsqueeze(2)
    torch.cuda.empty_cache()
    target_corr = torch.mean(target_corr * row_matrices, dim=[1,2])
    result = target_corr / denominator
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
    
    # Calculate number of unique cultivars in the phenotype file
    unique_cultivar_list = list(phenotype_file_df["genotype"].unique())
    print("Number of unique cultivars/genptypes in phenotype file:" + str(len(unique_cultivar_list)))

    row_matrices_df_list = []
    tar = tarfile.open(os.path.join(curr_wd, args.row_matrices_path, args.row_matrices))
    row_matrices_index_list = tar.getnames()
    row_matrices_index_list = [item.split('.')[5] for item in row_matrices_index_list]
    print("Number of row matrices:"+str(len(row_matrices_index_list)))
    for row_matrix in tar.getmembers():
        f = tar.extractfile(row_matrix)
        row_matrices_df_list.append(pd.read_csv(f, sep ='\t').set_index("Unnamed: 0"))
    tar.close()
    
    row_matrices_cultivar_list = list(row_matrices_df_list[0].index)

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
    cultivar_estimated_corr_df_pivot = cultivar_estimated_corr_df.pivot_table(index="CULTIVAR_1", columns="CULTIVAR_2", values="estimated_corr", fill_value=1)
    cultivar_estimated_corr_df_pivot.to_csv(os.path.join(args.target_path,"target_correlation_matrix_"+args.phenotype+".csv"))

    for i in range(0, len(row_matrices_index_list)):
        if list(cultivar_estimated_corr_df_pivot.index) != list(row_matrices_df_list[i].index):
            print("Index Not Matched -" + row_matrices_index_list[i])
    
    cultivar_estimated_corr = torch.FloatTensor(np.array(cultivar_estimated_corr_df_pivot)).to(device)
    del cultivar_estimated_corr_df_pivot
    del cultivar_estimated_corr_df
    row_matrices = torch.FloatTensor(np.array(row_matrices_df_list)).to(device)
    start_time = time.time()
    snp_list = []
    snp_list_ind = []
    corr_coeff = []
    corr_coeff_temp = correlation_coefficient(cultivar_estimated_corr, row_matrices)
    snp_list.append(row_matrices_index_list[torch.argmax(corr_coeff_temp).item()])
    snp_list_ind.append(torch.argmax(corr_coeff_temp).item())
    corr_coeff.append(torch.max(corr_coeff_temp).item())
    
    del row_matrices
    del corr_coeff_temp
    
    if args.num_iterations is None:
        args.num_iterations = len(row_matrices_index_list)
    
    while len(snp_list) < args.num_iterations:
        row_matrices_processed = torch.FloatTensor(np.array([row_matrices_df_list[i] for i in snp_list_ind])).to(device)
        row_matrices_processed = row_matrices_processed.sum(dim=0)
        row_matrices_left_ind = list(set(np.arange(0,len(row_matrices_index_list))) - set(snp_list_ind))
        row_matrices_left = torch.FloatTensor(np.array([row_matrices_df_list[i] for i in row_matrices_left_ind])).to(device)
        row_matrices_left = row_matrices_left + row_matrices_processed.repeat(row_matrices_left.size(0),1,1) 
        del row_matrices_processed
        corr = correlation_coefficient(cultivar_estimated_corr, row_matrices_left)
        snp_list.append(row_matrices_index_list[row_matrices_left_ind[torch.argmax(corr).item()]])
        snp_list_ind.append(row_matrices_left_ind[torch.argmax(corr).item()])
        corr_coeff.append(torch.max(corr).item())
    elapsed = time.time() - start_time
    elapsed = str(datetime.timedelta(seconds=elapsed))
                
    print("Execution Time:"+elapsed)
    print("Highest Correlation Achieved:" + str(round(max(corr_coeff),2)))
    snp_df = pd.DataFrame(snp_list) 
    snp_df.to_csv(os.path.join(curr_wd, args.output_path, args.row_matrices.split(".")[4] + '_full_SNP_list.csv'), index=False, header=['SNP_ID']) 
    max_ind = torch.argmax(torch.tensor(corr_coeff, device="cpu"))
    filtered_snp_df = pd.DataFrame(snp_list[0:max_ind])
    filtered_snp_df.to_csv(os.path.join(curr_wd, args.output_path, args.row_matrices.split(".")[4] + '_significant_SNP_list.csv'), index=False, header=['SNP_ID']) 
    print("Number of Effective SNPs:" + str(len(snp_list[0:max_ind])))
    
    plt.figure(figsize=(30,10))
    plt.plot(corr_coeff, marker='x')
    plt.title("Trend of Correlation Coefficient While Adding Features", fontsize=36)
    plt.xlabel("Index", fontsize=30)
    plt.tick_params(axis='both', labelsize=28)
    plt.ylabel("Correlation Coeeficient", fontsize=30)
    plt.savefig(os.path.join(curr_wd, args.save_plot_path, args.row_matrices.split(".")[4]+"_correlation_trend.png"))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="phenotype-to-genotype")
    parser.add_argument('--cuda', default= True)
    parser.add_argument('--gpu_device', default="7", help="GPU devices to be used")
    # input details
    parser.add_argument('--phenotype_file', type=str, default='https://raw.githubusercontent.com/genophenoenvo/JAGS-logistic-growth/main/data_clean/mac_growth_rate_modeled_season6.csv')
    parser.add_argument('--phenotype', type=str, default='max_growth_cm_gdd', choices=['max_growth_cm_gdd','max_height_cm'], help='PhenoCam Site ID')
    parser.add_argument('--row_matrices_path', type=str, default= 'row_matrices_v0.0.4',help = 'full path of row_matrices given in tar.gz format')
    parser.add_argument('--row_matrices', type=str, default= 'sorghum.filtered.season4.season6.max_growth_cm_gdd_p0005_qtl.row_matrices.tar.gz',help = 'full path of row_matrices given in tar.gz format')
    parser.add_argument('--num_iterations', type=int, default=None)
    # output path details
    parser.add_argument('--target_path', type=str, default='target_correlation_matrix')
    parser.add_argument('--output_path', type=str, default='SNP_outputs_v0.0.4')
    parser.add_argument('--save_plot_path', type=str, default='correlation_trend_plots_v0.0.4')
    args = parser.parse_args()
    main(args)
