import numpy as np
import pandas as pd
from operator import itemgetter
import math
from collections import Counter
from itertools import chain
from sklearn import preprocessing
# import operator

import itertools as irt
import cProfile
# print(int(math.factorial(9)/9))
# exit()
# index_check = [0, 0, 0, 2, 6, 24, 120, 720, 5040, 40320, 362880]
# print(index_check[9-1])
# exit()
# k = 1
# l = 0
# for i in range (0, 10, k+l):
#     print(i, l)
#     l = l + 2
# i = 0
# while i < 10:
#     print(i)
#     i = i +2
# # myiter = iter(range(0, 10))
# # for i in myiter:
# #     print(i)
# #     # next(myiter, None)
# #     i = i+2
# exit()

#==========================================================================================
#
#==========================================================================================
# reference samples cluster results
#==========================================================================================

# reference_dataFrame = pd.read_csv('reference_samples_80.txt', sep="\t", dtype=str, header=None)
reference_dataFrame = pd.read_csv('sample1.txt', sep="\t", dtype=str, header=None)
reference_Set1 = reference_dataFrame.drop(0, axis=0)  # drop first row
reference_gene_names = reference_Set1.values
reference_gene_names = reference_gene_names[:, 0]  # take only gene name column

total_gene_size = len(reference_gene_names)
print(total_gene_size, "total gene size")
# exit()

# dataFrame_all_samples_cls_res = pd.read_csv('exp_80.txt', sep="\t", header=None)
dataFrame_all_samples_cls_res = pd.read_csv('exp_reference.txt', sep="\t", header=None)
# exit()
dataFrame_all_samples_cls_res = dataFrame_all_samples_cls_res.T
all_samples = dataFrame_all_samples_cls_res.values
reference_samples = all_samples[:-1]  # DELETE LAST ROW OF PANDAS DATAFRAME
print(len(reference_samples))
print(reference_samples.shape)
# exit()

master_gene_class = [[], [], [], []]

for i_genes in range(len(reference_gene_names)):
    if (int(reference_samples[i_genes]) == 1):
        master_gene_class[0].append(reference_gene_names[i_genes])

    if (int(reference_samples[i_genes]) == 2):
        master_gene_class[1].append(reference_gene_names[i_genes])

    if (int(reference_samples[i_genes]) == 3):
        master_gene_class[2].append(reference_gene_names[i_genes])

    if (int(reference_samples[i_genes]) == 4):
        master_gene_class[3].append(reference_gene_names[i_genes])


#==========================================================================================
# random samples cluster results
#==========================================================================================

main_diff = []
for main_itr in range(1, 2):
    sample_name = 'exp_rand_sample_%d.txt' % main_itr
    print(sample_name)
    # exit()
    # dataFrame_cls_res = pd.read_csv('exp_rand_8.txt', sep="\t", header=None)
    dataFrame_cls_res = pd.read_csv(sample_name, sep="\t", header=None)
    total_iteration_time = len(dataFrame_cls_res)
    #for random samples
    dataFrame_cls_res = dataFrame_cls_res.T
    ll = dataFrame_cls_res.values
    random_samples_cluster_results = ll[:-1] # DELETE LAST ROW OF PANDAS DATAFRAME
    # print(len(random_samples_cluster_results), random_samples_cluster_results.shape)
    # exit()


    gene_names = [[], [], [], []]

    sub_diff = []

    for j_samples in range(total_iteration_time): # random samples repetition time

        gene_names = [[], [], [], []]

        # print(gene_names, sub_diff)

        for i_genes in range(total_gene_size):
                if (int(random_samples_cluster_results[i_genes,j_samples]) == 1):
                    gene_names[0].append(reference_gene_names[i_genes])

                if (int(random_samples_cluster_results[i_genes, j_samples]) == 2):
                    gene_names[1].append(reference_gene_names[i_genes])

                if (int(random_samples_cluster_results[i_genes,j_samples]) == 3):
                    gene_names[2].append(reference_gene_names[i_genes])

                if (int(random_samples_cluster_results[i_genes, j_samples]) == 4):
                    gene_names[3].append(reference_gene_names[i_genes])

        buffer_count = []

        all_permutation_result = list(irt.permutations([gene_names[0], gene_names[1], gene_names[2], gene_names[3]]))
        reference_cluster_result = list([master_gene_class[0], master_gene_class[1], master_gene_class[2], master_gene_class[3]])
        bu = []
        k = 0
        l_l = 0
        # print(len(all_permutation_result))
        min_distance_diff = 0
        pmut = 0
        wh_range = len(all_permutation_result)
        # print(wh_range)
        # exit()
        # for pmut in range(len(all_permutation_result)):
        while pmut < wh_range:
            print(pmut)

            if pmut == 0:
                single_permuted_res = np.array(all_permutation_result[pmut])
                k = 0
                # print(len(all_permutation_result[pmut]))
                # len(all_permutation_result[pmut]) - com_m = new_var; new_var!/new_var; new_var = new_var -1; pmut = pmut+ new_var
                for com_m in range(len(all_permutation_result[pmut])):
                    aa = single_permuted_res[com_m]
                    bb = reference_cluster_result[com_m]
                    # print(aa, bb)
                    diff_number = np.setdiff1d(bb, aa)
                    bu.append(diff_number.size)
                    k += diff_number.size
                    l_l = l_l + 1
                # print(l_l)
                min_distance_diff = k

            else:
                single_permuted_res = np.array(all_permutation_result[pmut])
                k = 0
                # print(len(all_permutation_result[pmut]))
                com_m = 0
                nes_wh = len(all_permutation_result[pmut])
                last_two_comb = nes_wh -2
                while com_m < nes_wh: #for com_m in range(len(all_permutation_result[pmut])):
                    aa = single_permuted_res[com_m]
                    bb = reference_cluster_result[com_m]
                    # print(aa, bb)
                    diff_number = np.setdiff1d(bb, aa)
                    bu.append(diff_number.size)
                    k += diff_number.size
                    l_l = l_l + 1

                    if (k > min_distance_diff) and (com_m < last_two_comb):
                        index_check = [0, 0, 0, 2, 6, 24, 120, 720, 5040, 40320, 362880]
                        # new_var! / new_var;
                        # new_var = new_var - 1;
                        index_var = len(all_permutation_result[pmut]) - com_m
                        new_var_1 = int(math.factorial(index_var) / index_var)
                        new_var = index_check[index_var]
                        pmut = pmut + (new_var -1)
                        print("test test ", new_var, new_var_1, index_var, com_m)
                        com_m = len(all_permutation_result[pmut])

                    com_m = com_m + 1
                # print(l_l)
                if min_distance_diff > k:
                    min_distance_diff = k
            pmut += 1   #================================================

        sub_diff.append(min_distance_diff)
        print("total ",l_l)
    # z = np.sum(X_normalized, axis=None)
    # print(sum(int(i) for i in X_normalized))
    print(sub_diff)
    norm = [(float(i) - min(sub_diff)) / (max(sub_diff) - min(sub_diff)) for i in sub_diff]
    # norm1 = [(float(i) - min(sub_diff)) / max(sub_diff) for i in sub_diff]
    # print(norm, norm1)
    norm1 = [1.15, 2, .99]
    main_diff.append(sum(i for i in norm))
print(main_diff)
f = open("myfile.txt", "w")
# mylist = [1, 2 ,6 ,56, 78]
f.write("\t".join(map(lambda x: str(x), main_diff)))
f.close()

exit()
