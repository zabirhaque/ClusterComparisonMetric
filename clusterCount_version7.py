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
# a=[]
# b=[]
# for i in range(1000000):
#    a.append(i)
#
# for i in range(0,1000000,2):
#    b.append(i)
# # for i in range(10):
#     # k = [c for c in b if c not in a]
# # diff_number = np.setdiff1d(b, a)
# cProfile.run('np.setdiff1d(b, a)')
# cProfile.run('pass')
# cProfile.run('for c in a:\n  if c not in b:\n    pass')
# cProfile.run('set(a) ^ set(b)')
# exit()
# # vv = list([a, b])
# # vv1 = np.array(vv)
# # print(vv1, vv1[0])
# a = [1, 2, 3]
# b = [4, 2]
# diff_number = np.setdiff1d(a, b)
# print(diff_number)
# # exit()
# a = [1, 2, 3]
# b = [4, 2]
# res1 = [9, 0]
# res2 = [8, 1, 9]
# print ([c for c in b if c not in a])
# print(set(a) ^ set(b))

# exit()

# x = list(irt.permutations([res1, res2]))
# # ref = list(irt([a, b]))
# df1 = pd.DataFrame({0: {0: a}, 1: {0: b}})  # data frame format ({  column: {row: value}, column: {row: value}    })
# df2 = pd.DataFrame({0: {0: res1}, 1: {0: res2}})
# df = pd.DataFrame(x)
# df1.columns = ['V2',	'V3']
# df2.columns = ['ref1', 'ref2']
# df3 = df1
# df3['tp1-tp2'] = list(df1.V3) - list(df2.ref2)
# print(df3)
# # df2 = df.values

# df.columns = ['V2',	'V3']
# print(df.values[:, 0] - df1.values[:, 0])
# v = df.groupby('Id').Day.transform('last') - df.Day
# df['Length'] = v.mask(v == 0)  # or v.mask(df.Status.eq('End'))

# exit()
# # y = np.array(x[1])
# # print(y, y[0])
# test = list([a, b])
# all_permutation_result = list(irt.permutations([res1, res2, res1]))
# reference_cluster_result = list([a, b, res2])
# bu = []
# k = 0
# min_distance_diff = 0
# for pmut in range(len(all_permutation_result)):
#     print(pmut)
#     single_permuted_res = np.array(all_permutation_result[pmut])
#     k = 0
#     for com_m in range(len(all_permutation_result[pmut])):
#
#         aa = single_permuted_res[com_m]
#         bb = reference_cluster_result[com_m]
#         print(aa, bb)
#         diff_number = np.setdiff1d(aa, bb)
#         bu.append(diff_number.size)
#         k = k + diff_number.size
#     # print(k, min_distance_diff)
#     if pmut == 0:
#         min_distance_diff = k
#
#     if min_distance_diff > k:
#         min_distance_diff = k
#
# print(bu, min_distance_diff)
# exit()
# buffer_count = []
# for main_matrix in range(len(master_gene_class)):
#     for sub_matrix in range(len(gene_names)):
#         s0 = np.array(master_gene_class[main_matrix])
#         s1 = np.array(gene_names[sub_matrix])
#         diff_number = np.setdiff1d(s0, s1)
#         buffer_count.append(diff_number.size)
#     min_count = min(int(i) for i in buffer_count)
#     sub_diff.append(min_count)
# print(sub_diff)
#
#
#
# exit()
# # c=[]
# # a=[1,2]
# # b=[3,4]
# # c.append([a,b])
# # print(len(c))
# arr= []
# arr.append([1, 2, 3])
# arr.append([4, 5, 6])
# arr1 = np.array(arr, np.float64)
#
#
# print(arr1.dtype)
# min_max_scaler = preprocessing.MinMaxScaler()
# X_train_minmax = min_max_scaler.fit_transform(arr1)
#
# X_normalized = preprocessing.normalize(arr1, norm='l2')
#
# print(X_train_minmax)
# print(X_normalized)
# z = np.sum(X_train_minmax, axis=None)
# print(z)
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

# print("mainsamples")
# print(len(master_gene_class[0]))
# print(len(master_gene_class[1]))
# print(len(master_gene_class[2]))
# print(len(master_gene_class[3]))
# exit()



#==========================================================================================
# random samples cluster results
#==========================================================================================

main_diff = []
for main_itr in range(1, 4):
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

        # print("subsamples")
        # print(len(gene_names[0]))
        # print(len(gene_names[1]))
        # print(len(gene_names[2]))
        # print(len(gene_names[3]))
        # s0 = np.array(gene_names[1])
        # s1 = np.array(gene_names[0])
        # diff_number = np.setdiff1d(s0, s1)
        # print(diff_number.size)
        # exit()
        buffer_count = []

        all_permutation_result = list(irt.permutations([gene_names[0], gene_names[1], gene_names[2], gene_names[3]]))
        reference_cluster_result = list([master_gene_class[0], master_gene_class[1], master_gene_class[2], master_gene_class[3]])
        bu = []
        k = 0
        l_l = 0
        # print(len(all_permutation_result))
        min_distance_diff = 0
        for pmut in range(len(all_permutation_result)):
            print(pmut)
            # pmut = pmut + 2
        # exit()

            single_permuted_res = np.array(all_permutation_result[pmut])
            k = 0
            print(len(all_permutation_result[pmut]))
            exit()

            for com_m in range(len(all_permutation_result[pmut])):

                aa = single_permuted_res[com_m]
                bb = reference_cluster_result[com_m]
                # print(aa, bb)
                diff_number = np.setdiff1d(bb, aa)
                bu.append(diff_number.size)
                k += diff_number.size
                l_l = l_l + 1
            # print(l_l)
            if pmut == 0:
                min_distance_diff = k

            if min_distance_diff > k:
                min_distance_diff = k
        # print(min_distance_diff)
        print(l_l)
        sub_diff.append(min_distance_diff)
    # arr1 = np.array(sub_diff, np.float64)

    # print(arr1, sub_diff)
    # min_max_scaler = preprocessing.MinMaxScaler()
    # X_train_minmax = min_max_scaler.fit_transform(arr1)

    # X_normalized = preprocessing.normalize(sub_diff, norm='l1')

    # print(X_train_minmax.reshape(-1, 1))
    # print(X_normalized)
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


    # for main_matrix in range(len(master_gene_class)):
    #     for sub_matrix in range(len(gene_names)):
    #         s0 = np.array(master_gene_class[main_matrix])
    #         s1 = np.array(gene_names[sub_matrix])
    #         diff_number = np.setdiff1d(s0, s1)
    #         buffer_count.append(diff_number.size)
    #     min_count = min(int(i) for i in buffer_count)
    #     sub_diff.append(min_count)
    # print(sub_diff)


    # arr.append([sub_diff])
    # main_diff.append(sum(int(i) for i in sub_diff))
    #
    # print(main_diff)
    #
    # if(jkl<3):
    #     jkl+=1
    # else:
    #     print(arr)
    #     print(len(arr))
    #     exit()

    # exit()

# arr = []
# arr.append([1,2,3])
# arr.append([4,5,6])
# exit()
print(len(main_diff))
exit()
total_main_diff = sum(int(i) for i in main_diff)
total_difference = total_main_diff/total_iteration_time
print(total_iteration_time, "times of repetition")
print(total_difference, "total diff")

exit()

# dataFrame = pd.read_csv('clusterCount.txt', sep="\t", header=None)
# dataFrame = pd.read_csv('random_v_20-30_40-50_60-70.txt', sep="\t", header=None)
# dataFrame.columns = ['gene_names', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9', 'V10', 'V11']

#FOR RANDOM SAMPLES
# dataFrame = pd.read_csv('random_v_10-20_30-40_50-60.txt', sep="\t", dtype=str, header=None)
# dataFrame.columns = ['gene_names', 'V20',	'V21',	'V22',	'V23',	'V24',	'V25',	'V26',	'V27',	'V28',	'V29',	'V30',	'V40',	'V41',	'V42',	'V43',	'V44',	'V45',	'V46',	'V47',	'V48',	'V49',	'V50',	'V60',	'V61',	'V62',	'V63',	'V64',	'V65',	'V66',	'V67',	'V68',	'V69',	'V70']
#
# #FOR V 72 TO V 79 SAMPLES
# # dataFrame.columns = ['gene_names', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9']
#
# Set1 = dataFrame.drop('gene_names', axis=1)  # DROP THE FIRST COLUMN AS IT IS NONE VALUE
# # df1 = dataFrame.dropna()
# trainData = Set1.drop(0, axis=0)
# # trainData = dataFrame.drop(0,axis=0)
# kk = trainData.values
#
