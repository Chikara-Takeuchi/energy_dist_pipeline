import pandas as pd
import numpy as np
from sklearn.metrics import pairwise_distances
import multiprocessing as mp
from itertools import combinations

from tqdm import tqdm
from tqdm.contrib.concurrent import process_map

import torch
import gc

import scanpy as sc
import pickle
import os
import torch

def pairwise_torch(X,cell_id_list,device,batch_size=100,norm_p=2.0,vardose=False,failsafe=True,downsampling=5000):
    cell_num_len = np.array([len(x) for x in cell_id_list])
    max_cell_num = max(cell_num_len)
    res = []

    test_tensor_nonpad = []
    test_tensor = []
    musk_tensor = []
    for index,cell_names in enumerate(cell_id_list):
        pad_num = max_cell_num - cell_num_len[index]
        test_tensor += [torch.tensor(np.pad(X.loc[cell_names,:],((pad_num,0),(0,0))))]
        musk_tensor += [torch.tensor([[0]*pad_num+[1]*cell_num_len[index]])]
    test_tensor = torch.stack(test_tensor).to(device)
    musk_tensor = torch.stack(musk_tensor).to(device)

    res = []
    cell_num_len_tensor = torch.tensor(cell_num_len).to(device)
    
    combis = np.array(list(combinations(range(len(cell_id_list)), 2)) 
                      + [(x,x) for x in range(len(cell_id_list))])
    total_len_combis = np.arange(len(combis))
    calc_list = np.array_split(total_len_combis,len(combis) // batch_size)
    
    if vardose:
        pbar = tqdm(calc_list)
    else:
        pbar = calc_list
    
    for combis_idx in pbar:
        
        batch_sample_size = len(combis_idx)
        tensor_0 = test_tensor[combis[combis_idx][:,0]]
        tensor_1 = test_tensor[combis[combis_idx][:,1]]

        musk_2D = torch.mul(
            musk_tensor[combis[combis_idx][:,0]].view(batch_sample_size,-1,1),
            musk_tensor[combis[combis_idx][:,1]]
        )
        sum_cross = torch.cdist(tensor_0,tensor_1,p=norm_p).pow(2).mul(musk_2D).sum((1,2))
        sum_cross = sum_cross.mul((cell_num_len_tensor[combis[combis_idx][:,0]] 
                                   * cell_num_len_tensor[combis[combis_idx][:,1]]).reciprocal())
        res += [sum_cross]
    combi_1 = combis[:,0]
    combi_2 = combis[:,1]
    
    return_res = torch.concat(res).to("cpu").tolist()
    
    del test_tensor,musk_tensor,cell_num_len_tensor,tensor_0,tensor_1,musk_2D,sum_cross
    gc.collect()
    torch.cuda.empty_cache()
    
    return list(zip(combi_1,combi_2,return_res))


def permutation_test(X, test_cell1, test_cell2, device, batch_num=10, total_permute=1000,
                     norm_p=2.0, return_permute=True):
    test_cells_concat = np.concatenate((test_cell1, test_cell2)).tolist()
    repeat_num = int(total_permute/batch_num)
    num_cell1 = len(test_cell1)
    num_cell2 = len(test_cell2)
    num_all = num_cell1 + num_cell2
    
    # Initialization
    all_cells_dist = None
    pairwise_dist = None
    pairwise_dist_tmp = None
    
    try:
        all_cells_dist = torch.Tensor(X.loc[test_cells_concat,:].to_numpy()).to(device)
        pairwise_dist_tmp = []

        for index_arr in np.array_split(np.arange(num_all), 1):
            pairwise_dist_tmp.append(torch.cdist(all_cells_dist[index_arr], all_cells_dist, p=norm_p).pow(2))

        pairwise_dist = torch.cat(pairwise_dist_tmp).repeat(batch_num, 1, 1)

        sum_target = pairwise_dist[0][num_cell1:][:,num_cell1:].sum() / (num_cell2*num_cell2)
        sum_non_target = pairwise_dist[0][:num_cell1][:,:num_cell1].sum() / (num_cell1*num_cell1)
        sum_cross = pairwise_dist[0][:num_cell1][:,num_cell1:].sum() / (num_cell1*num_cell2)
        obs_edist = 2*sum_cross - sum_target - sum_non_target
            
    except Exception as e:
        _cleanup_tensors_on_error(locals(), device)
        raise e

    if not return_permute:
        obs_edist = obs_edist.to("cpu")
        _cleanup_tensors_success(pairwise_dist_tmp, pairwise_dist, sum_target, 
                               sum_non_target, sum_cross, all_cells_dist, device)
        return obs_edist
    
    # If return_permute=True
    e_dist_list = []
    np.random.seed(0)
    
    group_target = None
    group_non_target = None
    extracted = None
    
    try:
        for i in range(repeat_num):
            random_id = np.array([np.random.permutation(num_all) for i in range(batch_num)])

            group_target = torch.tensor(random_id[:,:num_cell1]).to(device)
            group_non_target = torch.tensor(random_id[:,num_cell1:]).to(device)

            extracted = pairwise_dist.gather(1, group_target.unsqueeze(2).expand(-1, -1, num_all))
            extracted = extracted.gather(2, group_target.unsqueeze(1).expand(-1, num_cell1, -1))
            sum_target = extracted.sum((1,2)) / (num_cell1*num_cell1)

            extracted = pairwise_dist.gather(1, group_non_target.unsqueeze(2).expand(-1, -1, num_all))
            extracted = extracted.gather(2, group_non_target.unsqueeze(1).expand(-1, num_cell2, -1))
            sum_non_target = extracted.sum((1,2)) / (num_cell2*num_cell2)

            extracted = pairwise_dist.gather(1, group_target.unsqueeze(2).expand(-1, -1, num_all))
            extracted = extracted.gather(2, group_non_target.unsqueeze(1).expand(-1, num_cell1, -1))
            sum_cross = extracted.sum((1,2)) / (num_cell2*num_cell1)

            e_dist_list.append(2*sum_cross - sum_target - sum_non_target)

        e_dist_list = torch.cat(e_dist_list)
        
    except Exception as e:
        _cleanup_permutation_error(locals(), device)
        raise e
    
    e_dist_list = e_dist_list.to("cpu")
    obs_edist = obs_edist.to("cpu")
    
    _cleanup_permutation_success(pairwise_dist, pairwise_dist_tmp, sum_target, 
                                sum_non_target, sum_cross, all_cells_dist, 
                                group_target, group_non_target, extracted, device)

    return (obs_edist, e_dist_list)


def _cleanup_tensors_on_error(local_vars, device):
    tensor_names = ['pairwise_dist', 'all_cells_dist', 'pairwise_dist_tmp']
    
    for name in tensor_names:
        if name in local_vars and local_vars[name] is not None:
            del local_vars[name]
    
    gc.collect()
    if 'cuda' in str(device) and torch.cuda.is_available():
        torch.cuda.empty_cache()


def _cleanup_tensors_success(pairwise_dist_tmp, pairwise_dist, sum_target, 
                           sum_non_target, sum_cross, all_cells_dist, device):
    del pairwise_dist_tmp, pairwise_dist, sum_target, sum_non_target, sum_cross, all_cells_dist


def _cleanup_permutation_error(local_vars, device):
    tensor_names = ['pairwise_dist', 'all_cells_dist', 'e_dist_list', 
                    'group_target', 'group_non_target', 'extracted']
    
    for name in tensor_names:
        if name in local_vars and local_vars[name] is not None:
            del local_vars[name]
    
    gc.collect()
    if 'cuda' in str(device) and torch.cuda.is_available():
        torch.cuda.empty_cache()


def _cleanup_permutation_success(pairwise_dist, pairwise_dist_tmp, sum_target, 
                                sum_non_target, sum_cross, all_cells_dist, 
                                group_target, group_non_target, extracted, device):
    del (pairwise_dist, pairwise_dist_tmp, sum_target, sum_non_target, 
         sum_cross, all_cells_dist, group_target, group_non_target, extracted)



def disco_test(X, test_cell_list, device, batch_num=5, total_permute=1000, norm_p=2.0):
    # Initialization
    all_cells_dist = None
    pairwise_dist = None
    pairwise_dist_tmp = None
    group_cell_id = None
    extracted = None
    
    try:
        alpha = 1
        test_cells_concat = np.concatenate(test_cell_list).tolist()
        num_of_samples = len(test_cell_list)
    
        length_sample_array = [len(i) for i in test_cell_list]
        sample_index_array = [[i] * length_sample_array[i] for i in range(num_of_samples)]
        sample_index_array = np.concatenate(sample_index_array)
        num_all = len(test_cells_concat)
    
        repeat_num = int(total_permute / batch_num)
        all_cells_dist = torch.Tensor(X.loc[test_cells_concat, :].to_numpy()).to(device)
        pairwise_dist_tmp = []
    
        for index_arr in np.array_split(np.arange(num_all), 1):
            pairwise_dist_tmp.append(torch.cdist(all_cells_dist[index_arr], all_cells_dist, p=norm_p).pow(alpha))

        pairwise_dist = torch.cat(pairwise_dist_tmp).repeat(batch_num, 1, 1)
        total_disp = (num_all / 2) * (pairwise_dist[0].sum() / num_all / num_all)
    
        within_disp_list = []
        for i in range(num_of_samples):
            cell_index = (sample_index_array == i)
            tmp = (length_sample_array[i] / 2) * pairwise_dist[0][cell_index, :][:, cell_index].sum() / length_sample_array[i] / length_sample_array[i]
            within_disp_list.append(tmp)
        within_disp = torch.stack(within_disp_list).sum()
    
        between_disp = total_disp - within_disp
        F_value = (between_disp / (num_of_samples - 1)) / (within_disp / (num_all - num_of_samples))

        np.random.seed(0)
        F_value_permute_list = []
    
        for i in range(repeat_num):
            random_id = np.array([np.random.permutation(num_all) for i in range(batch_num)])
    
            group_cell_id = []
            for k in range(num_of_samples):
                cell_index = (sample_index_array == k)
                group_cell_id.append(torch.tensor(random_id[:, cell_index]).to(device))
    
            within_disp_list_permute = []
            for k in range(num_of_samples):
                extracted = pairwise_dist.gather(1, group_cell_id[k].unsqueeze(2).expand(-1, -1, num_all))
                extracted = extracted.gather(2, group_cell_id[k].unsqueeze(1).expand(-1, length_sample_array[k], -1))
                extracted = (length_sample_array[k] / 2) * extracted
                within_disp_list_permute.append(extracted.sum((1, 2)) / (length_sample_array[k] * length_sample_array[k]))
            within_disp_permute = torch.sum(torch.stack(within_disp_list_permute), 0)
            total_disp_permute = total_disp.expand(batch_num)
            between_disp_permute = total_disp_permute - within_disp_permute
    
            F_value_permute_list.append((between_disp_permute / (num_of_samples - 1)) / (within_disp_permute / (num_all - num_of_samples)))
    
        F_value_permute_list = torch.cat(F_value_permute_list).to("cpu")
        F_value = F_value.to("cpu")
        
        _cleanup_disco_success(all_cells_dist, between_disp_permute, within_disp_permute, 
                              total_disp_permute, pairwise_dist, device)
    
        return (F_value, F_value_permute_list)

    except Exception as e:
        print(f"An error occurred in disco_test: {e}")
        print("Cleaning up GPU memory...")
        _cleanup_disco_error(locals(), device)
        raise e


def _cleanup_disco_success(all_cells_dist, between_disp_permute, within_disp_permute, 
                          total_disp_permute, pairwise_dist, device):
    del all_cells_dist, between_disp_permute, within_disp_permute, total_disp_permute, pairwise_dist


def _cleanup_disco_error(local_vars, device):
    tensor_names = ['pairwise_dist', 'all_cells_dist', 'pairwise_dist_tmp', 
                    'group_cell_id', 'extracted']
    
    for name in tensor_names:
        if name in local_vars and local_vars[name] is not None:
            del local_vars[name]
            
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()