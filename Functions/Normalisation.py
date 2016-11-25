import pandas as pd
import numpy as np
import re
from statsmodels.robust import mad


def rank_sort(vector, ranks):
    '''
    Rank sort an input vector (pandas Series/numpy ndarray)
    '''
    
    order_vec = vector.sort_values(ascending=True)
    
    rank_series = pd.Series(dict(zip(order_vec.index, ranks)))
    
    return rank_series


def sort_values(vector):
    '''
    Sort the values of a vector, return the sorted vector
    from low to high
    '''
    
    order_vec = vector.sort_values(ascending=True)
    
    return order_vec.values


def calculate_rank_values(ranked_matrix):
    '''
    Calculate the average value for each rank,
    get new ranking and values attached
    '''
    
    av_values = ranked_matrix.median(axis=1)
    sorted_values = av_values.sort_values(ascending=True)
    ranks = [ri for ri in range(1, len(sorted_values) + 1)]
    
    rank_values = pd.Series(dict(zip(ranks, sorted_values.values)))
    return rank_values


def substitute_ranks_for_values(vector, rank_dict, gene_names):
    '''
    Substitute ranks for normalized values
    '''

    # pandas operations are FAST! make use of them
    # the rank_dict and vector share ranks, just convert
    # the vector to have the ranks as indices too then merge and munge
    rank_flip = pd.Series(dict(zip(vector, gene_names)))
    merge_rank = pd.concat([rank_flip, rank_dict], axis=1)
    merge_rank.columns = ['gene', vector.name]
    merge_rank.set_index('gene', drop=True, inplace=True)

    # the merge_rank is a dataframe, apply expects a series to be 
    # returned - stick some exceptions in here to catch any bugs
    return merge_rank.iloc[:, 0]


def bounded_normalize(rank_matrix):
    '''
    Normalize ranked values, bounded [0, 1]
    '''
    
    max_rank = np.max(rank_matrix)
    
    return rank_matrix/max_rank
    

def rank_normalize(rank_dict, rank_matrix):
    '''
    Substitute rank values into a ranked matrix
    '''
    gene_names = rank_matrix.index
    norm_matrix = rank_matrix.apply(substitute_ranks_for_values, 
                                    rank_dict=rank_dict,
                                   gene_names=gene_names)
    
    return norm_matrix


def rank_normalization(dataframe, bounded=False):
    '''
    Perform rank normalization on the input matrix
    against an average overall all matrix columns
    '''
    
    # there is a bottleneck in here somewhere - find it!
    rank_vals = [rx for rx in range(1, dataframe.shape[0] + 1)]
    rank_matrix = dataframe.apply(rank_sort, ranks=rank_vals)
    sorted_matrix = dataframe.apply(sort_values)
    rank_dict = calculate_rank_values(sorted_matrix)
    
    if bounded:
        norm_matrix = bounded_normalize(rank_matrix)
    else:
        norm_matrix = rank_normalize(rank_dict, rank_matrix)
    
    return norm_matrix


def pool_rank_normalise(dataframe, mapping, pool_rank_values):
    '''
    Rank normlase within a pool of samples against a reference
    using a reference quantile normalised value
    '''
    
    ranks = [ri for ri in range(1, len(pool_rank_values) + 1)]
    genes = dataframe.index
    _df_list = []
    for pool in mapping.keys():
        samples = mapping[pool]
        pool_df = dataframe.loc[:, samples]
        rank_matrix = pool_df.apply(rank_sort, ranks=ranks)
        rank_dict = pool_rank_values[pool].sort_values()
        rank_dict.index = ranks

        norm_matrix = rank_matrix.apply(substitute_ranks_for_values,
                                        rank_dict=rank_dict,
                                        gene_names=genes)
        _df_list.append(norm_matrix)
        
    normed_df = _df_list.pop(0)
    for _df in _df_list:
        df = pd.merge(normed_df, _df, left_index=True,
                      right_index=True, how='outer')
        normed_df = df
        
    return normed_df


def pool_normalise(dataframe, mapping, pool_rank_values,
                  cell_norm=False):
    '''
    Normalise sample-specific gene expression values against
    a pool-specific reference
    '''
    
    def div_pool(vector, denom):
        norm_vec = vector/denom
        norm_vec = norm_vec.fillna(0.0)
        norm_vec[np.isinf(norm_vec)] = 0

        return norm_vec
    
    def mul_pool(vector, factor):
        norm_vec = vector * factor
        norm_vec = norm_vec.fillna(0.0)
        norm_vec[np.isinf(norm_vec)] = 0

        return norm_vec
    
    def cell_normalisation(vector):
        '''
        Normalise within-cell to the maximal value
        '''
        
        max_val = np.max(vector)
        return vector/max_val
    
    _df_list = []
    for pool in mapping.keys():
        samples = mapping[pool]
        pool_df = dataframe.loc[:, samples]
        denom_array = pool_rank_values[pool]
        
        if cell_norm:
            _pool = pool_df.apply(cell_normalisation)
            normed_vals = _pool.apply(mul_pool, factor=denom_array)
        else:
            normed_vals = pool_df.apply(div_pool, denom=denom_array)
       
        _df_list.append(normed_vals)
        
    normed_df = _df_list.pop(0)
    for _df in _df_list:
        df = pd.merge(normed_df, _df, left_index=True,
                      right_index=True, how='outer')
        normed_df = df
        
    return normed_df


def pool_standard_normal(dataframe, mapping, pool_rank_values,
                        pool_rank_sd):
    '''
    Normalise within-pool values to standard normal using the
    pool reference mean and standard deviation
    '''
    
    def stand_norm(vector, mean, sd):
        s_norm = (vector - mean)/sd
        
        return s_norm
        
    _df_list = []
    for pool in mapping.keys():
        samples = mapping[pool]
        pool_df = dataframe.loc[:, samples]
        means = pool_rank_values[pool]
        sds = pool_rank_sd[pool]
        normed_vals = pool_df.apply(stand_norm, mean=means, sd=sds)
        _df_list.append(normed_vals)
    
    normed_df = _df_list.pop(0)
    for _df in _df_list:
        df = pd.merge(normed_df, _df, left_index=True,
                      right_index=True, how='outer')
        normed_df = df
        
    return normed_df


def standard_normal(dataframe, design_map, blocking_factors, parametric=False):
    '''
    Calculate Z-scores within pools of cells defined by multiple
    factors of a design matrix
    '''
    
    def stand_norm(vector, mean, sd):
        s_norm = (vector - mean)/sd
        
        return s_norm
        
    def mean_func(vector, parametric=False):
        if parametric:
            return np.mean(vector)
        else:
            return np.median(vector)
        
    def dev_func(vector, parametric=False):
        if parametric:
            return np.std(vector)
        else:
            return mad(vector, c=1)

    _df_list = []
        
    for pool, samp in design_map.groupby(blocking_factors):
        samples = samp.index
        pool_df = dataframe.loc[:, samples]
        means = pool_df.apply(mean_func, parametric=parametric, axis=1)
        sds = pool_df.apply(dev_func, parametric=parametric, axis=1)
        normed_vals = pool_df.apply(stand_norm, mean=means, sd=sds)
        _df_list.append(normed_vals)
    
    normed_df = _df_list.pop(0)
    for _df in _df_list:
        df = pd.merge(normed_df, _df, left_index=True,
                      right_index=True, how='outer')
        normed_df = df
    normed_df = normed_df.fillna(0.0)

    return normed_df
