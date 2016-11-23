import pandas as pd
import numpy as np
import re


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

def rank_normalize(rank_dict, rank_matrix):
    '''
    Substitute rank values into a ranked matrix
    '''
    gene_names = rank_matrix.index
    norm_matrix = rank_matrix.apply(substitute_ranks_for_values, 
                                    rank_dict=rank_dict,
                                   gene_names=gene_names)
    
    return norm_matrix


def rank_normalization(dataframe):
    '''
    Perform rank normalization on the input matrix
    against an average overall all matrix columns
    '''
    
    # there is a bottleneck in here somewhere - find it!
    rank_vals = [rx for rx in range(1, dataframe.shape[0] + 1)]
    rank_matrix = dataframe.apply(rank_sort, ranks=rank_vals)
    sorted_matrix = dataframe.apply(sort_values)
    rank_dict = calculate_rank_values(sorted_matrix)
    
    # this was the slow part
    norm_matrix = rank_normalize(rank_dict, rank_matrix)
    
    return norm_matrix
