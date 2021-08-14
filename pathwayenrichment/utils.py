import pickle
import random
import numpy as np


def saveToPickleFile(python_object, path_to_file='object.pkl'):
    """
    Save python object to pickle file
    """
    out_file = open(path_to_file, 'wb')
    pickle.dump(python_object, out_file)
    out_file.close()


def readFromPickleFile(path_to_file='object.pkl'):
    """
    Load python object from pickle file.
    Returns python object.
    """
    in_file = open(path_to_file, 'rb')
    python_object = pickle.load(in_file)
    return python_object


def computeSamplePvalue(distribution: list, statistic_value: float) -> float:
    """
    Compute sample p-value given null distribution and statistic value
    """
    N_equal_more_extreme = len(
        np.where(np.array(distribution) >= statistic_value)[0]
        )
    N_total = len(distribution)
    return N_equal_more_extreme / N_total


def copyProxyDictionary(proxydict) -> dict:
    """
    Convert ProxyDict, as returned by multiprocessing.Manager to regular dict.
    """
    return {k: v.copy() for k, v in proxydict.items()}


def randomPartition(elems: list, bin_sizes: list) -> list:
    """
    Randomly partition list elements into
    bins of given sizes
    """
    def shuffleList(elems):
        random.shuffle(elems)
        return elems

    elems = shuffleList(elems)
    partition = []
    start, end = 0, 0
    for bin_size in bin_sizes:
        end += bin_size
        partition.append(elems[start:end])
        start += bin_size
    return partition


def getCounts(array_like: list, sort_by_value=True) -> list:
    """
    Get dict of array item counts. Dict sorted
    by keys unless sort_by_value set to True.
    """
    unique_cond, counts_cond = np.unique(array_like, return_counts=True)
    counts = dict(zip(unique_cond, counts_cond))
    if sort_by_value:
        return dict(sorted(
            counts.items(), key=lambda item: item[1], reverse=True))
    else:
        return counts
