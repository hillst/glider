def calc_information(motif):
    import numpy as np 
    """
    Shannon entropy :)
    """
    return -1 * np.multiply(motif, np.log2(motif, where = motif != 0))

def error_correction(info, epsilon=0.001):
    return 2 - (info + epsilon)

def get_scores(motif):
    info = calc_information(motif)
    scores = motif * error_correction(info)
    return scores 


