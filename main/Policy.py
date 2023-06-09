import numpy as npfrom arg_parse import get_argsargs = get_args()def cosine_similarity(a, b):    """    Calculates the cosine similarity between matrix 'a' and each matrix in 'B'.    Args:        a (array-like): Matrix 'a'.        b (array-like): Matrix 'b'.    Returns:        similarities (list): List of cosine similarity values.    """        dot_product = np.sum(a * b)    norm_a = np.linalg.norm(a)    norm_b = np.linalg.norm(b)    similarity = dot_product / (norm_a * norm_b)    return similaritydef calculate_similarity(new_peptide, databank_file):    """    Calculates the cosine similarity between 'new_peptide' and each sequence in 'databank_file'.    Args:        new_peptide (array-like): One-hot matrix representing the new peptide sequence.        databank_file (str): File name of the databank containing peptide sequences.    Returns:        similarities (list): List of cosine similarity values.    """    similarities = []    with open(databank_file, 'r') as f:        databank = f.read().splitlines()        for sequence in databank:            sequence_matrix = one_hot_encoding(sequence)            similarity = cosine_similarity(new_peptide, sequence_matrix)            similarities.append(similarity)    return similaritiesdef one_hot_encoding(sequence, position = args.TCR_loss_pos):    """    Converts a peptide sequence into a one-hot matrix.    Args:        sequence (str): Peptide sequence.        position (str): TCR_loss_positions, should be the same as the TCR_loss_pos    Returns:        matrix (numpy array): One-hot matrix representation of the sequence.    """    start = int(position.split('_')[0])    end = int(position.split('_')[-1])    amino_acids = "ARNDCQEGHILKMFPSTWYV"    num_positions = len(sequence)    num_amino_acids = len(amino_acids)    matrix = np.zeros((num_positions, num_amino_acids))    for i, char in enumerate(sequence):        if char in amino_acids:            j = amino_acids.index(char)            matrix[i, j] = 1    return matrix[start - 1 : end]