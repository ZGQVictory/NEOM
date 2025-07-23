# This is the python script about hallucinate the peptide on the MHC
# A very usual one, which is just a little try
# Two directions need to change. IEDB & PANDORA directions. Search '##change' 

import numpy as np
import os
### argparse part ###
from arg_parse import get_args
args = get_args(); print('#', args)
### Policy module ###
from Policy import calculate_similarity, one_hot_encoding
from PANDORA.Database import Database
### Main Process ###
from process import init_frequency_matrix, print_initial_frequency_matrix, Simulated_Annealing

### Additional part ###
from addition_func import get_freq
# Generate target directory
if not os.path.exists(args.o):
    os.makedirs(args.o)    


###### Global Parameters
AA_freq = {'A': 0.07421620506799341,
 'R': 0.05161448614128464,
 'N': 0.044645808512757915,
 'D': 0.05362600083855441,
 'C': 0.02468745716794485,
 'Q': 0.03425965059141602,
 'E': 0.0543119256845875,
 'G': 0.074146941452645,
 'H': 0.026212984805266227,
 'I': 0.06791736761895376,
 'L': 0.09890786849715096,
 'K': 0.05815568230307968,
 'M': 0.02499019757964311,
 'F': 0.04741845974228475,
 'P': 0.038538003320306206,
 'S': 0.05722902947649442,
 'T': 0.05089136455028703,
 'W': 0.013029956129972148,
 'Y': 0.03228151231375858,
 'V': 0.07291909820561925}

########### Changing Part ##############################
dir_tcr = '../data/pmhc/' ###change 1###
file_tcr = 'tcr-specific--peptide.txt'
# IEDB——HLA0201-Database
dir_iedb = '../data/pmhc/'  ###Change 2###
file_iedb = 'pHLA-A0201--peptide.txt'
# PANDORA Init
###Change 3 in the main function###
########################################################


def main():
    # Global parameters
    global AA_freq, dir_tcr, file_tcr, dir_iedb, file_iedb, db
    
    # Set the log file
    logfile_dir = os.path.join(args.o, '..')
    logfilename = os.path.join(logfile_dir, args.o.split('/')[-1] + '_Details_logfile.txt')
    # Set the Frequency matrix file
    freqmat_dir = logfile_dir
    freqmatfilename = os.path.join(freqmat_dir, args.o.split('/')[-1] + '_Frequencymap.txt')
    
    # Calculate the freqmat from the dataset
    iedb_freqmat, num_pep = get_freq(dir_iedb, file_iedb, args.l)
    prob_ref_freqmat, num_pep = get_freq(dir_tcr, file_tcr, args.l)
    
    # Initialize the frequency matrix
    AA_freq_mat = init_frequency_matrix(args, AA_freq, iedb_freqmat)
    
    # Print in the Details file
    print_initial_frequency_matrix(args, AA_freq_mat, AA_freq.keys(), logfilename)
    
    # Initialize the peptide object
    pep = peptide(length=args.l, aa_freq_mat=AA_freq_mat)
    # PANDORA Init 
    db = Database.load('../data/pmhc/pandora_Database.pkl') ###Change 3###
    
    ## Simulated Annealing part##
    Simulated_Annealing(args, pep, db, iedb_freqmat, prob_ref_freqmat, logfilename, freqmatfilename)
    
    
    
    
#######Background#########
class peptide():
    def __init__(self, length, aa_freq_mat):
        global AA_freq, dir_tcr, file_tcr
        self.len = length
        self.freq = aa_freq_mat
        self.init_sequences = {}
        self.dictionary = dict(zip(np.linspace(0,19,20,dtype = int), AA_freq.keys()))
        self.seq = self.generate_loc(loc = np.linspace(0,8,9,dtype = int)) #a dict
    
    def generate_loc(self, loc, seq = None, booster_num = args.booster_num):
        #loc should be the position, like [1,3,8]
        #input seq is dict
        if seq == None:
            seq = []
            if args.start_pep == None:
                
                for i in loc: #i is the peptide amino acid positions
                    seq.append(self.dictionary[np.random.choice(20, p = self.freq[:,i])])
            else:
                for i in args.start_pep:
                    seq.append(i)
                
        else:
            if args.add_booster:
                seq_init = list(seq.values())
                seq = list(seq.values())
                for cc in range(1000):
                    for i in loc: #i is the peptide amino acid positions
                        seq_i = seq[i]
                        while seq_i == seq_init[i] or seq_i == seq[i]:
                            if args.rest_mut_pos != None:
                                if i in [ _ - 1 for _ in list(map(int, args.rest_mut_pos.split('_')))]:
                                    seq_i = self.dictionary[np.random.choice([0,7,9,10,12,13,14,17,19])]
                                else:
                                    seq_i = self.dictionary[np.random.choice(20, p = self.freq[:,i])]
                            else:
                                seq_i = self.dictionary[np.random.choice(20, p = self.freq[:,i])]
                        seq[i] = seq_i
                        seq_peptide = one_hot_encoding(seq)  # One-hot matrix of the new peptide sequence
                        databank_file = os.path.join(dir_tcr, file_tcr)
                        target_similarities = max(calculate_similarity(seq_peptide, databank_file))
                    if target_similarities >= booster_num/len(args.TCR_loss_pos.split('_')):
                        break
            else:
                seq = list(seq.values())
                for i in loc: #i is the peptide amino acid positions
                    seq_i = seq[i]
                    while seq_i == seq[i]:
                        if args.rest_mut_pos != None:
                            if i in [ _ - 1 for _ in list(map(int, args.rest_mut_pos.split('_')))]:
                                seq_i = self.dictionary[np.random.choice([0,7,9,10,12,13,14,17,19])]
                            else:
                                seq_i = self.dictionary[np.random.choice(20, p = self.freq[:,i])]
                        else:
                            seq_i = self.dictionary[np.random.choice(20, p = self.freq[:,i])]
                    seq[i] = seq_i
        return dict(zip(np.linspace(0,8,9,dtype = int), seq)) #a dict

    def select_position(self, n_mutations):
        
        if args.fre_mut_pos == None:
            raise TypeError("fre_mut_pos should not be None")
        else:
            frepos = list(map(int, args.fre_mut_pos.split('_')))
        if args.rest_mut_pos == None:
            mutpos = frepos
        else:
            frepos.extend(list(map(int, args.rest_mut_pos.split('_'))))
            mutpos = frepos
        
        mut_loc = np.random.choice(list(set(mutpos)),n_mutations,replace = False) # mutation location  
        return [i - 1 for i in mut_loc]

    def mutate(self, mut_loc):
        return self.generate_loc(mut_loc, seq = self.seq) #a dictionary
    
    def find_key(self, val):
        return list(key for key, value in self.dictionary.items() if value == val)
    
    def change_freq(self, rate, loc, seq, accepted):
    #loc is mutation location
    #seq is a dictionary of peptide
        if accepted == True:
            for i in loc:
                key = self.find_key(seq[i])
                self.freq[key, loc] *= (1 + rate)
        else:
            for i in loc:
                key = self.find_key(seq[i])
                self.freq[key, loc] *= (1 - rate*0.25)   
        ###### IMPORTANT!! here is no Normalization!!!!!!!!###



if __name__ == '__main__':
    main()
