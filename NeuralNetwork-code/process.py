import numpy as np
import os
import sys
import csv

###argparse part####
from arg_parse import get_args
args = get_args()

###Loss module####
from Loss import compute_loss

###Additional part###
from addition_func import normalize


def load_matrix_from_csv(file_path):
    """
    Load a probability matrix from a CSV file.
    :param file_path: Path to the CSV file.
    :return: Loaded matrix as a NumPy array.
    """
    with open(file_path, 'r') as f:
        reader = csv.reader(f)
        matrix = np.array([list(map(float, row)) for row in reader])
    return matrix

def init_frequency_matrix(args, AA_freq, iedb_freqmat):

    if args.start_fre is None:
        # BLOSUM62
        sum_freq = np.sum(list(AA_freq.values()))
        adj_freq = [f/sum_freq for f in list(AA_freq.values())]
        AA_freq = dict(zip(AA_freq, adj_freq))
        AA_freq_mat = np.zeros((20, 9))
        for i in range(9):
            AA_freq_mat[:, i] = list(AA_freq.values())

        w_blosum, w_iedb = args.frequency_weight.split('_')
        w_blosum, w_iedb = float(w_blosum), float(w_iedb)
        AA_freq_mat = normalize(w_blosum * AA_freq_mat + w_iedb * iedb_freqmat, args.l)
    else:
        # You can change the frequency here.
        AA_freq_mat = load_matrix_from_csv("./start_freq.csv")
        AA_freq_mat = normalize(AA_freq_mat, args.l)

    # Return the initialized frequency matrix
    return AA_freq_mat

def print_initial_frequency_matrix(args, AA_freq_mat, keys, logfilename):
    
    logfile = open(logfilename, 'a')
    
    w_blosum, w_iedb = args.frequency_weight.split('_')
    w_blosum, w_iedb = float(w_blosum), float(w_iedb)
    print('We use, initially, the frequency matrix ({0}, {1}) of the amino acids is:\n'.format(20, args.l), file=logfile)
    print(AA_freq_mat, file=logfile)
    print('\n Note:\n', file=logfile)
    print('column --> location\n', file=logfile)
    print('row --> amino acids, {0}\n'.format(keys), file=logfile)
    if args.start_fre is None:
        print('It consists of two parts, one is BLOSUM62, whose weight is {}. The other is IEDB, whose weight is {}'.format(w_blosum, w_iedb), file=logfile)
    print('='*70, file=logfile)
    
    logfile.close()


def Simulated_Annealing(args, pep, iedb_freqmat, prob_ref_freqmat, logfilename, freqmatfilename, verbose = args.verbose):
    # Initialize the mutation rate
    Mi, Mf = args.mutation_rate.split('_')
    M = np.linspace(int(Mi), int(Mf), args.steps) # stepped linear decay of the mutation rate
    # Initialize the start loss score to infinity
    current_loss = np.inf
    # Some markers
    record = 0 #record the step corresponding to the acceptance
    t_seq = None #record the peptide sequence corresponding to the acceptance, and actually this is a dictionary.
    if args.nomemory is not True:
        Frequencymap = []# Use this for frequency memory
        Frequencymap.append(pep.freq)
    
    ##############################
    ### Simulated Annealing Start! 
    ##############################
    
    # Add the tolerance in optimizing #
    if args.tolerance is not None:
        rolling_window = []
        rolling_window_width = 100
    
    for i in range(args.steps): # i should be the step number!
        logfile = open(logfilename, 'a')
        if args.tolerance is not None and i > rolling_window_width: # check if change in loss falls under the tolerance threshold for terminating the simulation.
            if np.std(rolling_window[-rolling_window_width:]) < args.tolerance:
                print(f'The change in loss over the last 100 steps has fallen under the tolerance threshold ({args.tolerance}). Terminating the simulation...', file = logfile)
                sys.exit()
        else:
            T = args.T_init * (np.exp(np.log(0.5) / args.half_life) ** i)
            n_mutations = round(M[i])
            pos = None
            accepted = False # default
            
            # Treat the first step and the else differently
            if t_seq == None:  #As for the first step and we don't change the frequency matrix
                
                loss = compute_loss(iedb_freqmat, prob_ref_freqmat, pep.seq.values(), logfile=logfile)
                try_loss = loss
                delta = try_loss - current_loss
                assert delta < 0 # Of course the loss should lower than the infinity
                accepted = True
                print('-' * 100, file = logfile)
                print('Starting...', file = logfile)
                print('Loss {}:'.format(list(pep.seq.values())), try_loss, file = logfile)
                print(f'Step {i:05d}: change accepted >> LOSS {current_loss:2.3f} --> {try_loss:2.3f}', file = logfile)
                print('sequence --> target: ', list(pep.seq.values()), file = logfile)
                print('='*70, file = logfile)
                # Reset the current loss and extract the pdbfile
                current_loss = try_loss
                t_seq = pep.seq
                    
                assert record == 0 # Make sure that this is the first and only one time to execute
                
                
            else:
                pos = pep.select_position(n_mutations)
                print("########### Now Step-{} ###########".format(i), file = logfile)
                print("########### Mutation Positions are {} ###########".format(pos), file = logfile)
                loss_target = [] # the loss value, length is mutant times
                seq_target = [] # the sequence, length is mutant times
                metrorate = []
                
                # mutate the sequence
                for nmt in range(args.mutant_times):
                    # print("mutate", nmt)
                    nmt_seq = pep.mutate(pos)
                    # print("costtime_mutate")
                    seq_target.append(nmt_seq)
                    rate_element = 1
                    for _ in pos:
                        new_aa_num = pep.find_key(nmt_seq[_])
                        target_aa_num = pep.find_key(pep.seq[_])
                        if args.metropolis_hasting:
                            rate_element *= prob_ref_freqmat[target_aa_num, _]/pep.freq[new_aa_num, _]
                        else:
                            pass
                    metrorate.append(rate_element)
                    
                    if verbose:
                        print("\n", file = logfile)
                        print("########### Now Round{} ###########".format(nmt+1), file = logfile)
                        print('->', list(nmt_seq.values()), file = logfile)
                
                # Simultaneously construct the structure
                for r in range(args.mutant_times):
                    loss = compute_loss(iedb_freqmat, prob_ref_freqmat, list(seq_target[r].values()), logfile = logfile)
                    loss_target.append(loss)
                ##### Simulated Annealing main part is over
                
                # Change the frequency
                delta = [p - current_loss for p in loss_target]
                de_num = [] # delta decrease number, which means we must accept
                nde_num = [] # delta not decrease but accept number
                in_num = [] # not accept number
                print("Delta: ", delta, file = logfile)
                
                for i0 in range(len(delta)): # i0 is r above, indeed
                    delta0 = delta[i0]
                    rate0 = metrorate[i0]
                    if args.metropolis_hasting:
                        print("Rate: ", metrorate, file = logfile)
                        
                    if 1 < np.exp (-delta0 / T) * rate0:
                        de_num.append(i0)
                        # change frequency
                        if args.freqnotchange:
                            pass
                        else:
                            pep.change_freq(rate = args.frequence_change_rate, loc = pos, seq = seq_target[i0], accepted = True)
                        
                    else:
                        if np.random.uniform(0,1) < np.exp (-delta0 / T) * rate0:
                            nde_num.append(i0)
                            # change frequency
                            if args.freqnotchange:
                                pass
                            else:
                                pep.change_freq(rate = args.frequence_change_rate, loc = pos, seq = seq_target[i0], accepted=True)
                            
                        else:
                            in_num.append(i0)
                            # change frequency
                            if args.freqnotchange:
                                pass
                            else:
                                pep.change_freq(rate = args.frequence_change_rate, loc = pos, seq = seq_target[i0], accepted=False)
                pep.freq = normalize(pep.freq, args.l)
                # print("Freq changed")
                    
                ## Add the memory to the frequency matrix
                if args.nomemory is not True:
                    assert len(Frequencymap) >= 1
                    mem_freq = np.zeros((20,9))
                    for nF in range(len(Frequencymap)):
                        w_freq = np.exp(-nF/len(Frequencymap))
                        mem_freq = mem_freq + w_freq * Frequencymap[nF]
                        
                    pep.freq = pep.freq + 0.1*mem_freq/len(Frequencymap)		      
                    pep.freq = normalize(pep.freq, args.l)
                    Frequencymap.append(pep.freq)
                
                ## Save the frequency matrix to the frequency matrix file
                with open(freqmatfilename,'w') as freqmatfile:
                    print(pep.freq.tolist(), file = freqmatfile)
                    
                # Print the accepted sequence
                if len(de_num) != 0 or len(nde_num) != 0:
                    print("Accept {} times".format(len(de_num)+len(nde_num)), file = logfile)
                    
                    print("delta_decrease {} times:".format(len(de_num)), file = logfile)
                    for nn in de_num:
                        print(seq_target[nn].values(), file = logfile)
                    
                    print("delta_notdecrease {} times".format(len(nde_num)), file = logfile)
                    for nn in nde_num:
                        print(seq_target[nn].values(), file = logfile)
                    
                    print("Reject {} times".format(len(in_num)), file = logfile)
                    for nn in in_num:
                        print(seq_target[nn].values(), file = logfile)
                else:
                    print("All is rejected.", file = logfile)
                    for nn in in_num:
                        print(seq_target[nn].values(), file = logfile)
                        
                # Select the lowest loss and decide the acceptance
                try_loss = min(loss_target)
                index = loss_target.index(try_loss)
                
                if len(de_num) == 0 and len(nde_num) == 0:
                    accepted = False
                    
                    print(f'Step {i:05d}: change rejected >> LOSS {current_loss :2.3f} !-> {try_loss:2.3f}', file = logfile)
                    print('-'*70, file = logfile)
                
                else:
                    pep.seq = seq_target[index]
                    t_seq = pep.seq #mark the target sequence that is lower than before
                    record = i #mark the target step
                    if len(de_num) != 0:
                        accepted = True
                        print(f'Step {i:05d}: change accepted >> LOSS {current_loss:2.3f} --> {try_loss:2.3f}', file = logfile)
                        current_loss = try_loss
                        print('sequence --> target: ', list(pep.seq.values()), file = logfile)
                        print('='*70, file = logfile)
                    else:
                        accepted = True
                        print(f'Step {i:05d}: change accepted despite not improving the loss >> LOSS {current_loss:2.3f} --> {try_loss:2.3f}', file = logfile)
                        current_loss = try_loss
                        print('sequence --> target: ', list(pep.seq.values()), file = logfile)
                        print('='*70, file = logfile)
                        

                        
        print(pep.freq, file = logfile)
        print('Now the best target step is {}'.format(record), file = logfile)
        logfile.close()
        
        if args.tolerance is not None:
            rolling_window.append(current_loss)
            
    with open(logfilename, 'a') as f:
        print('#'*70, file = f)
