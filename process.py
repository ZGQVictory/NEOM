import numpy as np
import os
import sys

###argparse part####
from arg_parse import get_args
args = get_args()

## Structure_generation module ###
from Structure_generation import predict_structure

###Loss module####
from Loss import compute_loss

###Additional part###
from addition_func import normalize, extract_pdbfile


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
        AA_freq_mat = np.array([[0.05445672, 0.04283333, 0.05416201, 0.03832072, 0.0417253, 0.04806577, 0.04272074, 0.04776715, 0.05068996],
                                [0.06527764, 0.03411898, 0.05428123, 0.0410749, 0.14058102, 0.08669664, 0.0693538, 0.08691186, 0.0462773],
                                [0.04625828, 0.03358262, 0.04090161, 0.04353651, 0.0446419, 0.051113, 0.04256141, 0.05134544, 0.03637675],
                                [0.04599089, 0.04713368, 0.19022338, 0.3332639, 0.05574351, 0.06614387, 0.05945745, 0.05245689, 0.03801516],
                                [0.03531699, 0.03195011, 0.02984686, 0.02404189, 0.03113122, 0.03792322, 0.02945648, 0.03067051, 0.04168797],
                                [0.0515302, 0.03824897, 0.05156555, 0.03916856, 0.050257, 0.04462986, 0.0570856, 0.07839896, 0.03248476],
                                [0.03614986, 0.02537132, 0.04110387, 0.05984382, 0.04879752, 0.05414784, 0.06304654, 0.06933161, 0.04014434],
                                [0.04685121, 0.03549217, 0.03950241, 0.04812439, 0.04155226, 0.04356847, 0.04560972, 0.05807302, 0.04325099],
                                [0.04281785, 0.02183424, 0.03727409, 0.03318576, 0.0404686, 0.03882474, 0.05987389, 0.0456511, 0.02529345],
                                [0.04978223, 0.03750688, 0.04397547, 0.0270582, 0.05059307, 0.06178099, 0.04955908, 0.03841228, 0.04965674],
                                [0.04809141, 0.33157502, 0.03244782, 0.02298641, 0.04984459, 0.0674829, 0.04333029, 0.04542962, 0.07636724],
                                [0.06973224, 0.03377837, 0.02443732, 0.0269464, 0.06302107, 0.05093015, 0.03871595, 0.06316715, 0.03969051],
                                [0.0381406, 0.0561661, 0.05278914, 0.0239436, 0.02772596, 0.03726892, 0.04117237, 0.0298189, 0.046453],
                                [0.06472726, 0.02823674, 0.04640838, 0.03235684, 0.04933927, 0.04586214, 0.05585968, 0.03697165, 0.03792922],
                                [0.03758478, 0.03668568, 0.03267762, 0.03580889, 0.03618991, 0.04823634, 0.05752547, 0.03705552, 0.02916212],
                                [0.05931315, 0.03373177, 0.05622882, 0.03964355, 0.04996495, 0.03748624, 0.04934019, 0.04164926, 0.0375242],
                                [0.05816498, 0.03257287, 0.03627517, 0.04545204, 0.05050443, 0.04501964, 0.05077999, 0.06708078, 0.04058439],
                                [0.03020014, 0.02609972, 0.03483754, 0.02594058, 0.0347238, 0.0389703, 0.02980376, 0.03482518, 0.03728486],
                                [0.06071727, 0.02807822, 0.05551976, 0.03224669, 0.03984764, 0.0367096, 0.05825986, 0.05749681, 0.04108073],
                                [0.05889629, 0.04500319, 0.04554194, 0.02705633, 0.053347, 0.05913938, 0.05648771, 0.02748632, 0.21004632]])

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


def Simulated_Annealing(args, pep, db, iedb_freqmat, prob_ref_freqmat, logfilename, freqmatfilename, verbose = args.verbose):
    # Initialize the mutation rate
    Mi, Mf = args.mutation_rate.split('-')
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
                #PANDORA part
                out_filename = 'TestCase_Step' + str(i)
                try:
                    predict_structure(pep.seq, out_filename, db = db, mutant_times = 1)
                except:
                    sys.exit('Step {} collapsed!!!'.format(i))
                    continue
                loss = compute_loss(iedb_freqmat, prob_ref_freqmat, out_filename, pep.seq.values(), logfile=logfile)
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
                if args.extract_pdbfile:
                    extract_pdbfile(out_filename, try_loss)
                assert record == 0 # Make sure that this is the first and only one time to execute
                
                # Remove the PANDORA output
                if args.total_structures is not True:
                    direction = os.path.join(args.o, [i for i in os.listdir(args.o) if out_filename in i][0])
                    os.system('rm -rf {}'.format(direction))
                        
                
            else:
                pos = pep.select_position(n_mutations)
                print("########### Now Step-{} ###########".format(i), file = logfile)
                print("########### Mutation Positions are {} ###########".format(pos), file = logfile)
                loss_target = [] # the loss value, length is mutant times
                seq_target = [] # the sequence, length is mutant times
                out_filename = [] # the out_filename, length is mutant times
                metrorate = []
                
                # mutate the sequence
                for nmt in range(args.mutant_times):
                    nmt_seq = pep.mutate(pos)
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
                    out_filename.append('TestCase_Step' + str(i) + '_NUM' + str(nmt+1))
                    
                    if verbose:
                        print("\n", file = logfile)
                        print("########### Now Round{} ###########".format(nmt+1), file = logfile)
                        print('->', list(nmt_seq.values()), file = logfile)
                
                # Simultaneously construct the structure
                try:
                    predict_structure(seq_target, out_filename, db = db, mutant_times=args.mutant_times)
                except:
                    sys.exit('Step {} collapsed!!!'.format(i))
                    continue
                
                for r in range(args.mutant_times):
                    loss = compute_loss(iedb_freqmat, prob_ref_freqmat, out_filename[r], list(seq_target[r].values()), logfile = logfile)
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
                        
                        if args.extract_pdbfile:
                            extract_pdbfile(out_filename[i0], loss_target[i0])
                                                
                    else:
                        if np.random.uniform(0,1) < np.exp (-delta0 / T) * rate0:
                            nde_num.append(i0)
                            # change frequency
                            if args.freqnotchange:
                                pass
                            else:
                                pep.change_freq(rate = args.frequence_change_rate, loc = pos, seq = seq_target[i0], accepted=True)
                            
                            if args.extract_pdbfile:
                                extract_pdbfile(out_filename[i0], loss_target[i0])
                        else:
                            in_num.append(i0)
                            # change frequency
                            if args.freqnotchange:
                                pass
                            else:
                                pep.change_freq(rate = args.frequence_change_rate, loc = pos, seq = seq_target[i0], accepted=False)
                pep.freq = normalize(pep.freq, args.l)
                
                    
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
                        
                # Remove the PANDORA output
                if args.total_structures is not True:
                    tsvname = out_filename[-1].split('_NUM')[0]
                    dir_base = os.path.join(args.o, tsvname)
                    os.system('rm -rf {}'.format(dir_base))
                        
        print(pep.freq, file = logfile)
        print('Now the best target step is {}'.format(record), file = logfile)
        logfile.close()
        
        if args.tolerance is not None:
            rolling_window.append(current_loss)
            
    with open(logfilename, 'a') as f:
        print('#'*70, file = f)
