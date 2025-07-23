from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Wrapper import Wrapper

import os
import numpy as np
from arg_parse import get_args


def predict_structure(sequence, out_filename, db, mutant_times):
    
    ############################### Origin ########################
    ##sequence is a dictionary
    ##out_filename is the directory name you wanna generate
    ##return the template id then you can change the file name
    #p_seq = ''
    #for i in  sequence.values():
    #    p_seq = p_seq + i
    #target_id = out_filename
    #target = PMHC.Target(target_id,
    #                     allele_type = [args.mhc],
    #                     peptide = p_seq,
    #                     anchors = [2,9])    
    #case = Pandora.Pandora(target, db, output_dir = args.o)   #defalut output directions is '/home/zgq/data/pmhc/pandora/PANDORA/PANDORA_files/data'
    #case.model(n_loop_models = args.step_model_num)
    #return case.template.id
    ###############################################################
    args = get_args()
    mt = mutant_times
    if mt == 1:
        p_seq = ''
        for i in  sequence.values():
            p_seq = p_seq + i
        target_id = out_filename
        target = PMHC.Target(target_id,
                             allele_type = [args.mhc],
                             peptide = p_seq,
                             anchors = [2,9])    
        case = Pandora.Pandora(target, db, output_dir = args.o)   #defalut output directions is '/home/zgq/data/pmhc/pandora/PANDORA/PANDORA_files/data'
        case.model(n_loop_models = args.step_model_num)
    else:
        # here sequence, out_filename are lists. 
        #
        ## generate the tsv file ###
        tsvname = out_filename[0].split('_NUM')[0]
        
        
        os.makedirs(os.path.join(args.o, tsvname), exist_ok=True)
        newpath = os.path.join(args.o, tsvname)
        
        with open(os.path.join(args.o, tsvname+'.tsv'), 'w') as tsvfile:
            print('', file=tsvfile)
            
            for nmt in range(mt):
                p_seq = ''
                for i in  sequence[nmt].values():
                    p_seq = p_seq + i
                
                print(p_seq, end = '\t', file = tsvfile)
                print(args.mhc, end = '\t', file = tsvfile)
                print(out_filename[nmt], file = tsvfile)
        ###########################
        
        wrap = Wrapper.Wrapper()
        wrap.create_targets(os.path.join(args.o, tsvname+'.tsv'), db, 'I', IDs_col=2)
        wrap.run_pandora(num_cores = args.num_cores, output_dir = newpath, n_loop_models = args.step_model_num)
        
        
        while 1:
            new_num = []
            
            for a,b,c in os.walk(newpath):
                if a == newpath:
                    name = b
                    break
            if len(name) == mt:
                break
                
            else:
                filenum = [int(i.split('NUM')[-1]) for i in name]
                for i in np.arange(1,11):
                    if i not in filenum:
                        new_num.append(i)
                        
                tsvname += '1'
                with open(os.path.join(args.o, tsvname+'.tsv'), 'w') as tsvfile:
                    print('', file=tsvfile)
                    
                    for j in new_num:
                        nmt = j - 1 #you need to minus 1
                        p_seq = ''
                        for i in  sequence[nmt].values():
                            p_seq = p_seq + i
                        
                        print(p_seq, end = '\t', file = tsvfile)
                        print(args.mhc, end = '\t', file = tsvfile)
                        print(out_filename[nmt], file = tsvfile)
                        assert out_filename[nmt].split('NUM')[-1] not in filenum
                        
                wrap = Wrapper.Wrapper()
                wrap.create_targets(os.path.join(args.o, tsvname+'.tsv'), db, 'I', IDs_col=2)
                wrap.run_pandora(num_cores = args.num_cores, output_dir = newpath, n_loop_models = args.step_model_num)

