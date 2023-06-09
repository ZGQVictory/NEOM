import numpy as npimport osfrom arg_parse import get_argsargs = get_args()def normalize(freq_mat, length):    #freq_mat: (aa_type, length)    for i in range(length):        loc_sum = np.sum(freq_mat[:,i])        for j in range(20):            freq_mat[j,i] = freq_mat[j,i] / loc_sum    return freq_matdef get_freq(direction, filename, length):#    test0 = 0#    test = 0    aa= ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']    num_pep = 0    pos = np.zeros((20,length)) ##line --> 20aa, column --> position    dictionary = dict(zip(aa, np.linspace(0,19,20,dtype = int)))    with open(os.path.join(direction, filename), 'r') as f:        norm = True # justify if the sequence is normal (length is 9, all aa in the dictionary)        for line in f:#            test0 += 1            line = line.replace('\n','')            line = line.replace(' ','')#            if len(line) == 9:#               test += 1            if len(line) == length and line.isupper() == True:                for i in range(length):                    if line[i] in dictionary.keys():                        norm = True                        continue                    else:                        norm = False                        break                if norm == True:                    num_pep += 1                    for i in range(length):                        nrow = dictionary[line[i]]                        pos[nrow][i] += 1                 else:                    continue#    print(test0, test)    return normalize(pos,length), num_pepdef extract_pdbfile(out_filename, loss):    #direction = args.o + '/' + out_filename + '_' +structure    #direction = os.path.join(args.o, [i for i in os.listdir(args.o) if out_filename in i][0])    if '_NUM' not in out_filename:        direction = os.path.join(args.o, [i for i in os.listdir(args.o) if out_filename in i][0])    else:        tsvname = out_filename.split('_NUM')[0]        dir_base = os.path.join(args.o, tsvname)        direction = os.path.join(dir_base, [i for i in os.listdir(dir_base) if out_filename in i][0])        filename = os.path.join(direction,'molpdf_DOPE.tsv')    #f = pd.read_csv(filename, sep='\t', header=None)    with open(filename, 'r') as f:        for line in f:            ff = line.split()            pdbfile = ff[0]            break    #pdbfile = f.values[0][0]    pdbdir = os.path.join(args.o, 'LOSS')    if not os.path.exists(pdbdir):        os.makedirs(pdbdir)     pdbfilename = os.path.join(pdbdir, 'loss{}_{}'.format(loss, pdbfile))        oldpath = os.getcwd()    os.chdir(direction)    os.system('cp {} {}'.format(pdbfile, pdbfilename))    os.chdir(oldpath)