import os
import numpy as np
import sys



def concludelog7(direction, fname, totalcon = True, acceptcon = True, pepcon = False, databank = None):
    loss = []
    step = []
    peptide = []
    mutant_times = 10
    aloss = []
    astep = []
    apeptide = []
    alosspack = []
    i = 0
    now_num, acc_num, suc_num = 0, 0, 0
    rnum = 10 #round number
    peppack = []
    apeppack = []
    apnum = []
    target_loss = []
    target_pep = []
    MARKER = False
    
    with open(os.path.join(direction, fname), 'r') as f:
        lines = f.readlines()
        
        for line in lines:
            i = i + 1 #line num starts from 1
            # Count how many rounds in a step
            #if "Round" in line:
                #if rnum < int(line.split("Round")[1].split()[0]):
                #    rnum = int(line.split("Round")[1].split()[0])
            #if rnum < mutant_times:
            #    rnum += 1
            if 'Now Step' in line: ########### Now Step-1 ###########
                ss = line.split("p-")[1].split()[0]
                step.append(ss)
                
            if '>> LOSS' in line and 'accept' in line:
                target_loss.append(float(line.split()[-1]))
                if int(line.split()[1].replace(':', '')) == 0: #step 0 situation
                    ss = 0
                    step.append(0)
                    aloss.append(float(line.split()[-1]))
                    loss.append(float(line.split()[-1]))
                else:
                    pass
            if 'sequence --> target' in line:
                target_pep0 = eval(line.split(' target: ')[-1])
                tar_pep = ''
                for _ in target_pep0:
                    tar_pep += _
                target_pep.append(tar_pep)
                astep.append(ss)
                
            if 'Loss [' in line and len(astep) == 0: #step 0 for the peptide
                pep = ''
                for element in eval(line.split(':')[0].split('s ')[1]):
                    pep += element
                peptide.append(pep)
                apeptide.append(pep)
            if 'RESULTS of the [' in line:
                    
                pep = ''
                for element in eval(line.split('RESULTS of the')[-1].split(':')[0]):
                    pep += element
                if len(peppack) < rnum:
                    peppack.append(pep)
                
                
            if 'Delta:' in line: #means peptide is collected over in one step
                peptide.append(peppack)
                peppack = []
                
                #calculate loss
                delta = eval(line.split("ta:")[1])
                loss_step = [_ + target_loss[-1] for _ in delta]
                loss.append(loss_step)
                
            if 'delta_decrease' in line:
                acc_num = int(line.split()[1])
                now_num = i
                suc_num = acc_num + i
                
            if 'delta_notdecrease' in line:
                acc_num = int(line.split()[1])
                if acc_num > 0 and i == suc_num + 1 and now_num < suc_num:
                    MARKER = True # means both have decrease and notdecrease
                now_num = i
                suc_num = acc_num + i

            if i>now_num and i <= suc_num:
                pep = ''
                if 'dict_values' in line:
                    line = line.replace('dict_values','')
                    line = line.replace(')','')
                    line = line.replace('(','')

                for element in eval(line):
                    pep += element
                if MARKER:
                    apeptide[-1].append(pep)
                    
                    for k in range(len(peptide[-1])):
                        if pep == peptide[-1][k]:
                            apnum.append(k)
                    for p in apnum:
                        aloss[-1].append(loss_step[p])
                        
                        
                    if i == suc_num:
                        MARKER = False
                        apnum = []
                        continue
                else:
                    apeppack.append(pep)
                
                if i == suc_num:
                    
                    apeptide.append(apeppack)
                    
                    for ap in apeppack:
                        for k in range(len(peptide[-1])):
                            if ap == peptide[-1][k]:
                                apnum.append(k)
                    for p in apnum:
                        alosspack.append(loss_step[p])
                    aloss.append(alosspack)
                    apnum = []
                    apeppack = []
                    alosspack = []

        




    #sl_file = open('{}_step_loss.log'.format(fname.split('.')[0]), 'w')
    #for i in step_loss.keys():
    #    print(i, step_loss[i], peptide[i], file = sl_file)
    #sl_file.close()
    #print(apeptide, aloss, astep)
    #print(len(apeptide), len(aloss), len(astep))
    if totalcon:
        with open('{}_step_loss.log'.format(fname.split('.')[0]), 'w') as f:
            nn = 0
            for i in range(len(step)):
                if i == 0:
                    print(step[i], loss[i], peptide[i], file = f)
                    nn += 1
                else:
                    
                    for j in range(len(peptide[i])):
                        print(step[i], loss[i][j], peptide[i][j], file = f)
                        nn += 1
            print("#Total {}".format(nn), file = f)
    
    #sl_file_accept = open('{}_step_loss_accept.log'.format(fname.split('.')[0]), 'w')
    #for i in astep_loss.keys():
    #    print(i, astep_loss[i], peptide[i], file = sl_file_accept)
    #sl_file_accept.close()
    if acceptcon:
        with open('{}_step_loss_accept.log'.format(fname.split('.')[0]), 'w') as f:
            nn = 0
            for i in range(len(astep)):
                if i == 0:
                    print(astep[i], aloss[i], apeptide[i], '*', file = f)
                    nn += 1
                else:
                    for j in range(len(apeptide[i])):
                        if apeptide[i][j] == target_pep[i]:
                            print(astep[i], aloss[i][j], apeptide[i][j], '*', file = f)
                        else:
                            print(astep[i], aloss[i][j], apeptide[i][j], file = f)
                        nn += 1
                        
            print("#Total {}".format(nn), file = f)
            
    if pepcon:
        with open(databank, 'a') as f:
            nn = 0
            print("# from {}".format(os.path.join(direction, fname)), file = f)
            for i in range(len(astep)):
                if i == 0:
                    print(apeptide[i], file = f)
                    nn += 1
                else:
                    for j in range(len(apeptide[i])):
                        if apeptide[i][j] == target_pep[i]:
                            print(apeptide[i][j], file = f)
                        else:
                            print(apeptide[i][j], file = f)
                        nn += 1
            print("# Totally, {} peptides.".format(nn), file = f)



    

    
    

def main():
    direction = sys.argv[1]
    fname = sys.argv[2]
    
    concludelog7(direction, fname)

    
if __name__ == '__main__':
    main()     

