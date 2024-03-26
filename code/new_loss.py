import os
import numpy as np


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

#know the pep-CDR residue-residue correlation & peptide-sequence & CDR residues 
#use peptide-sequence get the CDR-residue-probability to be the loss
#compare the loss map with the CDR residues


def old_pep_TCR_loss(pep_seq, cdr_seq):
    
    pep_cdr_freqmat = [[0.0, 0.0, 0.0, 0.04857142857142857, 0.017142857142857144, 0.0, 0.04, 0.08, 0.04285714285714286, 0.0, 0.0, 0.002857142857142857, 0.02, 0.07714285714285715, 0.0, 0.04857142857142857, 0.0, 0.0, 0.022857142857142857, 0.0], [0.16666666666666666, 0.05833333333333333, 0.0, 0.5583333333333333, 0.08333333333333333, 0.2916666666666667, 0.2333333333333333, 0.3916666666666666, 0.5083333333333333, 0.0, 0.0, 0.425, 0.08333333333333333, 0.08333333333333333, 0.17500000000000002, 0.05833333333333333, 0.008333333333333333, 0.0, 0.08333333333333333, 0.125], [0.3333333333333333, 0.022222222222222223, 0.0, 0.1, 0.0, 0.2, 0.0, 0.30000000000000004, 0.2666666666666667, 0.0, 0.0, 0.35555555555555557, 0.0, 0.0, 0.0, 0.0, 0.0, 0.03333333333333333, 0.0, 0.022222222222222223], [0.36363636363636365, 0.0, 0.0, 0.0, 0.0, 0.2090909090909091, 0.07272727272727272, 0.24545454545454548, 0.23636363636363636, 0.0, 0.0, 0.6363636363636365, 0.0, 0.0, 0.0, 0.018181818181818184, 0.0, 0.018181818181818184, 0.2, 0.0], [0.04285714285714286, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.09999999999999999, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.05333333333333334, 0.04, 0.18666666666666665, 0.23999999999999996, 0.31333333333333335, 0.0, 0.0, 0.21333333333333335, 0.0, 0.07333333333333333, 0.0, 0.14, 0.006666666666666667, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.00625, 0.0, 0.0625, 0.1, 0.0, 0.4312500000000001, 0.0, 0.0, 0.09375, 0.025, 0.00625, 0.0, 0.06875, 0.0, 0.049999999999999996, 0.0, 0.0], [0.0, 0.0, 0.0, 0.14032258064516132, 0.0, 0.01935483870967742, 0.01129032258064516, 0.06290322580645161, 0.1306451612903226, 0.0, 0.0, 0.02258064516129032, 0.004838709677419356, 0.0, 0.0032258064516129032, 0.0016129032258064516, 0.0, 0.016129032258064516, 0.0, 0.0], [0.1, 0.07692307692307693, 0.0, 0.03076923076923077, 0.05384615384615384, 0.03076923076923077, 0.0, 0.07692307692307693, 0.3692307692307692, 0.0, 0.0, 0.37692307692307697, 0.0, 0.023076923076923075, 0.0, 0.0, 0.03076923076923077, 0.10769230769230768, 0.07692307692307693, 0.04615384615384615], [0.0, 0.0, 0.0, 0.018181818181818184, 0.0, 0.2636363636363636, 0.045454545454545456, 0.175, 0.1954545454545455, 0.0, 0.0, 0.2590909090909091, 0.0659090909090909, 0.07954545454545454, 0.0, 0.24318181818181817, 0.10909090909090909, 0.05454545454545454, 0.0, 0.022727272727272728], [0.012264150943396227, 0.0, 0.0, 0.03490566037735849, 0.016037735849056607, 0.03207547169811321, 0.020754716981132078, 0.08490566037735849, 0.07641509433962264, 0.0, 0.0, 0.06792452830188679, 0.019811320754716984, 0.0018867924528301887, 0.009433962264150943, 0.020754716981132078, 0.008490566037735849, 0.009433962264150943, 0.002830188679245283, 0.009433962264150943], [0.0, 0.0, 0.0, 0.6599999999999999, 0.2, 0.5, 0.0, 0.42000000000000004, 0.06, 0.0, 0.0, 0.4, 0.27999999999999997, 0.16, 0.0, 0.12, 0.0, 0.0, 0.34, 0.22000000000000003], [0.2409090909090909, 0.0, 0.0, 0.09545454545454544, 0.0, 0.05909090909090909, 0.2545454545454546, 0.21818181818181823, 0.05, 0.0, 0.0, 0.45909090909090905, 0.022727272727272728, 0.0, 0.0, 0.031818181818181815, 0.0, 0.05454545454545454, 0.0, 0.08636363636363636], [0.12195121951219515, 0.0, 0.0, 0.058536585365853655, 0.0024390243902439024, 0.18780487804878046, 0.03170731707317073, 0.17804878048780487, 0.30000000000000004, 0.0, 0.0, 0.1048780487804878, 0.05365853658536586, 0.0024390243902439024, 0.06829268292682926, 0.01951219512195122, 0.04878048780487805, 0.046341463414634146, 0.06097560975609756, 0.08536585365853659], [0.08888888888888889, 0.0, 0.0, 0.10740740740740742, 0.003703703703703704, 0.22962962962962963, 0.07777777777777778, 0.2185185185185185, 0.12962962962962962, 0.0, 0.0, 0.2740740740740741, 0.025925925925925925, 0.0, 0.0, 0.1888888888888889, 0.0, 0.018518518518518517, 0.07777777777777778, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.007692307692307693, 0.0, 0.06923076923076923, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.007692307692307693, 0.06923076923076923, 0.0, 0.0], [0.10930232558139535, 0.0, 0.0069767441860465115, 0.20697674418604653, 0.046511627906976744, 0.0, 0.2558139534883721, 0.10465116279069768, 0.046511627906976744, 0.0, 0.0, 0.17906976744186048, 0.046511627906976744, 0.10232558139534882, 0.16511627906976745, 0.030232558139534883, 0.013953488372093023, 0.0, 0.023255813953488372, 0.013953488372093023], [0.3333333333333333, 0.03333333333333333, 0.0, 0.0, 0.2111111111111111, 0.011111111111111112, 0.13333333333333336, 0.0, 0.2777777777777778, 0.0, 0.0, 0.9222222222222223, 0.0, 0.1888888888888889, 0.1, 0.1111111111111111, 0.0, 0.0, 0.1, 0.1111111111111111], [0.12903225806451613, 0.0, 0.0, 0.08709677419354839, 0.1870967741935484, 0.3548387096774194, 0.14193548387096772, 0.10967741935483871, 0.22580645161290322, 0.0, 0.0, 0.061290322580645165, 0.25483870967741934, 0.01935483870967742, 0.0, 0.22580645161290322, 0.041935483870967745, 0.03225806451612903, 0.0, 0.16129032258064518], [0.014925373134328358, 0.0, 0.0, 0.0029850746268656717, 0.0, 0.035820895522388055, 0.034328358208955224, 0.013432835820895522, 0.13880597014925375, 0.0, 0.0, 0.12388059701492536, 0.0208955223880597, 0.02835820895522388, 0.07611940298507462, 0.029850746268656716, 0.0208955223880597, 0.0, 0.0, 0.013432835820895522]]
#    pep_cdr_discrete = [0.4000000000000001, 3.3333333333333335, 1.6333333333333335, 2.0, 0.14285714285714285, 1.2666666666666664, 0.8437500000000001, 0.4129032258064516, 1.4, 1.5318181818181817, 0.4273584905660377, 3.36, 1.5727272727272728, 1.3707317073170735, 1.4407407407407407, 0.15384615384615385, 1.3511627906976744, 2.5333333333333337, 2.032258064516129, 0.5537313432835821]
    dictionary = dict(zip(AA_freq.keys(), np.linspace(0,19,20,dtype = int)))
    
    
    pep_seq_num = []
    for i in pep_seq:
        pep_seq_num.append(dictionary[i])
        
#    if cdr_seq == None:
#        loss = 0
#        for i in pep_seq_num:
#            if pep_cdr_discrete[i] == 0:
#                loss += -np.log(1e-20)
#            else:
#                loss += - np.log(pep_cdr_discrete[i])
#
#    else:
    #cdr-str 2 num
    cdr_seq_num = []
    for i in cdr_seq:
        cdr_seq_num.append(dictionary[i])
    
    #make the score_map
    #calculate frequency
    cdr_freq = np.zeros(20)
    for j in pep_seq_num:
        cdr_freq += pep_cdr_freqmat[j]

    cdr_freq = cdr_freq / len(pep_seq)

    #calculate the loss
    loss = 0
    for i in cdr_seq_num:
        if cdr_freq[i] == 0:
            loss += -np.log(1e-20)
        else:
            loss += -np.log(cdr_freq[i])
    return loss



def IEDB_loss(iedb_freqmat, pep_seq):
    #iedb_freqmat line --> 20aa, column --> 9 position
    dictionary = dict(zip(AA_freq.keys(), np.linspace(0,19,20,dtype = int)))
    
    pep_seq_num = []
    for i in pep_seq:
        pep_seq_num.append(dictionary[i])
    
    iedb_prob = []
    for i in range(len(pep_seq_num)):
        iedb_prob.append(iedb_freqmat[pep_seq_num[i], i])
    
    loss = 0
    for i in iedb_prob:
        loss += -np.log(i)
    
    return loss

def pep_TCR_loss(pep_seq, cdr_seq, TCRpos):
    #tcr_freqmat line --> 20aa, column --> 9 position
    if cdr_seq == None:
        tcr_freqmat = np.array([[1.70045361e-01, 3.73620924e-01, 1.55740846e-01, 0.00000000e+00,
            0.00000000e+00, 3.00243131e-01, 2.24856037e-01, 5.77705508e-02,
            0.00000000e+00],
           [4.90819838e-01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
            3.33524892e-01, 0.00000000e+00, 5.40854722e-01, 0.00000000e+00,
            0.00000000e+00],
           [3.84963088e-02, 0.00000000e+00, 0.00000000e+00, 6.10911195e-02,
            2.61592467e-02, 0.00000000e+00, 1.27262034e-02, 0.00000000e+00,
            0.00000000e+00],
           [0.00000000e+00, 0.00000000e+00, 1.45648112e-01, 1.12161445e-01,
            0.00000000e+00, 0.00000000e+00, 2.33649240e-02, 0.00000000e+00,
            0.00000000e+00],
           [0.00000000e+00, 0.00000000e+00, 1.76929193e-05, 0.00000000e+00,
            0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
            1.91338751e-03],
           [0.00000000e+00, 0.00000000e+00, 4.93333531e-02, 2.84931611e-02,
            0.00000000e+00, 1.58511194e-02, 5.93555606e-03, 3.65994543e-02,
            0.00000000e+00],
           [0.00000000e+00, 0.00000000e+00, 3.64531426e-03, 1.40360115e-02,
            0.00000000e+00, 0.00000000e+00, 5.26304396e-03, 0.00000000e+00,
            0.00000000e+00],
           [4.04252577e-03, 6.83244033e-04, 1.28162125e-03, 5.92174477e-03,
            1.47915481e-03, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
            0.00000000e+00],
           [4.84851469e-02, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
            0.00000000e+00, 0.00000000e+00, 1.60283379e-02, 8.64788171e-02,
            0.00000000e+00],
           [7.93882597e-02, 3.13975087e-01, 0.00000000e+00, 5.03935864e-02,
            1.07892790e-02, 2.52311520e-01, 0.00000000e+00, 0.00000000e+00,
            3.93158314e-01],
           [3.10301687e-03, 2.34838731e-02, 6.15763554e-03, 0.00000000e+00,
            7.02860321e-04, 2.43506138e-03, 6.83868580e-04, 1.05420691e-03,
            1.96358927e-01],
           [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 5.31828295e-01,
            1.13864566e-01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
            0.00000000e+00],
           [1.71840240e-02, 2.64295223e-01, 0.00000000e+00, 7.27197870e-02,
            2.33539858e-02, 9.10236862e-02, 0.00000000e+00, 1.75141009e-02,
            0.00000000e+00],
           [0.00000000e+00, 0.00000000e+00, 1.71926839e-01, 0.00000000e+00,
            9.27703904e-02, 0.00000000e+00, 1.12829593e-01, 0.00000000e+00,
            0.00000000e+00],
           [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.11809795e-01,
            0.00000000e+00, 2.56580181e-01, 1.74687525e-02, 0.00000000e+00,
            0.00000000e+00],
           [1.28681032e-04, 0.00000000e+00, 2.20980294e-05, 1.70173638e-05,
            0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 4.91822936e-05,
            0.00000000e+00],
           [0.00000000e+00, 2.39416485e-02, 5.98792896e-02, 1.15280376e-02,
            2.22134005e-02, 5.77187870e-02, 3.60219667e-02, 1.99904757e-01,
            0.00000000e+00],
           [0.00000000e+00, 0.00000000e+00, 9.86667062e-02, 0.00000000e+00,
            1.95212564e-01, 0.00000000e+00, 0.00000000e+00, 1.46397817e-01,
            0.00000000e+00],
           [1.48306838e-01, 0.00000000e+00, 3.05619763e-01, 0.00000000e+00,
            1.76362143e-01, 0.00000000e+00, 0.00000000e+00, 4.53466706e-01,
            0.00000000e+00],
           [0.00000000e+00, 0.00000000e+00, 2.06072924e-03, 0.00000000e+00,
            3.56751758e-03, 2.38365137e-02, 3.96699537e-03, 7.64407455e-04,
            4.08569371e-01]])
        dictionary = dict(zip(AA_freq.keys(), np.linspace(0,19,20,dtype = int)))
        TCRpos = list(map(int, TCRpos.split('_')))
        
        pep_seq_num = []
        for i in pep_seq:
            pep_seq_num.append(dictionary[i])
        
        tcr_prob = []
        for i in range(len(pep_seq_num)):
            tcr_prob.append(tcr_freqmat[pep_seq_num[i], i])
        
        loss = 0
        n = 0
        for i in tcr_prob:
            n += 1
            if n not in TCRpos:
                loss += 0
            else:
                if i == 0:
                    loss += -np.log(1e-20)
                else:
                    loss += -np.log(i)
        
        return loss
    else:
        return old_pep_TCR_loss(pep_seq, cdr_seq)

 