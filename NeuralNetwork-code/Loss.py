import numpy as np
import os
from arg_parse import get_args
from addition_func import normalize

args = get_args()

AA_freq = {
    "A": 0.07421620506799341,
    "R": 0.05161448614128464,
    "N": 0.044645808512757915,
    "D": 0.05362600083855441,
    "C": 0.02468745716794485,
    "Q": 0.03425965059141602,
    "E": 0.0543119256845875,
    "G": 0.074146941452645,
    "H": 0.026212984805266227,
    "I": 0.06791736761895376,
    "L": 0.09890786849715096,
    "K": 0.05815568230307968,
    "M": 0.02499019757964311,
    "F": 0.04741845974228475,
    "P": 0.038538003320306206,
    "S": 0.05722902947649442,
    "T": 0.05089136455028703,
    "W": 0.013029956129972148,
    "Y": 0.03228151231375858,
    "V": 0.07291909820561925,
}


################################################################
# Write your New Loss function from the Neural Network

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils import data

import numpy as np
import os
import collections
from datetime import datetime
import time

from sklearn.model_selection import KFold


def try_gpu(i=0):
    if torch.cuda.device_count() >= i + 1:
        return torch.device(f'cuda:{i}')
    return torch.device('cpu')

class Vocab:
    def __init__(self, tokens=None, min_freq=0, reserved_tokens=None):
        if tokens is None:
            tokens = []
        if reserved_tokens is None:
            reserved_tokens = []
        
        counter = counter_corpus(tokens)
        self._token_freqs = sorted(counter.items(), key=lambda x: x[1], reverse=True)
        self.idx_to_token = ['<unk>'] + reserved_tokens
        self.token_to_idx = {token: idx for idx, token in enumerate(self.idx_to_token)}
        
        for token, freq in self._token_freqs:
            if freq < min_freq:
                break
            if token not in self.token_to_idx:
                self.idx_to_token.append(token)
                self.token_to_idx[token] = len(self.idx_to_token) - 1
    
    def __len__(self):
        return len(self.idx_to_token)
    
    def __getitem__(self, tokens):
        if not isinstance(tokens, (list, tuple)):
            return self.token_to_idx.get(tokens, self.unk)
        return [self.__getitem__(token) for token in tokens]
    
    def to_tokens(self, indices):
        if not isinstance(indices, (list, tuple)):
            return self.idx_to_token[indices]
        return [self.idx_to_token[index] for index in indices]
    
    @property
    def unk(self):
        return 0
    
    @property
    def token_freqs(self):
        return self._token_freqs
        
        
def counter_corpus(tokens):
    if len(tokens) == 0 or isinstance(tokens[0], list):
        tokens = [token for line in tokens for token in line]
    
    return collections.Counter(tokens)

def read_sequence(directory, file):
    with open(os.path.join(directory, file), 'r') as f:
        return f.readlines()
    
def tokenize(lines):
    alltoken = []
    # 读取这个文件，并将数字和seqeunce内容分开
    for line in lines:
        if list(line)[0] == '#':
            continue
        else:
            alltoken.append(line.strip().split())
    return alltoken

def load_seq(batch_size, is_train=False, directory='../data', file='All_loss_databank.log', max_tokens=-1, device=try_gpu()):
    lines = read_sequence(directory, file)
    tokens = tokenize(lines)
    digits = torch.tensor([float(d) for d in np.array(tokens)[:, 0]], dtype=torch.float32).unsqueeze(1).to(device)
    tokens = [list(_) for _ in np.array(tokens)[:, 1]]
    vocab = Vocab(tokens)
    corpus = torch.tensor([[vocab[token] for token in line]for line in tokens], dtype=torch.long).to(device)
    
    if max_tokens > 0:
        corpus = corpus[:max_tokens]

    num_batches = len(digits) // batch_size

    # Construct dataset
    # if you want: y_scale = log(y + 1e-9)
    # dataset = (corpus[:num_batches * batch_size], torch.log(digits[:num_batches * batch_size] + 1e-9))
    dataset = (corpus[:num_batches * batch_size], digits[:num_batches * batch_size] + 1e-9)
    dataset = data.TensorDataset(*dataset)
    data_iter = data.DataLoader(dataset, batch_size, shuffle=is_train)

    return data_iter, vocab

class PositionWiseFFN(nn.Module):
    def __init__(self, ffn_num_input, ffn_num_hiddens, ffn_num_outputs):
        super(PositionWiseFFN, self).__init__()
        self.dense1 = nn.Linear(ffn_num_input, ffn_num_hiddens)
        self.relu = nn.ReLU()
        self.dense2 = nn.Linear(ffn_num_hiddens, ffn_num_outputs)
        
    def forward(self, X):
        return self.dense2(self.relu(self.dense1(X)))

class RNNModel(nn.Module):
    def __init__(self, rnn_layer, vocab_size, ffn_num_hiddens, 
                 linear1_num_hiddens, linear2_num_hiddens, linear3_num_hiddens,
                 dropout, **kwargs):
        
        super(RNNModel, self).__init__()
        self.vocab_size = vocab_size
        self.rnn = rnn_layer
        self.num_hiddens = self.rnn.hidden_size
        self.relu = nn.ReLU()
        
        self.ffn = PositionWiseFFN(self.num_hiddens, ffn_num_hiddens, self.num_hiddens)
        self.dropout = nn.Dropout(dropout)
        
        self.linear1 = nn.Linear(self.num_hiddens, linear1_num_hiddens)
        
        self.bn = nn.BatchNorm2d(linear1_num_hiddens)
        self.linear2 = nn.Linear(linear1_num_hiddens,linear2_num_hiddens)
        self.linear3 = nn.Linear(linear2_num_hiddens,linear3_num_hiddens)
        self.output = nn.Linear(linear3_num_hiddens,1)
    
    def forward(self, inputs, state):
        # inputs --> (batch_size, num_steps)
        X = F.one_hot(inputs.T.long(), self.vocab_size)
        # X --> (num_steps, batch_size, features)
        X = X.to(torch.float32)
        Y, state = self.rnn(X, state)
        # Y --> (num_steps, batch_size, hidden_size)
        Y = Y.permute(1, 0, 2)
        # Y --> (batch_size, num_steps, hidden_size)
        Y = Y.sum(dim = 1) #在num_steps维度上合并
        
        output = self.linear1(self.dropout(self.ffn(Y))+Y).unsqueeze(2).unsqueeze(3)
        output = self.linear2(self.relu(self.bn(output).squeeze(3).squeeze(2)))
        output = self.output(self.relu(self.linear3(output)))
        
        return output, state
    
    def begin_state(self, device, batch_size=1):
        if not isinstance(self.rnn, nn.LSTM):
            # GRU & RNN
            return torch.zeros((self.rnn.num_layers, 
                                batch_size, self.num_hiddens), device = device)
        
        else:
            # LSTM
            return (torch.zeros((
                self.rnn.num_layers,
                batch_size, self.num_hiddens), device=device),
                    torch.zeros((
                        self.rnn.num_layers,
                        batch_size, self.num_hiddens), device=device))

num_hiddens, num_layers = 64, 1

directory = './datanet/'
file = 'student_train.log' 
lines = read_sequence(directory, file)
tokens = tokenize(lines)
tokens = [list(_) for _ in np.array(tokens)[:, 1]]
vocab = Vocab(tokens)

#len(vocab) = 21
rnn_layer = nn.LSTM(len(vocab), num_hiddens, num_layers)
ffn_num_hiddens, linear1_num_hiddens, linear2_num_hiddens, linear3_num_hiddens = 32, 32, 32, 8

clone = RNNModel(rnn_layer, len(vocab), ffn_num_hiddens, 
               linear1_num_hiddens, linear2_num_hiddens, linear3_num_hiddens, dropout=0.5)
clone.load_state_dict(torch.load('./RNN_model.params'))
def predict_single(seq, net, vocab, device):
    net.eval()
    state = net.begin_state(batch_size=1, device = device)
    seq = torch.tensor(vocab[[i for i in seq]], device = device).unsqueeze(0)
    # print(seq)
    with torch.no_grad():
        y, _ = net(seq,state)
    return y.squeeze(0,1).numpy()

# Example
# a = predict_single('SLASYLRSM', clone, vocab, try_gpu())

################################################################


# know the pep-CDR residue-residue correlation & peptide-sequence & CDR residues
# use peptide-sequence get the CDR-residue-probability to be the loss
# compare the loss map with the CDR residues


def explicit_pep_TCR_loss(pep_seq, cdr_seq):
    pass


def IEDB_loss(iedb_freqmat, pep_seq):
    # iedb_freqmat line --> 20aa, column --> 9 position
    dictionary = dict(zip(AA_freq.keys(), np.linspace(0, 19, 20, dtype=int)))

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


def pep_TCR_loss(prob_ref_freqmat, pep_seq, cdr_seq, TCRpos):
    # prob_ref_freqmat line --> 20aa, column --> 9 position
    if cdr_seq == None:
        pep_cdr_discrete = [
            0.4000000000000001,
            3.3333333333333335,
            1.6333333333333335,
            2.0,
            0.14285714285714285,
            1.2666666666666664,
            0.8437500000000001,
            0.4129032258064516,
            1.4,
            1.5318181818181817,
            0.4273584905660377,
            3.36,
            1.5727272727272728,
            1.3707317073170735,
            1.4407407407407407,
            0.15384615384615385,
            1.3511627906976744,
            2.5333333333333337,
            2.032258064516129,
            0.5537313432835821,
        ]
        tcr_freqmat = np.zeros((args.l, 20))
        for i in range(args.l):
            tcr_freqmat[i] = prob_ref_freqmat.T[i] * pep_cdr_discrete
        tcr_freqmat = normalize(tcr_freqmat.T, args.l)
        dictionary = dict(zip(AA_freq.keys(), np.linspace(0, 19, 20, dtype=int)))
        TCRpos = list(map(int, TCRpos.split("_")))

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
        return explicit_pep_TCR_loss(pep_seq, cdr_seq)


def compute_loss(
    iedb_freqmat,
    prob_ref_freqmat,
    pep_seq,
    logfile,
    w1=10,
    w2=-0.1,
    w3=None,
    w4=args.weight_iedb,
):
    if args.cdr_sequence == None:
        w3 = args.weight_cdr_dis
    else:
        w3 = args.weight_cdr

    cdr_loss = pep_TCR_loss(
        prob_ref_freqmat,
        pep_seq=pep_seq,
        cdr_seq=args.cdr_sequence,
        TCRpos=args.TCR_loss_pos,
    )
    iedb_loss = IEDB_loss(iedb_freqmat, pep_seq)

    print("", file=logfile)
    print("LOSS RESULTS of the {}:".format(pep_seq), file=logfile)

    if args.cdr_sequence == None:
        print("****cdr loss --> {} given no cdr".format(w3 * cdr_loss), file=logfile)
    else:
        print(
            "****cdr loss --> {} with cdr sequence as {}".format(
                w3 * cdr_loss, args.cdr_sequence
            ),
            file=logfile,
        )
    print("****iedb loss --> {}".format(w4 * iedb_loss), file=logfile)

    # return (
    #     w3 * cdr_loss + w4 * iedb_loss
    # )  # DOPE is negative, thus w2 should be negative
    return predict_single(pep_seq, clone, vocab, try_gpu())