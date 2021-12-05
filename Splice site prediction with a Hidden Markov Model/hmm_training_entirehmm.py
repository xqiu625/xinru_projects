import csv
import numpy as np

DNA_file = open('DNA_sequence.txt','r')
all_sequence = {}
DNA_key = []

for i in DNA_file:
    if i[0] == '>':
        all_sequence[i[1:].replace('\n', '')] = ''
        key = i[1:].replace('\n', '')
        DNA_key.append(i[1:].replace('\n', ''))
    else:
        all_sequence[key] += i.replace('\n', '')

DNA_file.close()

exon_info = {}
EXON_file = open('exon_file.txt','r')
for i in EXON_file:
    line = i.replace('\r\n', '').split(' ')
    exon_info[line[0]] = {'exon':[],'intron':[]}
    for i in range(len(line[1:]) -1):
        if i%2 == 0:
            exon_info[line[0]]['exon'].append((int(line[1:][i]) -1,int(line[1:][i+1])))
        else:
            exon_info[line[0]]['intron'].append((int(line[1:][i]),int(line[1:][i+1]) -1))
    
    

#print all_sequence['ACU08131']
#print len(all_sequence.keys())
#print exon_info['ACU08131']
#print all_sequence['ACU08131'][520:641]
#print all_sequence['ACU08131'][641:1065]
#print all_sequence['ACU08131'][1362:1859]

#### check if any sequence need to be reversed:  None
for key in exon_info:
    for (x,y) in exon_info[key]['exon']:
        if x > y: print key


#print DNA_key
states = ('e', 'd1', 'd2', 'i', 'a1', 'a2')
states_freq = {'e':0, 'd1':0, 'd2':0, 'i':0, 'a1':0, 'a2':0}
transition_freq = {}
emission_freq = {}
transition_prob = {}
emission_prob = {}
for i in states:
    transition_freq[i] = {'e':0, 'd1':0, 'd2':0, 'i':0, 'a1':0, 'a2':0}
    emission_freq[i] = {'A':0, 'C':0, 'G':0, 'T':0}
    transition_prob[i] = {'e':0, 'd1':0, 'd2':0, 'i':0, 'a1':0, 'a2':0}
    emission_prob[i] = {'A':0, 'C':0, 'G':0, 'T':0}


DNA_train = DNA_key[100:]
for gene in DNA_train:
    state_sequence = ['e' for _ in range(len(all_sequence[gene]))]
    for (x,y) in exon_info[gene]['intron']:
        for i in range(x,y):
            state_sequence[i] = 'i'
        if all_sequence[gene][x: x+2] in ('GT','AC'):
            state_sequence[x] = 'd1'
            state_sequence[x+1] = 'd2'
            if all_sequence[gene][x: x+2] in ('AC'): print gene
        if all_sequence[gene][y-2: y] in ('AG','CT'):
            state_sequence[y-2] = 'a1'
            state_sequence[y-1] = 'a2'
    for k in range(min(min(exon_info[gene]['exon'])), max(max(exon_info[gene]['exon']))-1):
        #print gene
        states_freq[state_sequence[k]] += 1
        transition_freq[state_sequence[k]][state_sequence[k+1]] += 1
        if all_sequence[gene][k] not in ('A', 'C', 'G', 'T'):
            for s in ['A', 'C', 'G', 'T']:
                emission_freq[state_sequence[k]][s] += 0.25
        else:
            emission_freq[state_sequence[k]][all_sequence[gene][k]] += 1


for i in states:
    for j in states:
        transition_prob[i][j] = float(transition_freq[i][j])/states_freq[i]
    for k in ['A', 'C', 'G', 'T']:
        emission_prob[i][k] = float(emission_freq[i][k])/states_freq[i]



transition_prob = {
   'e'      : {'e': 0.90, 'd1': 0.10, 'd2': 0.00, 'i': 0.00, 'a1': 0.00, 'a2': 0.00},
   'd1'    : {'e': 0.00, 'd1': 0.00, 'd2': 1.00, 'i': 0.00, 'a1': 0.00, 'a2': 0.00},
   'd2'    : {'e': 0.00, 'd1': 0.00, 'd2': 0.00, 'i': 1.00, 'a1': 0.00, 'a2': 0.00},
   'i'    : {'e': 0.00, 'd1': 0.00, 'd2': 0.00, 'i': 0.90, 'a1': 0.10, 'a2': 0.00},
   'a1' : {'e': 0.00, 'd1': 0.00, 'd2': 0.00, 'i': 0.00, 'a1': 0.00, 'a2': 1.00},
   'a2' : {'e': 1.00, 'd1': 0.00, 'd2': 0.00, 'i': 0.00, 'a1': 0.00, 'a2': 0.00}
   }

emission_prob = {
   'e'      : {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25},
   'd1'    : {'A': 0.05, 'C': 0.00, 'G': 0.95, 'T': 0.00},
   'd2'    : {'A': 0.00, 'C': 0.05, 'G': 0.00, 'T': 0.95},
   'i'    : {'A': 0.40, 'C': 0.10, 'G': 0.10, 'T': 0.40},
   'a1' : {'A': 0.95, 'C': 0.00, 'G': 0.05, 'T': 0.00},
   'a2' : {'A': 0.05, 'C': 0.00, 'G': 0.95, 'T': 0.00}
   }


#print states_freq
#print transition_freq
print emission_freq
print transition_prob
print emission_prob

def viterbi(obs, states, start_p, trans_p, emit_p):
    V = [{}]
    path = {}
 
    # Initialize base cases (t == 0)
    for y in states:
        V[0][y] = start_p[y] * emit_p[y][obs[0]]
        path[y] = [y]
 
    # Run Viterbi for t > 0
    for t in range(1, len(obs)):
        V.append({})
        newpath = {}
 
        for y in states:
            (prob, state) = max((V[t-1][y0] * trans_p[y0][y] * emit_p[y][obs[t]], y0) for y0 in states)
            V[t][y] = prob
            newpath[y] = path[state] + [y]
        path = newpath
    n = 0           # if only one element is observed max is sought in the initialization values
    if len(obs) != 1:
        n = t
    #print_dptable(V)
    (prob, state) = max((V[n][y], y) for y in states)
    return (prob, path[state])
    #return path[state]




start_probability = {'e': 0.50, 'd1': 0.00, 'd2': 0.00, 'i': 0.50, 'a1': 0.00, 'a2': 0.00}

test = DNA_key[0]
observations = all_sequence[test][520:]
print 'ssss', ''.join(viterbi(observations,states,start_probability,transition_prob,emission_prob)[1])
predict = viterbi(observations,states,start_probability,transition_prob,emission_prob)[1]
for i in range(len(predict)):
    if predict[i] == 'd1': print i

symbols = ['A', 'C', 'G', 'T']
states = ('e', 'd1', 'd2', 'i', 'a1', 'a2')
symbol_index = {'A':0, 'C':1, 'G':2, 'T':3}
states_index = {'e':0, 'd1':1, 'd2':2, 'i':3, 'a1':4, 'a2':5}

start_prob = np.array([0.50, 0.00, 0.00, 0.50, 0.00, 0.00])

emission_mat = np.zeros((len(states),len(symbols)))
for i in states_index:
    for j in symbol_index:
        emission_mat[states_index[i]][symbol_index[j]] = emission_prob[i][j]

transition_mat = np.zeros((len(states),len(states)))
for i in states_index:
    for j in states_index:
        transition_mat[states_index[i]][states_index[j]] = transition_prob[i][j]

#print transition_mat
#print emission_mat

#import random
#import hmm
#import time
#import sys

#HMM = hmm.HMM(symbols, states, start_prob, transition_mat, emission_mat)
#predict =  HMM.viterbi(observations)
#print predict.index('d1')

