import csv
import numpy as np
import sys

####  Read Data, for both gene sequence and exon position

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
    line = i.replace('\n', '').split(' ')
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


#### Viterbi algorithm, given transition probability, emission probability, observation, predict the most likely hidden states
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

        # Don't need to remember the old paths
        path = newpath
    n = 0           # if only one element is observed max is sought in the initialization values
    if len(obs) != 1:
        n = t
    #print_dptable(V)
    (prob, state) = max((V[n][y], y) for y in states)
    return (prob, path[state])
    #return path[state]

#def print_dptable(V):
#    s = "    " + " ".join(("%7d" % i) for i in range(len(V))) + "\n"
#    for y in V[0]:
#        s += "%.5s: " % y
#        s += " ".join("%.7s" % ("%f" % v[y]) for v in V)
#        s += "\n"
#    print(s)
#


#####  Donor Model

#print DNA_key
#states = ('e', 'd1', 'd2', 'i', 'a1', 'a2')
states = ('e', 'd1', 'd2', 'i')
#states_freq = {'e':0, 'd1':0, 'd2':0, 'i':0, 'a1':0, 'a2':0}
states_freq = {'e':0, 'd1':0, 'd2':0, 'i':0}
transition_freq = {}
emission_freq = {}
transition_prob = {}
emission_prob = {}
for i in states:
    #transition_freq[i] = {'e':0, 'd1':0, 'd2':0, 'i':0, 'a1':0, 'a2':0}
    transition_freq[i] = {'e':0, 'd1':0, 'd2':0, 'i':0}
    emission_freq[i] = {'A':0, 'C':0, 'G':0, 'T':0}
    #transition_prob[i] = {'e':0, 'd1':0, 'd2':0, 'i':0, 'a1':0, 'a2':0}
    transition_prob[i] = {'e':0, 'd1':0, 'd2':0, 'i':0}
    emission_prob[i] = {'A':0, 'C':0, 'G':0, 'T':0}


DNA_train = DNA_key[100:]
for gene in DNA_train:
    for (x,y) in exon_info[gene]['exon']:
        state_sequence = ['e','e','e','i','i','i','i','i','i']
        if y < max(max(exon_info[gene]['exon'])):
            if all_sequence[gene][y: y+2] in ('GT','AC'):
                state_sequence[3] = 'd1'
                state_sequence[4] = 'd2'
            for k in range(y-3, y+5):
                states_freq[state_sequence[k-y+3]] += 1
                transition_freq[state_sequence[k-y+3]][state_sequence[k-y+4]] += 1
                if all_sequence[gene][k] not in ('A', 'C', 'G', 'T'):
                    for s in ['A', 'C', 'G', 'T']:
                        emission_freq[state_sequence[k-y+3]][s] += 0.25
                else:
                    emission_freq[state_sequence[k-y+3]][all_sequence[gene][k]] += 1


for i in states:
    for j in states:
        transition_prob[i][j] = float(transition_freq[i][j])/states_freq[i]
    for k in ['A', 'C', 'G', 'T']:
        emission_prob[i][k] = float(emission_freq[i][k])/states_freq[i]


#### for testing a given probability
#transition_prob = {
#   'e'      : {'e': 0.90, 'd1': 0.10, 'd2': 0.00, 'i': 0.00, 'a1': 0.00, 'a2': 0.00},
#   'd1'    : {'e': 0.00, 'd1': 0.00, 'd2': 1.00, 'i': 0.00, 'a1': 0.00, 'a2': 0.00},
#   'd2'    : {'e': 0.00, 'd1': 0.00, 'd2': 0.00, 'i': 1.00, 'a1': 0.00, 'a2': 0.00},
#   'i'    : {'e': 0.00, 'd1': 0.00, 'd2': 0.00, 'i': 0.90, 'a1': 0.10, 'a2': 0.00},
#   'a1' : {'e': 0.00, 'd1': 0.00, 'd2': 0.00, 'i': 0.00, 'a1': 0.00, 'a2': 1.00},
#   'a2' : {'e': 1.00, 'd1': 0.00, 'd2': 0.00, 'i': 0.00, 'a1': 0.00, 'a2': 0.00}
#   }

#emission_prob = {
#   'e'      : {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25},
#   'd1'    : {'A': 0.05, 'C': 0.00, 'G': 0.95, 'T': 0.00},
#   'd2'    : {'A': 0.00, 'C': 0.05, 'G': 0.00, 'T': 0.95},
#   'i'    : {'A': 0.40, 'C': 0.10, 'G': 0.10, 'T': 0.40},
#   'a1' : {'A': 0.95, 'C': 0.00, 'G': 0.05, 'T': 0.00},
#   'a2' : {'A': 0.05, 'C': 0.00, 'G': 0.95, 'T': 0.00}
#   }


#print states_freq
#print transition_freq
#print emission_freq
print 'Donor model prob'
print transition_prob
print emission_prob


#Below is testing
start_probability = {'e': 0.50, 'd1': 0.00, 'd2': 0.00, 'i': 0.50, 'a1': 0.00, 'a2': 0.00}


real_t = 0
real_f = 0
tp = 0
fp = 0

for test in DNA_key[0:100]:
    donor_pos = []
    for x in range(len(all_sequence[test])-9):
        observations = all_sequence[test][x:x+9]
        try:
            predict = viterbi(observations,states,start_probability,transition_prob,emission_prob)[1]
            if predict[3] == 'd1' and predict[4] == 'd2': donor_pos.append(x+3)
        except KeyError:
            pass

    #print 'predicted donor position', donor_pos

    actual_donor_pos = [y for (x,y) in exon_info[test]['exon'] if y < max(max(exon_info[test]['exon'])) ]
    #print 'actual donor position', actual_donor_pos

    for i in donor_pos:
        if i >  min(min(exon_info[test]['exon'])) and i <  max(max(exon_info[test]['exon'])):
            if i in actual_donor_pos:
                tp += 1
            else:
                fp += 1

    real_t += len(actual_donor_pos)
    real_f += max(max(exon_info[test]['exon'])) - min(min(exon_info[test]['exon'])) - len(actual_donor_pos)
fn = real_t - tp
tn = real_f - fp

print 'Donor Model Confusion Matrix'
print tp, fn
print fp, tn




####### Acceptor Model

#print DNA_key
#states = ('e', 'd1', 'd2', 'i', 'a1', 'a2')
states = ('e', 'a1', 'a2', 'i')
#states_freq = {'e':0, 'd1':0, 'd2':0, 'i':0, 'a1':0, 'a2':0}
states_freq = {'e':0, 'a1':0, 'a2':0, 'i':0}
transition_freq = {}
emission_freq = {}
transition_prob = {}
emission_prob = {}
for i in states:
    #transition_freq[i] = {'e':0, 'd1':0, 'd2':0, 'i':0, 'a1':0, 'a2':0}
    transition_freq[i] = {'e':0, 'a1':0, 'a2':0, 'i':0}
    emission_freq[i] = {'A':0, 'C':0, 'G':0, 'T':0}
    #transition_prob[i] = {'e':0, 'd1':0, 'd2':0, 'i':0, 'a1':0, 'a2':0}
    transition_prob[i] = {'e':0, 'a1':0, 'a2':0, 'i':0}
    emission_prob[i] = {'A':0, 'C':0, 'G':0, 'T':0}


DNA_train = DNA_key[100:]
for gene in DNA_train:
    for (x,y) in exon_info[gene]['intron']:
        state_sequence = ['i','i','i','i','i','i','i','i','i','i','i','i','i','i','i','e','e','e']
        if all_sequence[gene][y-2: y] in ('AG','CT'):
            state_sequence[13] = 'a1'
            state_sequence[14] = 'a2'
        #print gene, all_sequence[gene][y-15:y+3]
        #print gene,state_sequence

        for k in range(y-15, y+2):
        #print gene
            states_freq[state_sequence[k-y+15]] += 1
            transition_freq[state_sequence[k-y+15]][state_sequence[k-y+16]] += 1
            if all_sequence[gene][k] not in ('A', 'C', 'G', 'T'):
                for s in ['A', 'C', 'G', 'T']:
                    emission_freq[state_sequence[k-y+15]][s] += 0.25
            else:
                emission_freq[state_sequence[k-y+15]][all_sequence[gene][k]] += 1



for i in states:
    for j in states:
        transition_prob[i][j] = float(transition_freq[i][j])/states_freq[i]
    for k in ['A', 'C', 'G', 'T']:
        emission_prob[i][k] = float(emission_freq[i][k])/states_freq[i]

#print states_freq
#print transition_freq
#print emission_freq
print 'Acceptor model prob'
print transition_prob
print emission_prob

start_probability = {'e': 0.50, 'd1': 0.00, 'd2': 0.00, 'i': 0.50, 'a1': 0.00, 'a2': 0.00}


real_t = 0
real_f = 0
tp = 0
fp = 0

for test in DNA_key[0:100]:
    acceptor_pos = []
    for x in range(len(all_sequence[test])-18):
        observations = all_sequence[test][x:x+18]
        try:
            predict = viterbi(observations,states,start_probability,transition_prob,emission_prob)[1]
            if predict[13] == 'a1' and predict[14] == 'a2': acceptor_pos.append(x+14)
        except KeyError:
            pass

#    print 'predicted acceptor position', acceptor_pos

    actual_acceptor_pos = [y-1 for (x,y) in exon_info[test]['intron'] if y < max(max(exon_info[test]['exon'])) ]
#    print 'actual acceptor position', actual_acceptor_pos

    for i in acceptor_pos:
        if i >  min(min(exon_info[test]['exon'])) and i <  max(max(exon_info[test]['exon'])):
            if i in actual_acceptor_pos:
                tp += 1
            else:
                fp += 1

    real_t += len(actual_acceptor_pos)
    real_f += max(max(exon_info[test]['exon'])) - min(min(exon_info[test]['exon'])) - len(actual_acceptor_pos)
fn = real_t - tp
tn = real_f - fp

print 'Acceptor Model Confusion Matrix'
print tp, fn
print fp, tn

        
