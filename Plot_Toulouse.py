#!/usr/bin/env python
import pickle

import pandas as pd
import numpy as np
import math as m
import networkx as nx
import matplotlib.pyplot as plt
from tqdm import tqdm

## Excel, découpage,...
data = pd.read_excel(r"C:\Users\steun\PycharmProjects\UCIS-Project\filaire-de-voirie.xlsx")
T = pd.DataFrame(data, columns=['Geo Shape', 'longueur'])
T = np.array(T)
A = ['a', 'z', 'e', 'r', 't', 'y', 'u', 'i', 'o', 'p', 'q', 's', 'd', 'f', 'g', 'h', 'j', 'k', 'l', 'm', 'w', 'x', 'c',
     'v', 'b', 'n', 'M', 'L', 'S', '{', '}', '"', ':']
N = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0']

T = list(T)
print("Formatting data...")

B = [list(k) for k in tqdm(T)]

L = []

print("Cleaning bad cells...")
for i in tqdm(B):
    if type(i[0]) == float or i[0] == 'nan':
        B.remove(i)

print("Removing unwanted characters...")
for i in tqdm(B):
    s = []
    for j in list(i[0]):
        if j not in A:
            s.append(j)
    s = "".join(s)
    L.append([s, i[1]])

F = []


def fonctiontropfortetropopti(L):
    i = 0
    M = ''
    C = [[], []]
    i = 3
    while L[0][i] != ',':
        M += L[0][i]
        i += 1
    M = list(M)
    M = ''.join(M)
    C[0].append(float(M))
    M = ''
    i += 1
    while L[0][i] != ']':
        M += L[0][i]
        i += 1
    M = list(M)
    M = ''.join(M)
    C[0].append(float(M))
    M = ''
    i = len(L[0]) - 5
    while L[0][i] != ',':
        M += L[0][i]
        i -= 1
    M = list(M)
    M.reverse()
    M = ''.join(M)
    C[1].append(float(M))
    i -= 1
    M = ''
    while L[0][i] != '[':
        M += L[0][i]
        i -= 1
    M = list(M)
    M.reverse()
    M = ''.join(M)
    C[1].append(float(M))
    C.append(float(L[1]))
    C[1][0], C[1][1] = C[1][1], C[1][0]
    return C


print("Removing intermediate coordinates...")
for j in tqdm(L):
    F.append(fonctiontropfortetropopti(j))
print("Done!")
print("-" * 25 + "\n")

##Mat Array
C = []  # Juste les coordonnées des sommets
print("Getting coordinates...")
for i in tqdm(F):
    if i[0] not in C:
        C.append(i[0])
    if i[1] not in C:
        C.append(i[1])
n = len(C)  # nbre de sommets distincts

print("Computing matrix cells...")
M = np.zeros([n, n])
print("Cleaning matrix...")

for i in tqdm(F):
    a = C.index(i[0])
    b = C.index(i[1])
    if a != b:  # Routes qui font des boucles??
        M[a][b] = i[2]
        M[b][a] = i[2]
# matrice d'adjacence
print("Done!")
print("-" * 25 + "\n")

## Affichage
P = [k for k in range(n)]
G = nx.Graph()
G.add_nodes_from(P)
D = []
print("Creating visualization...")
M = list(M)
for i in tqdm(range(n)):
    for j in range(i, n):
        if M[i][j] != 0:
            D.append((i, j, M[i][j]))

G.add_weighted_edges_from(D)

labels_edges = {}
print("Matching labels and edges...")
labels_edges = {edge: G.edges[edge]['weight'] for edge in tqdm(G.edges)}
posi = {}
for k in P:
    posi[k] = C[k]
##
# nx.draw_networkx_edges(G,posi)
print("Done!")
# plt.savefig('Graphe_Toulouse')
# plt.show()

pickle.dump(G, open('C:\\Users\\steun\\PycharmProjects\\UCIS-Project\\Graphe', 'wb'))
pickle.dump(posi, open('/Visualisation/Posi', 'wb'))
