"""
gcpy: python utilities for chromatographs
Copyright (C) 2021 David Ollodart <LICENSE>

Provided a networkx representation of the structure, determine the
chemical identity of groups.  Most groups are defined by at most two
edges of separation. For example a carboxyl group is

  O-H
  |
R-C=0

which has two edges of separation--you need to know the oxygen in the
alcohol is not part of, e.g., an ether. This classifies so-called
moeities which are at most needing structure extending to second nearest
neighbors.

"""

from pysmiles.read_smiles import read_smiles
from collections import defaultdict

ary_dct = {0: '', 1: 'primary', 2: 'secondary', 3: 'tertiary'}
alk_dct = {2: 'alkyne', 3: 'alkene', 4: 'alkane'}


def smiles2carbontypes(smiles):
    graph = read_smiles(smiles, explicit_hydrogen=True)
    dldct = graph.nodes(data='element')
    count = 0
    for node0 in graph.nodes:

        if dldct[node0] != 'C':
            count += 1
            continue
            # only C atoms needed for ECN

        nbors1 = ''
        nbors2 = ''
        edge01 = []
        edge12 = []
        for node1 in graph.neighbors(node0):
            nbors1 += dldct[node1]
            edge01.append(node1)
            tmp = []
            for node2 in graph.neighbors(node1):
                nbors2 += dldct[node2]
                if node2 != node0:
                    tmp.append(node2)
            edge12.append(tmp)

        nbors1 = ''.join(sorted(nbors1))
        nbors2 = ''.join(sorted(nbors2))

        ary = nbors1.count('C')
        hyd = nbors1.count('H')

        ln1 = len(nbors1)
        if ln1 == ary + hyd:
            yield ary_dct[ary], alk_dct[ln1]

        if nbors1 == 'CCO':
            yield '', 'ketone'
        if nbors1 == 'CHO':
            yield '', 'aldehyde'

        if nbors1 == 'CHHO' or nbors1 == 'CCHO':
            # look through (could be an ether)
            # won't detect epoxides
            is_ether = True
            for i in range(len(edge01)):
                if dldct[edge01[i]] == 'O':
                    n = edge12[i]
                    if len(n) == 1 and dldct[n[0]] == 'H':
                        is_ether = False
                        yield ary_dct[ary], 'alcohol'
            if is_ether:
                yield '', 'ether'

        if nbors1 == 'COO':
            is_anhydride = True
            for i in range(len(edge01)):
                if dldct[edge01[i]] == 'O':
                    n = edge12[i]
                    if len(n) == 1 and dldct[n[0]] == 'H':
                        is_anhydride = False
                        yield '', 'carboxylic acid'
            if is_anhydride:
                yield '', 'acid anhydride or ester'
                # requires 3rd nearest neighbors to discriminate

        if 'N' in nbors1:
            if nbors1 == 'CNO':
                yield '', 'amide'
            elif nbors1 == 'CN':
                yield '', 'nitrile'
            else:
                for i in range(len(edge01)):
                    if dldct[edge01[i]] == 'N':
                        n = edge12[i]
                        is_amine = True
                        amine_ary = 1
                        for j in n:
                            if dldct[j] == 'C':
                                is_amine = False
                                yield '', 'amine, or carbon adjacent to amide'
                            elif dldct[j] == 'C':
                                amine_ary += 1
                        if is_amine:
                            yield ary_dct[amine_ary], 'amine'

        if 'S' in nbors1:
            for i in range(len(edge01)):
                if dldct[edge01[i]] == 'S':
                    n = edge12[i]
                    if len(n) == 1 and dldct[n[0]] == 'H':
                        yield '', 'thiol'
                    else:
                        yield '', 'sulfide'


class_dct = {'alkane': 'aliphatic',
             'alkene': 'olefinic',
             'alkyne': 'acetylinic',
             'carboxylic acid': 'carboxyl',
             'aldehyde': 'carbonyl',
             'ketone': 'carbonyl',
             'ester': 'carboxyl',
             'nitrile': 'nitrile',
             'ether': 'ether',
             'amine': 'amine',
             'alcohol': 'alcohol',
             'acid anhydride': 'ester',  # may be carboxyl
             'amide': 'aliphatic',
             'acid anhydride or ester': 'ester'  # ambiguous class
            }


def factory():
    return 0


class_dct = defaultdict(factory, class_dct)
u = 0

# ecn appears to be based on functional groups defined not just by
# carbons but by heteroatoms or even terminal groups add 1 to all cases
# where a functional group is defined as something other than a carbon
# atom giving carbon atom classing (may be incorrect, articles are
# paywalled)

ecn_dct = {"aliphatic": 1.00,
           "aromatic": 1.00,
           "olefinic": 0.95,
           "acetylinic": 1.30,
           "carbonyl": u,
           "carboxyl": u,
           "nitrile": 0.30,
           "ether": -1.00 + 1.5,  # since 2 carbon atoms will be ether classed
           "primary alcohol": -0.50 + 1,
           "secondary alcohol": -0.75 + 1,
           "tertiary alcohol": -0.25 + 1,
           "amine": u + 1
          }

ecn_dct = defaultdict(factory, ecn_dct)

if __name__ == '__main__':
    principal_functional_groups = {  # from Carey, Giuliano 9th ed
        'Ethane': 'CC',
        'Ethene': 'C=C',
        'Ethyne': 'C#C',
        '1,3-Butadiene': 'C=CC=C',
        'Benzene': 'C1=CC=CC=C1',
        'Chloroethane': 'CCCl',
        'Chlorobenzene': 'C1=CC=C(Cl)C=C1',
        'Ethanol': 'CCO',
        'Phenol': 'C1=CC=C(O)C=C1',
        'Ethoxyethane': 'CCOCC',
        'Epoxyethane': 'C1OC1',
        'Ethanal': 'CC=O',
        '2-Propanone': 'CC(=O)C',
        'Ethanoic acid': 'CC(=O)O',
        'Ethanoyl chloride': 'CC(=O)Cl',
        'Ethanoic anhydride': 'CC(=O)OC(=O)C',
        'Ethyl ethanoate': 'CC(=O)OCC',
        'N-Methylethanamide': 'CC(=O)NC',
        'Ethanamine': 'CCN',
        'Ethanenitrile': 'CC#N',
        'Nitrobenzene': 'C1=CC=C(N(O)O)C=C1',
        'Ethanethiol': 'CCS',
        'Diethyl sulfide': 'CCSCC'}

    for name in principal_functional_groups:
        ecn = 0
        smiles = principal_functional_groups[name]
        for ary, typ in smiles2carbontypes(smiles):
            if typ == 'alcohol':
                typ = ary + ' ' + typ
            ecn += ecn_dct[class_dct[typ]]
        print(name, ecn)
