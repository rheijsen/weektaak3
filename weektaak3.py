import matplotlib.pyplot as plt
import pandas as pd
from textwrap import wrap


def lees_inhoud(bestand_naam):
    with open(bestand_naam) as bestand:
        header = bestand.readline()
        seq = ''.join(line.rstrip() for line in bestand)
    return seq


def dna_to_mrna(seq):
    seq = seq.upper()
    mrna = ""
    for i in seq:
        if i == 'T':
            mrna += 'U'
        else:
            mrna += i
    return mrna


def split_in_codons(mrna):
    codon_list = wrap(mrna, 3)
    return codon_list


def count_codons(codon_list):
    count_set = [[i, codon_list.count(i)] for i in set(codon_list)]
    return count_set


def add_count_to_dictonary(count_set):
    aa3 = {"Ala": {"GCU": 0, "GCC": 0, "GCA": 0, "GCG": 0},
           "Arg": {"CGU": 0, "CGC": 0, "CGA": 0, "CGG": 0, "AGA": 0, "AGG": 0},
           "Asn": {"AAU": 0, "AAC": 0},
           "Asp": {"GAU": 0, "GAC": 0},
           "Cys": {"UGU": 0, "UGC": 0},
           "Gln": {"CAA": 0, "CAG": 0},
           "Glu": {"GAA": 0, "GAG": 0},
           "Gly": {"GGU": 0, "GGC": 0, "GGA": 0, "GGG": 0},
           "His": {"CAU": 0, "CAC": 0},
           "Ile": {"AUU": 0, "AUC": 0, "AUA": 0},
           "Leu": {"UUA": 0, "UUG": 0, "CUU": 0, "CUC": 0, "CUA": 0, "CUG": 0},
           "Lys": {"AAA": 0, "AAG": 0},
           "Met": {"AUG": 0},
           "Phe": {"UUU": 0, "UUC": 0},
           "Pro": {"CCU": 0, "CCC": 0, "CCA": 0, "CCG": 0},
           "Ser": {"UCU": 0, "UCC": 0, "UCA": 0, "UCG": 0, "AGU": 0, "AGC": 0},
           "Thr": {"ACU": 0, "ACC": 0, "ACA": 0, "ACG": 0},
           "Trp": {"UGG": 0},
           "Tyr": {"UAU": 0, "UAC": 0},
           "Val": {"GUU": 0, "GUC": 0, "GUA": 0, "GUG": 0},
           "Stop": {"UAG": 0, "UGA": 0, "UAA": 0}
           }
    for value in count_set:
        for i in aa3:
            for j, k in aa3[i].items():
                if value[0] == j:
                    aa3[i][j] = value[1]
    return aa3


def bar_graph(aa3, naam):
    bar = pd.DataFrame(aa3)
    bar.plot(kind="bar", stacked=True, title="Codon usage\n" + naam, figsize=(6, 6))
    plt.xlabel("Codon")
    plt.ylabel("Frequentie")
    plt.legend(loc=(0.97 , 0.08))
    plt.show()


if __name__ == '__main__':
    bestand_naam = 'SIVmnd-2'
    seq = lees_inhoud(bestand_naam)
    mrna = dna_to_mrna(seq)
    codon_list = split_in_codons(mrna)
    count_set = count_codons(codon_list)
    aa3 = add_count_to_dictonary(count_set)
    bar_graph(aa3, bestand_naam)
