#!/usr/bin/env

"""
 Run this script from terminal with 3 input paths (in this specific order):
 1. Path from N0.tsv file from OrthoFinder's result without the second and third columns
        (If you used the prepare_file.py, take the output tsv)
 2. Path from the folder with the chosen proteomes files of species
 3. Path to save the fasta files of HOGs
"""

import sys
import os.path
import os, glob
import pandas as pd
from Bio import SeqIO
from collections import defaultdict


def parse_input_file(species_per_gene_family):
    # Parses the input tsv and returns a dictionary with HOGs, species and genes
    hogsTSV=open(species_per_gene_family, 'r')

    lines = hogsTSV.readlines()[0:] # All lines of the file

    species_names = [sname.rstrip() for sname in lines[0].split("\t")[1:]] # Keep species names in a list

    genes = []

    hogs_species = {}

    for idx in range(1,len(lines)):

        lines_noTab = lines[idx].split("\t") # Remove \t char between cols

        gene_species = lines_noTab[1:] # Get only columns with species

        hog = lines_noTab[0][3:] # Get the HOG name without the file's identifier

        genes = [species.strip().split(", ") for species in gene_species] # List of lists for each column


        # Create a nested dictionary
        # Keys = gene family name
        # Values: dictionary with species names as keys and values lists of genes per species per gene family
        hogs_species[hog] = {}
        for index in range(len(species_names)):
            if genes[index] == []:
                continue
            hogs_species[hog][species_names[index]] = genes[index]
    return hogs_species


def rename_genes(hogs_species):
    # Rename gene names
    # Concatenate keys with values of hogs_species' second level dictionary

    rename_genes = defaultdict(lambda: defaultdict(dict))  #dict(dict())

    for v in hogs_species.values():

        for species,genes2 in v.items():
            species_replaced = species.replace('_', '') # No _ in species names

            new_genes = []
            for gene in genes2:
                if genes2 == [""]:
                    continue
                gene_replaced = gene.replace('_', '')# No _ in gene names
                gene_new = species_replaced + "_" + gene_replaced
                new_genes.append(gene_new)
                rename_genes[species][gene] = gene_new

            v[species] = new_genes

    return rename_genes

def fasta2dict(fil):
    # By: https://gist.github.com/jacob-ogre/5318981
    """
    Convert ONE fasta file to one dictionary with headers as keys and sequences as values

    Read fasta-format file fil, return dict of form scaffold:sequence.
    Note: Uses only the unique identifier of each sequence, rather than the
    entire header, for dict keys.
    """

    dic = {}
    cur_scaf = ''
    cur_seq = []
    for line in open(fil):
        if line.startswith(">") and cur_scaf == '':
            cur_scaf = line[1:].split('\n')[0]
        elif line.startswith(">") and cur_scaf != '':
            dic[cur_scaf] = ''.join(cur_seq)
            cur_scaf = line[1:].split('\n')[0]
            cur_seq = []
        else:
            cur_seq.append(line.rstrip())
    dic[cur_scaf] = ''.join(cur_seq)
    return dic


def ret_new_gene(gene,updated_genes):
    # Returns the gene with new name from the dictionary of rename_genes function

    for species,genes in updated_genes.items():
        for old_gene,new_gene in genes.items():
            if gene==old_gene:
                return new_gene


def ret_old_gene(gene,updated_genes_df):
    # Returns the gene with old name from the dictionary of rename_genes function

    for species,genes in updated_genes_df.items():
            for old_gene,new_gene in genes.items():
                if gene==new_gene:
                    return old_gene


def proteomes2dict(proteomes_path):
    # Creates a nested 2 levels dictionary from proteomes files
    # Keys: species names
    # Values: gene names as keys and sequences as values

    #folder_path = '/Users/klara_el/bioinfo/Thesis/Duplications/teleost_proteomes/'
    folder_path = proteomes_path
    fasta_paths = glob.glob(os.path.join(folder_path, '*.fasta'))

    proteomes_fasta = dict(dict())


    for fasta in fasta_paths:
        fasta_name = os.path.basename(fasta)
        split_fastaName = fasta_name.split('.')
        species_name = split_fastaName[0]

        proteomes_fasta[species_name] = fasta2dict(fasta)

    return proteomes_fasta


def gene_families_file(hogs_species,updated_genes,proteomes_fasta,hogs_path):
    # Creates one fasta file per gene family

    for hogs,species in hogs_species.items():

        hog_fasta = open(hogs_path + hogs +".fa" , "w")
        for species,genes in hogs_species[hogs].items(): # Iteration in each gene family
            for i in genes: # Iteration in each gene list per gene family
                if ret_old_gene(i,updated_genes) != None:
                   old = ret_old_gene(i,updated_genes)
                   hog_fasta.write(">" + i + "\n")
                   hog_fasta.write(proteomes_fasta[species][old] + "\n") # Returns the sequence
                                                                         # using the initial gene name
    return

def main():
    # Call all the above functions

    species_per_gene_family = sys.argv[1]
    proteomes_path = sys.argv[2]
    hogs_path = sys.argv[3]


    hogs_species = parse_input_file(species_per_gene_family)
    updated_genes = rename_genes(hogs_species)
    proteomes_fasta = proteomes2dict(proteomes_path)
    gene_families_file(hogs_species,updated_genes,proteomes_fasta,hogs_path)


    return

main()
