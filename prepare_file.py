#!usr/bin/env

"""
    This file takes as input the N0.tsv file of Phylogenetical_Hierarchical_Orthogroups
    folder generated from OrthoFinder2 pipeline and creates the input in the script hogsToFasta.py
    Run this script from terminal and give the path of initial file as external variable in terminal
"""

import sys
import pandas as pd

tsv_file = sys.argv[1]
tsv_file2df = pd.read_csv(tsv_file, sep='\t', header=0, index_col=False) # Initial input file to dataframe
tsv_file2df = tsv_file2df.set_index('HOG') # Set the first column as index (HOG name)
tsv_file2df = tsv_file2df.drop(['OG', 'Gene Tree Parent Clade'], axis = 1) # Drop the columns that are not needed


species_to_delete = ["Sphaeramia_orbicularis","Myripristis_murdjan", \
                    "Gadus_morhua","Danio_rerio", "Triplophysa_tibetana", \
                    "Astyanax_mexicanus", "Electrophorus_electricus",      \
                    "Pangasianodon_hypophthalmus", "Lepisosteus_oculatus"] # The species to delete from the initial analysis.
                                                                           # We kept the 24 closest to Sparus aurata

species_of_interest = tsv_file2df.drop(species_to_delete, axis = 1) # Final dataframe from initial file

saurata_no_null = species_of_interest.dropna(axis=0, subset=['Sparus_aurata_KLR_BA']) # Delete the rows that Sparus aurata has null values
                                                                                      # Change accordingly to your species of interest

perc = 80.0 # Percentage of species present in each gene family (HOG)
min_count =  int(((100-perc)/100)*saurata_no_null.shape[1] + 1)
hogs_80 = saurata_no_null.dropna( axis=0,thresh=saurata_no_null.shape[1]-min_count) # HOGs that contain values in more than 80% of species (19 or more of 24 species in total)

# Count how many species present per gene family (HOG)

# hogs_count = hogs_80.apply(lambda x: x.count(), axis=1)
# hogs_count.to_excel('/path/to/folder/species_per_hog.xlsx')

# Write df to csv (csv file will be the input file to hogsToFasta.py)
#clean_file = "/Users/klara_el/bioinfo/Thesis/Duplications/test/N0_clean.tsv"
clean_file = "/path/to/folder/N0_clean.tsv"
new_tsv = hogs_80.to_csv(clean_file, index=True, sep="\t", na_rep='', float_format=None,header=True, index_label=None, \
                        mode='w',encoding=None, date_format=None, decimal='.')
