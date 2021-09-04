import glob
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import math


def violin_generax(path_to_speciesEventsCounts, path_to_generax_figure):
    # Place in a folder all the *_speciesEventCounts.txt files from generax output
    path = path_to_speciesEventsCounts

    txt_paths = glob.glob(os.path.join(path, '*.txt'))

    species = ['Acanthopagruslatus','ArgyrosomusregiusBA','Bettasplendens','Channaargus',
    'Collichthyslucidus','Dicentrarchuslabrax','Gasterosteusaculeatus','Gymnodracoacuticeps',
    'Kryptolebiasmarmoratus','Larimichthyscrocea','Latescalcarifer','Molamola',
    'Monopterusalbus','Oreochromisniloticus','Oryziaslatipes','Percafluviatilis',
    'Poeciliaformosa','Scophthalmusmaximus','Serioladumerili','SparusaurataKLRBA',
    'Takifugurubripes','Tetraodonnigroviridis','Xiphophorusmaculatus','Zebrambuna']



    dup_perSpecies_perFamily = dict()

    for txt in txt_paths:
        txt_name = os.path.basename(txt)
        split_txtName = txt_name.split('_')
        hog_name = split_txtName[0] #keep the name of the gene family

        with open(txt,"r") as f:

            lines = f.readlines()
            no_ofDups_list = []
            for line in lines:
                if line.startswith("node"): # Keep only the leaves, not the inner nodes of gene trees
                    continue

                species_name = line.split()[0].rstrip() #Keep species name for headers
                no_ofDups = line.split()[2].rstrip() #Keep number of duplications per species

                #if species_name in species:
                no_ofDups_list.append(int(no_ofDups)) #Keep the duplications for each HOG

        dup_perSpecies_perFamily[hog_name]= no_ofDups_list

    df = pd.DataFrame(data=dup_perSpecies_perFamily, index = species)
    df = df.T
    df.rename(columns = {'ArgyrosomusregiusBA': 'Argyrosomusregius','SparusaurataKLRBA':'Sparusaurata'}, inplace = True)

    species_rename = ['Acanthopagrus latus','Argyrosomus regius','Betta splendens','Channa argus',
        'Collichthys lucidus','Dicentrarchus labrax','Gasterosteus aculeatus','Gymnodraco acuticeps',
        'Kryptolebias marmoratus','Larimichthys crocea','Lates calcarifer','Mola mola',
        'Monopterus albus','Oreochromis niloticus','Oryzias latipes','Perca fluviatilis',
        'Poecilia formosa','Scophthalmus maximus','Seriola dumerili','Sparus aurata',
        'Takifugu rubripes','Tetraodon nigroviridis','Xiphophorus maculatus','Zebra mbuna']

    # Create the violin plot
    fig, axs = plt.subplots(figsize=(20, 30))
    count=1

    colors = sns.color_palette("pastel", 25)

    for species_name in df.columns:

        new_df = df[species_name]
        new_df = new_df.replace(0,np.nan)
        new_df = new_df.replace(1,np.nan)
        new_df = new_df.dropna() # Ignore the species with 0 or 1 duplication per HOG
        w = math.log(max(new_df.value_counts()),2)

        plt.subplot(12,2,count)
        sns.violinplot(data=new_df,inner="points",saturation=.25,color=colors[count],width=w,orient='h')


        for name in species_rename: # Convert the species name for y axis
            if name.split(" ")[0] + name.split(" ")[1] == species_name:
                plt.ylabel(name)

        plt.xlabel('# Duplications per HOG')
        plt.xlim(-40,100)
        plt.ylim(-5,5)
        count+=1

    plt.subplots_adjust(hspace=0.7, wspace=0.2)
    fig.suptitle('Gene duplication events per species',y=0.91, fontsize=22)
    plt.savefig('path_to_generax_figure')
    plt.show()
    return





def violin_cafe(path_to_fams_txt, path_to_cafe_figure):

    skip_row = []
    new_line = []


    # The summary_cafe_fams.txt output file from cafe reports
    cafe_fams = path_to_fams_txt
    with open(cafe_fams, "r") as cafe:
        lines = cafe.readlines()[2:]
        for line in lines:
            if line.startswith("<"):
                skip_row.append(line.rstrip().split(":\t"))
                continue
            new_line.append(line.rstrip().split(","))

    df = pd.read_csv(cafe_fams, sep="\t",index_col=0, names=["Families"]).iloc[2:]
    df_split = df.Families.str.split(",",expand=True)
    df_new = df_split.loc[(~df.index.str.startswith('<')),:].T

    replacements={'N0.':"", '\[': ':', '\]': ''}
    df_clean = df_new.replace(replacements,regex=True)

    #Sparus col in descending order
    df_sth = df_clean.set_index('Sparus<28>:')
    df_sth = df_sth.reindex(sorted(df_sth.index, key=lambda x: int(x.split(':')[-1].split("*")[0]),reverse=True)).reset_index()

    pd.set_option('display.max_columns', 500)
    pd.set_option('display.max_rows', 500)

    cafe_dict = df_sth.to_dict('list')

    species_dict={}
    for name,expansions in cafe_dict.items():
        expansions_dict = {}
        for expansion in expansions:
            if expansion==None:
                continue
            hog = expansion.split(":")[0]
            exp = int(expansion.split(":")[1].split("*")[0])
            if exp>0:
                expansions_dict[hog] = exp

        species_dict[name] = expansions_dict

    cafe_dataframe = pd.DataFrame.from_dict(species_dict)



    species_rename = ['Acanthopagrus latus','Argyrosomus regius','Betta splendens','Channa argus',
        'Collichthys lucidus','Dicentrarchus labrax','Gasterosteus aculeatus','Gymnodraco acuticeps',
        'Kryptolebias marmoratus','Larimichthys crocea','Lates calcarifer','Mola mola',
        'Monopterus albus','Oreochromis niloticus','Oryzias latipes','Perca fluviatilis',
        'Poecilia formosa','Scophthalmus maximus','Seriola dumerili','Sparus aurata',
        'Takifugu rubripes','Tetraodon nigroviridis','Xiphophorus maculatus','Zebra mbuna']

    fig, axs = plt.subplots(figsize=(20, 30))
    count=1

    colors = sns.color_palette("pastel", 25)

    for species_name in cafe_dataframe.columns:

        new_df = cafe_dataframe[species_name]
        new_df = new_df.replace(0,np.nan)
        new_df = new_df.replace(1,np.nan)
        new_df = new_df.dropna()
        w = math.log(max(new_df.value_counts()),2)

        plt.subplot(12,2,count)
        sns.violinplot(data=new_df,inner="points",saturation=.25,color=colors[count],width=w,orient='h')

        for name in species_rename:
            if species_name.split("<")[0] == name.split(" ")[0]:
                plt.ylabel(name)
        plt.xlabel('# Expansions per HOG')
        plt.xlim(-40,100)
        plt.ylim(-5,5)
        count+=1

    plt.subplots_adjust(hspace=0.7, wspace=0.2)
    fig.suptitle('Gene family expansions per species',y=0.91, fontsize=22)
    plt.savefig(path_to_cafe_figure)
    plt.show()

    return


def main():

    # Violin plot for GeneRax
    path_to_speciesEventsCounts='/path/to/speciesEventCounts/folder/'
    path_to_generax_figure = '/path/to/save/violin/plot/generax_violins.png'
    violin_generax(path_to_speciesEventsCounts, path_to_generax_figure)

    # Violin plot for CAFE
    path_to_fams_txt = "/path/to/CAFE_report/summary_cafe_fams.txt"
    path_to_cafe_figure = '/path/to/save/violin/plot/cafe_violins.png'
    violin_cafe(path_to_fams_txt, path_to_cafe_figure)

    return

main()
