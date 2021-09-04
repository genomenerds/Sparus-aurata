import pandas as pd
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2, venn2_circles


def parse_generax(path_to_speciesEventsCounts):
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
                if line.startswith("node"):
                    continue

                species_name = line.split()[0].rstrip() #Keep species name for headers
                no_ofDups = line.split()[2].rstrip() #Keep number of duplications per species


                no_ofDups_list.append(int(no_ofDups))

        dup_perSpecies_perFamily[hog_name]= no_ofDups_list

    df = pd.DataFrame(data=dup_perSpecies_perFamily, index = species)
    df = df.T

    first_col_saur = df.pop('SparusaurataKLRBA')
    df.insert(0,'SparusaurataKLRBA', first_col_saur)

    ascending_df = df.sort_values(by=['SparusaurataKLRBA'],ascending=False)
    pd.set_option('display.max_columns', 500)
    pd.set_option('display.max_rows', 500)

    df_not_null_sparus = ascending_df[:489]
    sparus_hogs_dups = list(df_not_null_sparus.index.values)


    return df_not_null_sparus,sparus_hogs_dups

def parse_cafe(path_to_fams_txt):

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
            #line = line.rstrip().split(",")

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


    sparus_list = cafe_dict['Sparus<28>:']
    hogs_expanded = []
    for orthogroup in sparus_list:
        hog = orthogroup.split(":")[0]
        number = int(orthogroup.split(":")[1].split("*")[0])
        if number>0:
            hogs_expanded.append(hog)

    expansions_dict = {}
    for expansion in sparus_list:
        key = expansion.split(":")[0]
        value = expansion.split(":")[1].split("*")[0]
        expansions_dict[key] = value

    expansions_df = pd.DataFrame.from_dict(expansions_dict, orient='index')

    return hogs_expanded, expansions_df

def parse_eggNOG(clean_orthofinder_file, eggnog_file):
    import pandas as pd

    # N0_clean.tsv file from prepare_file.py script
    species_per_gene_family = clean_orthofinder_file


    df = pd.read_csv(species_per_gene_family, sep='\t') # Clean file from orthofinder, genes per species per HOG
    df = df.set_index("HOG")
    df_Saur = df["Sparus_aurata_KLR_BA"] #dataframe with Sparus aurata genes per HOG

    #Convert df to dictionary
    saur_dict = df_Saur.to_dict()
    for k,v in saur_dict.items():
        #k=k.split(".")[1] #Remove N0 from each HOG name
        if v.find(",")<0:
            saur_dict[k] = v
        saur_dict[k]=v.split(",")

    #Parse eggNOG mapper output file
    gene_name = []
    description = []
    with open(eggnog_file, "r") as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("#") or line.startswith("query"): #Keep only lines with gene descriptions
                continue
            line = line.rstrip().split("\t")
            gene_name.append(line[0]) # list with all the annotated genes
            if line[7]==['-']:
                description.append("-")
            description.append(line[7]) #corresponding list with description of gene's function


    description_dict = dict()

    for key,val in saur_dict.items():
        res=[]
        for genes in val:
            for ele in range(len(gene_name)):
                if gene_name[ele] in genes:
                    res.append(description[ele]) #gene_name[ele] and description[ele]: corresponding gene name and description

            description_dict[key] = res # Dictionary with HOGs as keys and descriptions of genes (as eggNOG indicated)
                                        # per HOG as values


    saur_dataframe = pd.Series(saur_dict).to_frame('Genes')
    description_dataframe = pd.Series(description_dict).to_frame('Description')
    total_df = pd.concat([saur_dataframe , description_dataframe], axis=1)

    #Rename the index column
    rows =[]
    for row in total_df.index:
        rows.append(row.split(".")[1])

    total_df.index = rows
    #total_csv = "/hog_with_GO.tsv" # File with eggNOG description for each gene per HOG
    #total_df.to_csv(total_csv, sep="\t")
    return total_df



#Call the above functions
path_to_speciesEventsCounts = '/path/to/speciesEventsCounts/'
df_generax,sparus_generax_dups = parse_generax()

path_to_fams_txt = "/path/to/CAFE_report/summary_cafe_fams.txt"
sparus_cafe_hogs, cafe_expansions_df = parse_cafe(path_to_fams_txt)

clean_orthofinder_file = '/path/to/N0_clean.tsv'
eggnog_file =  /path/to/out.emapper.annotations'
eggNOG_df = parse_eggNOG(clean_orthofinder_file, eggnog_file)

#dups_generax = df_generax["SparusaurataKLRBA"]

#Keep the intersection

intersection = list(set(sparus_cafe_hogs).intersection(set(sparus_generax_dups))) # Intersection of CAFE and GeneRax HOGs
intersection_GO = eggNOG_df.loc[intersection,:] # Keep the HOGs that both GeneRax and CAFE indicated
intersection_GO["#Genes"] = intersection_GO['Genes'].str.len() # Count number of genes per family (initially from Orthofinder)
sorted_intersection = intersection_GO.sort_values(by="#Genes", axis=0, ascending=False) # Descending order of df based in # of genes per HOG

#Add results from GeneRax and CAFE
add_generax_dups = pd.concat([sorted_intersection,df_generax.loc[intersection,:]["SparusaurataKLRBA"]],axis=1)
final_df = pd.concat([add_generax_dups,cafe_expansions_df.loc[intersection,:]],axis=1) #Add column of Cafe

final_df = final_df.rename(columns={"SparusaurataKLRBA": "GeneRax_Duplications", 0: "CAFE_expansions"})

#complete_csv = "/complete_GO_CAFE_GeneRAX.tsv" #File with the union of GeneRax and CAFE HOGs
#final_df.to_csv(complete_csv, sep="\t")

# Creation of a venn diagramm for the intersected gene families
venn2([set(sparus_cafe_hogs), set(sparus_generax_dups)],set_labels = ('CAFE', 'GeneRax'))

plt.savefig("/path/to/venn.png")
plt.show()
