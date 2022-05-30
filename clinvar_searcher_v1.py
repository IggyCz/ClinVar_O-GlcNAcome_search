import vcf
import csv
import requests
import json
from time import sleep
import os

counter = 0
#Loads files necessary for analsysis, in order:
#Conor's list of O-GlcNAc modified proteins
#O-GlcNAcome published in PMID: 33479245
GlcNAcome = csv.reader(open(r"Conor_L_list.csv", "r"))
OVS_GlcNAcome = csv.reader(open(r"GlcNAcome.csv", "r"))

#Function to retrieve vcf entries from loaded file if in O-GlcNAcome gene list
#takes individual entries from vcf database, returns gene, variant ID, clinical significance
#and molecular consequence as list
def filter_vcf(record):
    summary = dict(record.INFO)
    if summary.get("GENEINFO"):
        end = summary.get("GENEINFO").find(":")
        if summary.get("GENEINFO")[:end] in gene_list:
            return([summary.get("GENEINFO")[:end], record.ID, summary.get("CLNSIG"), summary.get("MC")])


#function which can be used in filter to retain entries fulfiling criteria
def missense_filter(variant):
    if None in variant:
        return False
    elif 'SO:0001583|missense_variant' in variant[3]:
        return True
    else:
        return False


#api interactor takes list of variant IDs, returns list of lists of [variantID, protein consequence]
def protein_consequence(variants):
    global counter
    while True:
        api_response = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id="+variants+"&retmode=json")
        if str(api_response) == "<Response [200]>":
            api_response = api_response.text
            break
        else:
            print("api error")
            sleep(5)
            continue
    counter += 1
    print(counter)
    variant_json = json.loads(api_response)
    consequence_list = []
    for uid in variant_json["result"]["uids"]:
        try:
            pc = variant_json["result"][uid]["protein_change"].split(",")[0]
            consequence_list.append([uid, pc])
        except KeyError:
            with open("error_log.txt", "a") as errorLog:
                consequence_list.append([uid, "0"])
                errorLog.write(uid +"\n")
    sleep(0.5)
    return consequence_list


#Takes string and list of characters to remove from string, returns string without those
#characters
def replace(string, character_to_remove):
    new_string = string
    for character in character_to_remove:
        new_string = new_string.replace(character, "")
    return new_string

#Takes list including gene name and site information (as string e.g. with comma separated values)
#returns list [gene, [site, site, site]]
def protein_sites(row, sites_index, separator):
    gene = row[0].split(" ")[0]
    if gene != "Gene":
        return [gene ,row[sites_index].split(separator)]
    else:
        return [0]



#takes a csv file of in which first column is the gene name, and site information
#is found in column = site_index, should indicate the way multiple sites are split in one
#column. Appends global variable GlcNAcome_dict
def GlcNAcome_dict_maker(csv_reader_object, separator = " ", site_index = 1):
    global GlcNAcome_dict
    GlcNAcome = map(lambda x: protein_sites(x, site_index, separator), list(csv_reader_object))
    for row in GlcNAcome:
        if row[0] and row[0] in GlcNAcome_dict:
            for site in row[1]:
                if site not in GlcNAcome_dict[row[0]]["sites"] and site != "":
                    GlcNAcome_dict[row[0]]["sites"].append(site)
        elif row[0]:
            GlcNAcome_dict[row[0]] = {"sites": [], "mutations":{}}
            for site in row[1]:
                if site not in GlcNAcome_dict[row[0]]["sites"] and site != "":
                    GlcNAcome_dict[row[0]]["sites"].append(site)
        else:
            continue

#Filter function to retain pathogenic mutations
def pathogenicityCheck(mutation):
    for i in mutation[4]:
        if "Pathogenic" in i or "pathogenic" in i:
            return True
    return False

    
#Function which takes entries in dictionaries generated here and reformats them (and if
#specified filters based on pathogenicity) for saving as a csv 
def reformatForCsv(key, relativeSite):
    mutations = GlcNAcome_dict[key]["mutations"].keys()
    if not mutations:
        return False
    sites = [int(i[1:]) for i in GlcNAcome_dict[key]["sites"]]
    missense = [mutation for mutation in mutations if int(mutation[1:-1])-relativeSite in sites]
    pathogenicity = [[key, int(mutation[1:-1])-relativeSite, relativeSite, mutation, \
                      GlcNAcome_dict[key]["mutations"][mutation]["pathogenicity"], \
                      GlcNAcome_dict[key]["mutations"][mutation]["variantID"]] for mutation in missense]
    if relativeSite:
        pathogenicity = list(filter(pathogenicityCheck, pathogenicity))
        return pathogenicity
    else:
        return pathogenicity


#As intermediate .jsons are saved, this conditional allows for avoiding time consuming step if this file exists
if not os.path.exists('temp_var.json'):
    GlcNAcome = csv.reader(open(r"Conor_L_list.csv", "r"))
    OVS_GlcNAcome = csv.reader(open(r"GlcNAcome.csv", "r"))
    vcf_reader = vcf.Reader(open(r"clinvarCh38_25_12_2021.vcf", 'r'))
    #Makes list of O-GlcNAc modified proteins (genes) from loaded O-GlcNAcome data
    gene_list = [i[0] for i in GlcNAcome if i[0] != "Gene"]
    for gene in OVS_GlcNAcome:
        if gene[0] not in gene_list and gene[4] != "":
            gene_list.append(gene[0])

    #maps filter condition to vcf database, making a list of variants in genes in O-GlcNAcome
    filter1_list = map(filter_vcf, vcf_reader)
    filter1_list = [entry for entry in list(filter1_list) if entry != None]
    print(len(filter1_list))

    #makes new list with filtering old based on variant being
    #a missense mutation
    filter2_list = list(filter(missense_filter, filter1_list))
    print(len(filter2_list))
    temp_filter2_list = [i[1] for i in filter2_list]

    #chops up filtered list to allow for submission to clinvar api in groups of 100
    chopped_lists = [str(temp_filter2_list[x:x+100]).replace(r"'","").replace("[","").replace("]","") for x in range(0, len(temp_filter2_list), 100)]

    pc_list = map(protein_consequence, chopped_lists)
    pc_list = [item for sublist in list(pc_list) for item in sublist]

    for index, item in enumerate(filter2_list):
        item.append(pc_list[index])
        
    #Make dictionary from lists 
    variant_dict = {}
    for row in filter2_list:
        clean_p_consequence = row[4]
        variantID = clean_p_consequence[0]
        protein_consequence = clean_p_consequence[1]
        pathogenicity = row[2]
        if row[0] in variant_dict:
            variant_dict[row[0]][protein_consequence] = [variantID, pathogenicity]
        else:
            variant_dict[row[0]] = {protein_consequence: [variantID, pathogenicity]}
    
    with open("temp_var.json", "w") as temp_var:
        json.dump(variant_dict, temp_var)
else:
    variant_dict = json.load(open("temp_var.json", "r"))

#as previous conditional
if not os.path.exists('temp_Glc_dict.json'):
    GlcNAcome = csv.reader(open(r"Conor_L_list.csv", "r"))
    OVS_GlcNAcome = csv.reader(open(r"GlcNAcome.csv", "r"))
    GlcNAcome_dict = {}
    GlcNAcome_dict_maker(GlcNAcome)
    GlcNAcome_dict_maker(OVS_GlcNAcome, "/", 4)
    with open("temp_Glc_dict.json", "w") as temp_Glc_dict:
        json.dump(GlcNAcome_dict, temp_Glc_dict)      

else:
    GlcNAcome_dict = json.load(open("temp_Glc_dict.json", "r"))





tempAllProximal = []
#Creates a list of missense mutations +/-5 amino acids
#relative to O-GlcNAc sites, based on O-GlcNAcome dictionary
#and variant dictionary
for key in variant_dict:
    for protein_consequence in variant_dict[key]:
        for GlcSite in GlcNAcome_dict[key]["sites"]:
            GlcSite = GlcSite[1:]
            try:
                if str(protein_consequence)[1:-1] != "" and int(protein_consequence[1:-1]) <= int(GlcSite) + 5\
                   and int(protein_consequence[1:-1]) >= int(GlcSite) - 5:
                    row = [key, GlcSite, protein_consequence, variant_dict[key][protein_consequence][1], variant_dict[key][protein_consequence][0]]
                    tempAllProximal.append(row)
            except ValueError:
                     continue

#Places missense mutations from tempAllProximal list
#and adds them to the GlcNAcome dictionary
for row in tempAllProximal:
    if row[2] in GlcNAcome_dict[row[0]]["mutations"]:
        if row[4] in GlcNAcome_dict[row[0]]["mutations"][row[2]]["variantID"]:
            continue
        else:
            GlcNAcome_dict[row[0]]["mutations"]["pathogenicty"].append(row[3])
            GlcNAcome_dict[row[0]]["mutations"]["variantID"].append(row[4])
    else:
        GlcNAcome_dict[row[0]]["mutations"][row[2]] = {"pathogenicity":[row[3]], "variantID": [row[4]]}

with open("temp_Glc_dict.json", "w") as temp_Glc_dict:
        json.dump(GlcNAcome_dict, temp_Glc_dict)
    
genes = list(GlcNAcome_dict.keys())

#writes header for output csv
with open("output_v3.csv", "w") as output_:
        output_ = csv.writer(output_)
        header = ["Gene","Site", "Relative position of mutation", "Mutation", "Pathogenicity", "VariantID"]
        output_.writerow(header)
        
#writes final list of mutations to csv
for relativeToSite in range(-5,6):
    rows = [i for i in map(lambda x: reformatForCsv(x, relativeToSite), genes) if i]
    flat_list = [item for sublist in rows for item in sublist]
    with open("output_v3.csv", "a", newline = "") as output_:
        output_ = csv.writer(output_)
        output_.writerows(flat_list)




