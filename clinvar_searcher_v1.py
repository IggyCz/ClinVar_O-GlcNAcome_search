import vcf
import csv
import requests
import json
from time import sleep
import os

counter = 0
#changes working directory
os.chdir(r"C:\Users\130002112\OneDrive - University of Dundee\AllmyFiles\Conor_project")

#Loads files necessary for analsysis, in order:
#vcf of clinvar variants
#Conor's list of O-GlcNAc modified proteins
#O-GlcNAcome published in PMID: 33479245

vcf_reader = vcf.Reader(open(r"clinvarCh38_25_12_2021.vcf", 'r'))
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
def variant_filter(variant):
    if None in variant:
        return False
    if "Pathogenic" in variant[2] or "Likely_pathogenic" in variant[2] or\
    "Uncertain_significance" in variant[2] or "Conflicting_interpretations_of_pathogenicity" in variant[2]:
        if 'SO:0001583|missense_variant' in variant[3]:
            return True    
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
        pc = variant_json["result"][uid]["protein_change"].split(",")[0]
        consequence_list.append([uid, pc])
    sleep(5)
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
    return [row[0],row[sites_index].split(separator)]



#takes a csv file of in which first column is the gene name, and site information
#is found in colum = site_index, should indicate the way multiple sites are split in one
#column
GlcNAcome_dict = {}
def GlcNAcome_dict_maker(csv_reader_object, separator = ",", site_index = 1):
    with open(csv_file, "r") as GlcNAcome:
        GlcNAcome = csv_reader_object
        GlcNAcome = map(lambda x: protein_sites(x, site_index, separator), list(GlcNAcome))
        for row in GlcNAcome:
            if row[0] in GlcNAcome_dict:
                for site in row[1]:
                    if site not in GlcNAcome_dict[row[0]]:
                        GlcNAcome_dict[row[0]].append(site[1:])
            else:
                GlcNAcome_dict[row[0]] = []
                for site in row[1]:
                    GlcNAcome_dict[row[0]].append(site[1:])

#Makes list of O-GlcNAc modified proteins (genes) from loaded O-GlcNAcome data
gene_list = [i[0] for i in GlcNAcome if i[0] != "Gene name(s)"]
for gene in OVS_GlcNAcome:
    if gene[0] not in gene_list and gene[4] != "":
        gene_list.append(gene[0])

#maps filter condition to vcf database, making a list of variants in genes in O-GlcNAcome
filter1_list = map(filter_vcf, vcf_reader)
filter1_list = [entry for entry in list(filter1_list) if entry != None]
print(len(filter1_list))

#makes new list with filtering old based on pathogenicity and variant being
#a missense mutation
filter2_list = [i[1] for i in list(filter(variant_filter, filter1_list))]
print(len(filter2_list))

#chops up filtered list to allow for submission to clinvar api
chopped_lists = [str(filter2_list[x:x+100]).replace(r"'","").replace("[","").replace("]","") for x in range(0, len(filter2_list), 100)]


pc_list = map(protein_consequence, chopped_lists)
pc_list = [item for sublist in list(pc_list) for item in sublist]

for index, item in enumerate(filter2_list):
    item.append(pc_list[index])

#to remove replaced with filter2_list
with open(r"filtered_patho.csv", "w", newline = "") as output:
    output = csv.writer(output)
    output.writerows(filtered_list)


GlcNAcome_dict_maker(GlcNAcome)
GlcNAcome_dict_maker(OVS_GlcNAcome)
                

#Make dictionar from lists 
variant_dict = {}
for row in filter2_list:
    clean_p_consequence = replace(row[4], r"'[]' ")
    variantID = clean_p_consequence.split(",")[0]
    protein_consequence = clean_p_consequence.split(",")[1]
    pathogenicity = row[2]
    if row[0] in variant_dict:
        variant_dict[row[0]][protein_consequence] = [variantID, pathogenicity]
    else:
        variant_dict[row[0]] = {protein_consequence: [variantID, pathogenicity]}

#saves to output csv
with open(r"output_clinvar.csv", "w", newline = "") as output_csv:
    output_csv = csv.writer(output_csv)
    output_csv.writerow(["Gene", "Protein Consequence", "Pathogenicity", "Variant ID"])
    for key in variant_dict:
        for protein_consequence in variant_dict[key]:
            if str(protein_consequence)[1:-1] != "" and str(protein_consequence)[1:-1] in GlcNAcome_dict[key]:
                row = [key, protein_consequence, variant_dict[key][protein_consequence][1], variant_dict[key][protein_consequence][0]]
                output_csv.writerow(row)