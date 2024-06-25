from Bio import SeqIO
import csv

# lambda that generates tels output file
tels_file_path = (lambda self, sample_name, extension : 
             self.source_prefix + sample_name + self.source_suffix + extension)


def fileDict(sample_name):
    sample = ""
    seqPlatform = ""

    #Determine organism
    if sample_name[0] == 'B':
        sample = "Bovine fecal"
        sample_name = sample_name[2:]
    elif sample_name[0] == 'H':
        sample = "Human fecal"
        sample_name = sample_name[2:]
    elif sample_name[0] == 'M':
        sample = "Mock"
        sample_name = sample_name[2:]
    else: #sample_name[0] = 'S'
        sample = "Soil"
        sample_name = sample_name[1:]

    #Determine whether V2 or XT
    if sample_name[0] == 'V':
        sample = sample + " (V2)"
        sample_name = sample_name[2:]
    else: #sample_name[0] = 'X'
        sample = sample + " (XT)"
        sample_name = sample_name[2:]

    #Determine probes
    if sample_name[0:2] == "AM":
        sample = sample + " + ARG-MGE probe"
        seqPlatform = "TELSeq"
        sample_name = sample_name[2:]
    elif sample_name[0] == 'A':
        sample = sample + " + ARG probe"
        seqPlatform = "TELSeq"
        sample_name = sample_name[1:]
    elif sample_name[0] == 'M':
        sample = sample + " + MGE probe"
        seqPlatform = "TELSeq"
        sample_name = sample_name[1:]
    else: #sample_name[0] == 'N'
        seqPlatform = "PacBio"
        return (sample, seqPlatform)

    #Determine baits
    sample = sample + sample_name
    return (sample, seqPlatform)

def sort(list_of_pairs):
    if len(list_of_pairs) == 1:
        return list_of_pairs
    else:
        listA = sort(list_of_pairs[0:int(len(list_of_pairs)/2)])
        listB = sort(list_of_pairs[int(len(list_of_pairs)/2):])
        sorted_list = list()
        for i in range(0,len(list_of_pairs)):
            if len(listA) == 0:
                sorted_list.append(listB.pop(0))
            elif len(listB) == 0:
                sorted_list.append(listA.pop(0))
            elif listA[0][0] < listB[0][0]:
                sorted_list.append(listA.pop(0))
            else:
                sorted_list.append(listB.pop(0))
        return sorted_list

def megares_analyzer(megares_annot_file):   
    drug_mechanisms = list()
    other_mechanisms = list()
    with open(megares_annot_file, "r") as csv_file:
        csv_reader = csv.reader(csv_file)
        for row_num, row in enumerate(csv_reader):
            if row_num == 0: continue
            class_mech = (row[2], row[3])
            if row[2] == "betalactams": 
                class_mech = ("Betalactams", row[3])
            elif row[2] == "Quaternary Ammonium Compounds (QACs) resistance":
                class_mech = ("QACs", row[3])
            elif row[2] == "Drug and biocide and metal resistance":
                class_mech = ("Drug biocide\nand metal", row[3])
            elif row[2] == "Cationic antimicrobial peptides":
                class_mech = ("Cationic peptides", row[3])
            elif row[2] == "Phenolic compound resistance":
                class_mech = ("Phenolic cpd.", row[3])
            elif 'resistance' in row[2]:
                class_mech = (row[2][:-11], row[3])
            if row[1] == "Drugs":
                if drug_mechanisms.count(class_mech) == 0:
                    drug_mechanisms.append(class_mech)
            else:
                if other_mechanisms.count(class_mech) == 0:
                    other_mechanisms.append(class_mech)
    drug_mechanisms = sort(drug_mechanisms)
    other_mechanisms = sort(other_mechanisms)
    return (drug_mechanisms, other_mechanisms)

def mge_analyzer(mge_annot_file):
    mge_annot_list = list()
    with open(mge_annot_file, "r") as csv_file:
        csv_reader = csv.reader(csv_file)
        for row_num, row in enumerate(csv_reader):
            if row_num == 0: continue
            type_accession = (row[1], row[0])
            if mge_annot_list.count(type_accession) == 0:
                mge_annot_list.append(type_accession)
    mge_annot_list = sort(mge_annot_list)
    return mge_annot_list

def get_genes_length(fasta_file):
    gene_length_dict = dict()
    with open(fasta_file, "r") as fasta:
        fasta_reader = SeqIO.parse(fasta, 'fasta')
        for reference in fasta_reader:
            gene_length_dict.update({reference.id: len(reference.seq)})
    return gene_length_dict

# Returns sample name definition in tuple form:
# (Organism, Platform, Chemistry, Probe)
def get_sample_name_definition(sample_name, new_probe_name=False):
    triplicate = sample_name[-1]
    
    # Determine organism
    if sample_name[0] == 'B':
        if new_probe_name:
            organism = "BF"
        else:
            organism = "Bovine"
        sample_name = sample_name[2:]
    elif sample_name[0] == 'H':
        if new_probe_name:
            organism = "FMT"
        else:
            organism = "Human"
        sample_name = sample_name[2:]
    elif sample_name[0] == 'M':
        if new_probe_name:
            organism = "MOCK"
        else:
            organism = "Mock"
        sample_name = sample_name[2:]
    else: #sample_name[0] = 'S'
        if new_probe_name:
            organism = "PPS"
        else:
            organism = "Soil"
        sample_name = sample_name[1:]

    # Determine chemistry
    chemistry = sample_name[:2]
    sample_name = sample_name[2:]
    
    # Determine probe specificity and
    # sequencing platform
    if sample_name[0:2] == "AM":
        platform = "TELSeq"
        probe = "Combo"
    elif sample_name[0] == 'A':
        platform = "TELSeq"
        if new_probe_name:
            probe = "RES"
        else:
            probe = "ARG"
    elif sample_name[0] == 'M':
        platform = "TELSeq"
        if new_probe_name:
            probe = "MOB"
        else:
            probe = "MGE"
    else: #sample_name[0] == 'N'
        platform = "PacBio"
        probe = None
    #print(f'{organism} {platform} {chemistry} {probe} {triplicate}')
    # return (organism, platform, chemistry, probe)
    return (organism, platform, chemistry, probe, triplicate)

def getSampleAndIndex(sample_name):
    sub_table = ""
    index = 0

    #Determine organism
    if sample_name[0] == 'B':
        sub_table = "Bovine"
        sample_name = sample_name[2:]
    elif sample_name[0] == 'H':
        sub_table = "Human"
        sample_name = sample_name[2:]
    elif sample_name[0] == 'M':
        sub_table = "Mock"
        sample_name = sample_name[2:]
    else: #sample_name[0] = 'S'
        sub_table = "Soil"
        sample_name = sample_name[1:]

    #Determine whether V2 or XT
    if sample_name[0] == 'V':
        sub_table = sub_table + "+V2"
        sample_name = sample_name[2:]
    else: #sample_name[0] = 'X'
        sub_table = sub_table + "+XT"
        sample_name = sample_name[2:]

    #Determine probes
    if sample_name[0:2] == "AM":
        index = 2
    elif sample_name[0] == 'A':
        index = 1
    elif sample_name[0] == 'M':
        index = 3
    else: #sample_name[0] == 'N'
        index = 4
    return (sub_table,index)

def get_mge_annot_dict(filepath):
    mge_annot = dict()
    with open(filepath, "r") as csv_file:
        csv_reader = csv.reader(csv_file)
        for line_num, line in enumerate(csv_reader):
            if line_num == 0: continue
            mge_annot[line[0]] = line[1]
    return mge_annot