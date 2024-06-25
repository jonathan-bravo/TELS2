from bs4 import BeautifulSoup
from Bio import SeqIO
import requests
import csv
import sys

with open('annotations/TELS2_MGEs_Annotations.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    reader.__next__()
    class_and_db = dict()
    for row in reader:
        class_and_db[row[0]] = [row[1], row[2]]

gene_count = 0
already_found_organisms = dict()

# with open('databases/iceberg_db.fasta', 'r') as fastafile:
#     iceberg = dict()
#     for ref in SeqIO.parse(fastafile, 'fasta'):
#         if ref.id not in class_and_db:
#             continue
#         iceberg[str(ref.id)] = [ref.id.split('|')[2]]

#         accession = ref.id.split('|')[4]
#         URL = 'https://www.ncbi.nlm.nih.gov/nuccore/' + accession
#         r = requests.get(URL)
#         soup = BeautifulSoup(r.content, 'xml')
#         xml = soup.prettify().split('\n')
#         for line in xml:
#             if 'ORGANISM' in line:
#                 organism_accession = line[line.find('ORGANISM=')+len('ORGANISM='):]
#                 organism_accession = organism_accession[:organism_accession.find('&')]
#                 break
#         if organism_accession in already_found_organisms:
#             iceberg[str(ref.id)].append(already_found_organisms[organism_accession])
#         else:
#             URL = 'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=' + organism_accession
#             new_r = requests.get(URL)
#             new_soup = BeautifulSoup(new_r.content, 'xml')
#             new_xml = new_soup.prettify().split('\n')
#             i = -1
#             for new_line in new_xml:
#                 i += 1
#                 if 'Taxonomy ID: ' in new_line:
#                     organism_line = new_xml[i - 2]
#                     while organism_line[0] == ' ':
#                         organism_line = organism_line[1:]
#                     iceberg[str(ref.id)].append(organism_line)
#                     already_found_organisms[organism_accession] = organism_line
#                     break 
#                 elif 'Taxonomy browser' in new_line:
#                     organism_line = new_line[new_line.find('(')+1:new_line.find(')')]
#                     iceberg[str(ref.id)].append(organism_line)
#                     already_found_organisms[organism_accession] = organism_line
#                     break 
            
        
#         gene_count += 1
#         print(str(gene_count) + ' / ' + str(len(class_and_db)))
#         if organism_accession not in already_found_organisms:
#             print('still not fully working')
        
# with open('databases/plasmid_finder_db.fasta', 'r') as fastafile:
#     plasmid = dict()
#     for ref in SeqIO.parse(fastafile, 'fasta'):
#         if ref.id not in class_and_db:
#             continue
#         plasmid[str(ref.id)] = list()

# for gene in plasmid:
#     accession = gene[gene.rfind('_')+1:]
#     URL = 'https://www.ncbi.nlm.nih.gov/nuccore/' + accession
#     r = requests.get(URL)
#     soup = BeautifulSoup(r.content, 'xml')
#     xml = soup.prettify().split('\n')
#     i = 0
#     for line in xml:
#         i += 1
#         if '<div class="rprtheader">' in xml[i-3]:
#             plasmid[gene].append(line[11:(
#                 line.find(',') if ',' in ref.description
#                 else len(line))])
#             break

#     gene_count += 1
#     print(str(gene_count) + ' / ' + str(len(class_and_db)))

# with open('annotations/iceberg_annotations.csv', 'w') as csv_file:
#     for accession in iceberg:
#         csv_file.write(accession)
#         csv_file.write(',' + class_and_db[accession][0])
#         for annot in iceberg[accession]:
#             csv_file.write(',' + annot)
#         csv_file.write('\n')

# with open('annotations/plasmid_annotations.csv', 'w') as csv_file:
#     for accession in plasmid:
#         csv_file.write(accession)
#         csv_file.write(',' + class_and_db[accession][0])
#         for annot in plasmid[accession]:
#             csv_file.write(',' + annot)
#         csv_file.write('\n')

with open('databases/aclame_db.fasta', 'r') as fastafile:
    aclame = dict()
    for ref in SeqIO.parse(fastafile, 'fasta'):
        if ref.id not in class_and_db:
            continue
        aclame[str(ref.id)] = list()
        if 'MgeName: ' in ref.description:
            group = ref.description.split(': ')[5]
            aclame[ref.id].append(group[:group.find('#')-1])
        elif 'genbank:GeneID:' in ref.description:
            ID = ref.description.split(':')[-1]
            if ID[-1] == '\n':
                ID = ID[:-1]
            URL = 'https://www.ncbi.nlm.nih.gov/gene/' + ID
            r = requests.get(URL)
            soup = BeautifulSoup(r.content, 'xml')
            xml = soup.prettify().split('\n')
            i = 0
            inDiv = False
            for line in xml:
                i += 1
                if 'id="summaryDiv"' in line:
                    inDiv = True
                if inDiv:
                    if "</div>" in line:
                        inDiv = False
                    if 'Organism' in xml[i-5]:
                        aclame[ref.id].append(line)
                        break
        gene_count += 1
        print(str(gene_count) + ' / ' + str(len(class_and_db)))

with open('annotations/aclame_annotations.csv', 'w') as csv_file:
    for accession in aclame:
        csv_file.write(accession)
        csv_file.write(',' + class_and_db[accession][0])
        for annot in aclame[accession]:
            csv_file.write(',' + annot)
        csv_file.write('\n')