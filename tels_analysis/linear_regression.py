from Bio import SeqIO
import csv

def megares_full_reference_length(megares_file_name, reference_length_file_name):
    total_length = 0
    for gene in SeqIO.parse(megares_file_name, 'fasta'):
        if "RequiresSNPConfirmation" in gene.name:
            continue
        total_length += len(gene.seq)
    with open(reference_length_file_name, 'w') as reference_length_file:
        reference_length_writer = csv.writer(reference_length_file)
        reference_length_writer.writerow(("Megares", str(total_length)))

def mge_full_reference_length(sample_name, SOURCE_PREFIX, SOURCE_SUFFIX, OVERLAP_OUTPUT, aclame_file_name, iceberg_file_name, plasmid_finder_file_name, reference_length_file_name, mge_lock):
    total_length = 0
    overlap_mges = list()
    overlap_file_name = SOURCE_PREFIX + sample_name + SOURCE_SUFFIX + OVERLAP_OUTPUT
    with open(overlap_file_name, 'r') as overlap_file:
        overlap_file_reader = csv.reader(overlap_file, delimiter=',')
        for row in overlap_file_reader:
            overlap_mges.append(row[0])
    for gene in SeqIO.parse(aclame_file_name, 'fasta'):
        if gene.name in overlap_mges:
            continue
        total_length += len(gene.seq)
    for gene in SeqIO.parse(iceberg_file_name, 'fasta'):
        if gene.name in overlap_mges:
            continue
        total_length += len(gene.seq)
    for gene in SeqIO.parse(plasmid_finder_file_name, 'fasta'):
        if gene.name in overlap_mges:
            continue
        total_length += len(gene.seq)
    mge_lock.acquire()
    with open(reference_length_file_name, 'a') as reference_length_file:
        reference_length_writer = csv.writer(reference_length_file)
        reference_length_writer.writerow((sample_name, str(total_length)))
    mge_lock.release()