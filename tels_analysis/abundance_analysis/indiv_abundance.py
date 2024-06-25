from math import log10
import statistics

class IndivAbundance:

    def __init__(self, is_amr, mge_annot = dict()):
        # Object vars that are always used
        self.absolute_abundance = dict()
        self.relative_abundance = dict()
        self.is_amr = is_amr

        # Object vars used only in amr analysis
        self.accession_to_group = dict()
        
        # Object vars used only in mge analysis
        self.mge_annot = mge_annot
        

    def add_to_absolute(self, filepath):
        # Open file that has diversity information
        diversity_file = open(filepath, "r")

        # Loop through diversity file
        for line_num, line in enumerate(diversity_file):
            # Skips information on total on-target reads and richness
            if line_num < (19 if self.is_amr else 20):
                continue

            # Split line by comma-defined columns
            line_list = line.split(',')
            #print(line_list)
            if self.is_amr:
                accession = line_list[0].split("|")

                # If this is the first time we meet accession in sample group
                if line_list[0] not in self.absolute_abundance:
                    self.absolute_abundance[line_list[0]] = 0
                    self.accession_to_group[line_list[0]] = accession[4]
                
                # Add up abundance through all samples in sample group
                self.absolute_abundance[line_list[0]] += int(line_list[1])

            else:

                # If this is the first time we meet accession in sample group
                if line_list[0] not in self.absolute_abundance:
                    self.absolute_abundance[line_list[0]] = 0
                
                # Add up abundance through all samples in sample group
                self.absolute_abundance[line_list[0]] += int(line_list[1])

    def make_abundance_relative(self, file_size, gene_length):
        # In amr analysis, the sum of the intermediate values for each accession
        # in a group is logarithmized in order to find the relative abundance
        if self.is_amr:
            intermediate_values = dict()
            for accession, count in self.absolute_abundance.items():

                # If this is the first time we encounter arg group
                if self.accession_to_group[accession] not in intermediate_values:
                    intermediate_values[self.accession_to_group[accession]] = 0

                # Add abundance relative to gene length and sample group file size
                intermediate_values[self.accession_to_group[accession]] += (
                    (100.0 * float(count))/(float(gene_length[accession]*file_size)))
            for group, relative_abundance in intermediate_values.items():
                self.relative_abundance[group] = log10(relative_abundance)

        # # If mge analysis, take the logarithm for each type
        # else:
        #     intermediate_values = dict()
        #     for accession, count in self.absolute_abundance.items():

        #         # If this is the first time we encounter mge type
        #         if self.mge_annot[accession] not in intermediate_values:
        #             intermediate_values[self.mge_annot[accession]] = 0

        #         # Add abundance relative to gene length and sample group file size
        #         intermediate_values[self.mge_annot[accession]] += (
        #             (100.0 * float(count))/(float(gene_length[accession]*file_size))
        #         )
        #     for group, relative_abundance in intermediate_values.items():
        #         self.relative_abundance[group] = log10(relative_abundance)

        # If mge analysis, take the logarithm for each accession
        else:
            for accession, count in self.absolute_abundance.items():
                self.relative_abundance[accession] = log10(
                    (100.0 * float(count))/(float(gene_length[accession]*file_size))
                )
        

    def get_abundance(self):
        return self.relative_abundance