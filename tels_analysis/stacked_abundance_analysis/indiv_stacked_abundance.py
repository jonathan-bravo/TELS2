from math import log10
import csv

class IndivStackedAbundance:

    def __init__(self, is_amr, mutable_category_list, mges_annot = dict()):
        # Object vars that are always used
        self.category_list = mutable_category_list
        self.absolute_abundance = dict()
        self.relative_abundance = dict()
        self.stats_filepaths = dict()
        self.is_amr = is_amr

        # Object vars used only in amr analysis
        self.group_to_class = dict()

        # Object vars used only in mge analysis
        self.mge_annot = mges_annot.copy()
        

    def add_to_absolute(self, legend, filepath, stats):
        # If this is the first time we encounter 
        # this sequencing platform and probe type
        if legend not in self.stats_filepaths:
            self.stats_filepaths[legend] = list()
            self.absolute_abundance[legend] = dict()

        self.stats_filepaths[legend].append(stats)
        with open(filepath, "r") as csv_file:
            csv_reader = csv.reader(csv_file)
            # Loop through diversity file
            for line_num, line in enumerate(csv_reader):
                # Skips information on total on-target reads and richness
                if line_num < (19 if self.is_amr else 20):
                    continue

                # AMR analysis
                if self.is_amr:
                    # Count alignments by ARG group
                    arg_header = line[0].split("|")
                    if arg_header[4] not in list(self.absolute_abundance[legend].keys()):
                        self.absolute_abundance[legend][arg_header[4]] = int(line[1])
                    else:
                        self.absolute_abundance[legend][arg_header[4]] += int(line[1])

                    # Save AMR class information if needed
                    if arg_header[4] not in self.group_to_class:
                        if arg_header[1] != "Drugs":
                            class_name = arg_header[1] + " resistance"
                        elif arg_header[2] == "betalactams":
                            class_name = "Betalactams"
                        elif arg_header[2] == "Mycobacterium_tuberculosis-specific_Drug":
                            class_name = "M_tuberculosis-specific_Drug"
                        else:
                            class_name = arg_header[2]

                        self.group_to_class[arg_header[4]] = class_name

                    if self.group_to_class[arg_header[4]] not in self.category_list:
                        self.category_list.append(self.group_to_class[arg_header[4]])

                # MGE analysis
                else:
                    # Count alignments by MGE accession
                    if line[0] not in list(self.absolute_abundance[legend].keys()):
                        self.absolute_abundance[legend][line[0]] = int(line[1])
                    else:
                        self.absolute_abundance[legend][line[0]] += int(line[1])

                    # Save MGE type information if needed
                    if self.mge_annot[line[0]] not in self.category_list:
                        self.category_list.append(self.mge_annot[line[0]])
    
    def make_abundance_relative(self):
        # Get average read count of replicate for each sample group
        average_read_count = dict()
        for legend, filepath_list in self.stats_filepaths.items():
            for filepath in filepath_list:
                with open(filepath, "r") as stat_file:
                    stat_file.readline()                    # skip duplicate stats information

                    if legend not in average_read_count:    
                        average_read_count[legend] = int(stat_file.readline().split(',')[1])
                    else:
                        average_read_count[legend] += int(stat_file.readline().split(',')[1])

            average_read_count[legend] = average_read_count[legend]/3

        # Make abundance relative to average read count and save
        # total relative abundance by gene for sorting
        relative_total = dict()
        relative_by_probe_type = dict()
        for legend, alignment_count_dict in self.absolute_abundance.items():

            # This is the first time we encounter 
            # this sequencing platform and probe type
            if legend not in relative_by_probe_type:
                relative_by_probe_type[legend] = dict()
                
            for gene, count in alignment_count_dict.items():

                relative_by_probe_type[legend].update(
                    {gene:log10(count/average_read_count[legend]*1000000)})
                
                if gene not in relative_total: 
                    relative_total[gene] = relative_by_probe_type[legend][gene]
                else:
                    relative_total[gene] += relative_by_probe_type[legend][gene]

        # Sorting relative abundance
        sorted_relative_total = {}
        while len(relative_total) > 0:
            max = ('key', -1)
            for key,val in relative_total.items():
                if val > max[1]:
                    max = (key,val)
            sorted_relative_total.update({max[0]:max[1]})
            relative_total.pop(max[0])

        # AMR analysis
        if self.is_amr:
            group_to_class_unsorted_copy = self.group_to_class
            self.group_to_class = dict()
            for key in sorted_relative_total:
                self.group_to_class[key] = group_to_class_unsorted_copy[key]
                for legend in relative_by_probe_type:
                    if legend not in self.relative_abundance:
                        self.relative_abundance[legend] = dict()
                    if key not in relative_by_probe_type[legend]:
                        self.relative_abundance[legend].update({key:0})
                    else:
                        self.relative_abundance[legend].update({key:relative_by_probe_type[legend][key]})
            group_to_class_unsorted_copy.clear()
            sorted_relative_total.clear()

        # MGE analysis
        else:
            mge_annot_unsorted_copy = self.mge_annot
            self.mge_annot = dict()
            for key in sorted_relative_total:
                self.mge_annot[key] = mge_annot_unsorted_copy[key]
                for legend in relative_by_probe_type:
                    if legend not in self.relative_abundance:
                        self.relative_abundance[legend] = dict()
                    if key not in relative_by_probe_type[legend]:
                        self.relative_abundance[legend].update({key:0})
                    else:
                        self.relative_abundance[legend].update({key:relative_by_probe_type[legend][key]})
            mge_annot_unsorted_copy.clear()
            sorted_relative_total.clear()
        

    def get_abundance(self):
        return self.relative_abundance

    def get_categories(self):
        if self.is_amr:
            return self.group_to_class
        else:
            return self.mge_annot
        