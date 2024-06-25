from tels_analysis import fileDict
import numpy
import csv


class IndivComposition:
    def __init__(self, sample_name_def, is_amr_analysis, mge_dict = dict()):
        self.sample_name_definition = sample_name_def
        self.is_amr_analysis = is_amr_analysis

        if self.is_amr_analysis:
            self.element_types = [
                'Drugs', 'Metals', 'Biocides', 'Multi-compound']
        else:
            self.element_types = [
                'PLASMID', 'PROPHAGE', 'TE', 'IS', 'ICE', 'VIRUS', 'UNCLASSIFIED']

        # In order if amr analysis: drug, metal, biocide, and multi-compound
        # In order if mge analysis: plasmid, phage, TE, IS, ICE, virus, unclassified
        self.composition = [
            [0] * len(self.element_types),
            [0] * len(self.element_types),
            [0] * len(self.element_types)]

        # In order if amr analysis: class, mechanism, group
        # In order if mge analysis: type, accession
        self.richness = [
            [0] * (3 if self.is_amr_analysis else 2),
            [0] * (3 if self.is_amr_analysis else 2),
            [0] * (3 if self.is_amr_analysis else 2)]
        self.cr_index = 0

        self.mge_dict = mge_dict

    def get_richness(self):
        # Find median and standard dev across replicates 
        # for each annotation level in self.richness,
        # returns info as csv row
        summary = [0] * len(self.richness[0])
        for level in range(len(self.richness[0])):
            richness_for_level = [
                self.richness[0][level],
                self.richness[1][level],
                self.richness[2][level]
            ]
            median = numpy.median(richness_for_level)
            stdev = numpy.std(richness_for_level, ddof=1)
            summary[level] = str(median) + ' (' + "{0:.3f}".format(stdev) + ')'

        return ','.join(self.sample_name_definition) + ',' + ','.join(summary)
    
    def get_composition(self):
        # Return composition information for all three replicates
        return self.composition


    def add_to_data(self, filepath):
        # Keep track of all types, classes,
        # mechanisms, and groups encountered
        level_lists = list()
        for level in range(len(self.richness[0])):
            level_lists.append(list())

        with open(filepath, "r") as csv_file:
            csv_reader = csv.reader(csv_file)
            for row_num, row in enumerate(csv_reader):
                if row_num < (19 if self.is_amr_analysis else 20):
                    continue

                # If amr analysis, use accession header to calculate richness
                if self.is_amr_analysis:
                    arg_annotations = row[0].split('|')
                    type = arg_annotations[1]
                    for index in range(3):
                        if arg_annotations[2 + index] not in level_lists[index]:
                            level_lists[index].append(arg_annotations[2 + index])

                # If mge analysis, use mge dict and full
                # accession name to calculate richness
                else:
                    type = self.mge_dict[row[0]]
                    if type not in level_lists[0]:
                        level_lists[0].append(type)
                    if row[0] not in level_lists[1]:
                        level_lists[1].append(row[0])

                self.composition[self.cr_index][self.element_types.index(type)] += int(row[1])

        for level in range(len(level_lists)):
            self.richness[self.cr_index][level] = len(level_lists[level])

        self.cr_index += 1
  
