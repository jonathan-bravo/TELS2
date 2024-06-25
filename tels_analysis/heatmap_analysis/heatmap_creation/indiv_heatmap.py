import numpy
import csv

class IndivHeatmap:
    def __init__(self, organism, indiv_point_dict, antimicrobial = None):
        self.organism = organism
        self.antimicrobial = antimicrobial
        self.indiv_point_dict = indiv_point_dict
        self.indiv_point_list = [k for k in self.indiv_point_dict] # THIS IS THE VALUE THAT WILL MAKE THE RIGHT LABEL AXIS

        #print(self.indiv_point_dict)

        self.columns = [
            [0]*len(self.indiv_point_dict), [0]*len(self.indiv_point_dict), [0]*len(self.indiv_point_dict),  # V2 RES
            [0]*len(self.indiv_point_dict), [0]*len(self.indiv_point_dict), [0]*len(self.indiv_point_dict),  # V2 MOB
            [0]*len(self.indiv_point_dict), [0]*len(self.indiv_point_dict), [0]*len(self.indiv_point_dict),  # V2 Combo
            [0]*len(self.indiv_point_dict), [0]*len(self.indiv_point_dict), [0]*len(self.indiv_point_dict),  # XT RES
            [0]*len(self.indiv_point_dict), [0]*len(self.indiv_point_dict), [0]*len(self.indiv_point_dict),  # XT MOB
            [0]*len(self.indiv_point_dict), [0]*len(self.indiv_point_dict), [0]*len(self.indiv_point_dict),  # XT Combo
            [0]*len(self.indiv_point_dict), [0]*len(self.indiv_point_dict), [0]*len(self.indiv_point_dict),  # PacBio
            [0]*len(self.indiv_point_dict), [0]*len(self.indiv_point_dict), [0]*len(self.indiv_point_dict)]  # PacBio
        
        self.x_axis_list = [
            'V2 RES',   'V2 MOB',   'V2 Combo',
            'XT RES',   'XT MOB',   'XT Combo',
            'PacBio']
        # self.x_axis_list = [
        #     'XT-HS2 RES',   'XT-HS2 MOB', 'XT-HS2 Combo',
        #     'XT RES',   'XT MOB',   'XT Combo',
        #     'PacBio']
        
        self.dupl_list = ['A', 'B', 'C', 'D', 'E', 'F']
        
    def add_to_map(
            self, x_axis, filepath, duplicate, bool_dict, mge_annotations = dict()):
        with open(filepath, "r") as csv_file:
            csv_reader = csv.reader(csv_file)
            for row_num, row in enumerate(csv_reader):
                if row_num < 19: continue
                if row[0] == 'MGE Header': continue

                # AMR analysis
                if self.antimicrobial != None:
                    arg_annot = row[0].split('|')
                    if self.antimicrobial == "Drugs":
                        if arg_annot[1] != "Drugs":
                            continue
                    else:
                        if arg_annot[1] == "Drugs":
                            continue
                    arg_mechansism = arg_annot[3].replace('_', ' ')

                    index = self.indiv_point_list.index(arg_mechansism)
                    column = 3 * self.x_axis_list.index(x_axis) + self.dupl_list.index(duplicate)
                    self.columns[column][index] = 1
                    arg_class = arg_annot[2].replace('_', ' ')
                    if arg_class == "betalactams":
                        arg_class = "Betalactams"
                    if arg_class == "Quaternary Ammonium Compounds (QACs) resistance":
                        arg_class = "QACs"
                    if arg_class == "Drug and biocide and metal resistance":
                        arg_class = "Drug biocide\nand metal"
                    if arg_class == "Cationic antimicrobial peptides":
                        arg_class = "Cationic peptides"
                    if arg_class == "Phenolic compound resistance":
                        arg_class = "Phenolic cpd."
                    if 'resistance' in arg_class:
                        arg_class = arg_class[:-11]
                    bool_dict[(arg_class, arg_mechansism)] = True

                # MGE analysis
                else:
                    index = self.indiv_point_list.index(row[0])
                    column = 3 * self.x_axis_list.index(x_axis) + self.dupl_list.index(duplicate)
                    self.columns[column][index] = 1
                    bool_dict[(mge_annotations[row[0]], row[0])] = True

#    ^
#   / \
#  /   \
# /     \
# |  |  |
# |  |  |
# \  \  /
#  \  \/
#  /\  \
# /  \  \
# |  |  |
# |  |  |
# \     /
#  \   /
#   \ /
#    v
    
    def make_map(self, bool_dict, category_list):
        # Remove ARG mechansisms or MGE accessions
        # that are in none of the samples
        #print(category_list)
        #print(bool_dict)
        right_label = []
        #print(len(bool_dict), len(self.columns[0]))
        for tup in bool_dict:
            if not(bool_dict[tup]):
                index = self.indiv_point_list.index(tup[1])
                #print(index)
                for column in self.columns:
                    #print(column)
                    column.pop(index)
                self.indiv_point_list.pop(index)
            else:
                right_label.append(tup)

        # Assign heatmap color value
        for column in self.columns:
            for row_num in range(0, len(column)):
                if column[row_num] == 1:
                    color_value = 1 + category_list.index(
                        self.indiv_point_dict[self.indiv_point_list[row_num]])
                    column[row_num] = color_value
        vals = numpy.array(self.columns).transpose()
        
        #print(len(bool_dict), len(self.columns[0]))

        # for t in bool_dict:
        #     print(t)
        #print(right_label)
        #print(len(right_label))

        return vals, self.x_axis_list, right_label