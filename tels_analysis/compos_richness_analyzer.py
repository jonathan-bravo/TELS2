from tels_analysis.compos_richness_analysis.indiv_composition import IndivComposition
from tels_analysis import get_sample_name_definition
from tels_analysis import get_mge_annot_dict
from tels_analysis import tels_file_path
from matplotlib import pyplot
from PIL import Image
import seaborn
import numpy
import os


class ComposRichnessAnalyzer:

    def __init__(
            self, SOURCE_PREFIX, SOURCE_SUFFIX, AMR_DIV_EXT,
            MGE_EXT, MGE_CLASSIFICATION, AMR_ANALYSIS, MGE_ANALYSIS):
        
        # Tels filepaths variables
        self.source_prefix = SOURCE_PREFIX
        self.source_suffix = SOURCE_SUFFIX
        self.amr_reads_ext = AMR_DIV_EXT
        self.mge_reads_ext = MGE_EXT

        # Variables for AMR analysis
        self.amr_analysis = AMR_ANALYSIS
        self.amr_cr = dict()
        
        # Variables for MGE analysis
        self.mge_dict = get_mge_annot_dict(MGE_CLASSIFICATION)
        self.mge_analysis = MGE_ANALYSIS
        self.mge_cr = dict()
        

    def analyze_sample(self, sample_name):
        # Find sample name definition
        sample_name_def = get_sample_name_definition(sample_name)
        if sample_name_def[1] == 'PacBio':
            sample_name_def = (
                sample_name_def[0],
                sample_name_def[1] + (' A' if sample_name_def[2] == 'V2' else ' B'),
                '__',
                '__'
            )
        
        # Remove replicate info from sample name
        new_sample_name = (sample_name[:-1]
                           if 'NEGAM' not in sample_name
                           else sample_name[:-2])

        # Go through amr analysis
        if self.amr_analysis:
            if new_sample_name not in self.amr_cr:
                self.amr_cr[new_sample_name] = IndivComposition(sample_name_def, True)
            self.amr_cr[new_sample_name].add_to_data(
                tels_file_path(self, sample_name, self.amr_reads_ext))

        # Go through mge analysis
        if self.mge_analysis:
            if new_sample_name not in self.mge_cr:
                self.mge_cr[new_sample_name] = IndivComposition(
                    sample_name_def, False, self.mge_dict)
            self.mge_cr[new_sample_name].add_to_data(
                tels_file_path(self, sample_name, self.mge_reads_ext))

    def make_bar_charts(
            self, ext, 
            arg_output_folder = None,
            mge_output_folder = None):
        
        # Find prime numbers for factorization
        prime_numbers = [2, 3]
        for num in range(1, 19134):
            is_prime = True
            for prime_num in prime_numbers:
                if (6*num-1) % prime_num == 0:
                    is_prime = False
                    break
                if (prime_num*prime_num) > (6*num-1):
                    break
            if is_prime:
                prime_numbers.append((6*num-1))

            is_prime = True
            for prime_num in prime_numbers:
                if (6*num+1) % prime_num == 0:
                    is_prime = False
                    break
                if (prime_num*prime_num) > (6*num+1):
                    break
            if is_prime:
                prime_numbers.append((6*num+1))

        # Factorize number, returning a list of all prime factors
        def factorization(num):
            factors = list()
            remainder = num
            prime_numbers_index = 0
            while (remainder > pow(prime_numbers[prime_numbers_index], 2)):
                if remainder % prime_numbers[prime_numbers_index] == 0:
                    factors.append(prime_numbers[prime_numbers_index])
                    remainder = int(remainder / prime_numbers[prime_numbers_index])
                else:
                    prime_numbers_index += 1
            factors.append(remainder)
            return factors

        # Find least common multiple of composition sum across replicates
        def least_common_multiple(sum_list):
            # Retrieve non-zero and non-duplicate sums
            all_nums = []
            for value in sum_list:
                if value != 0 and value not in all_nums:
                    all_nums.append(value)

            # If only one unique value, LCM is self
            if len(all_nums) == 1:
                return all_nums[0]

            # Find LCM
            prev_LCM = all_nums[0]
            for num in all_nums[1:]:
                prev_factorization = factorization(prev_LCM)
                curr_factorization = factorization(num)
                greatest_common_denominator = 1
                next_curr_index = 0
                for prev_factor in prev_factorization:
                    while (next_curr_index < len(curr_factorization)
                           and curr_factorization[next_curr_index] < prev_factor):
                        next_curr_index += 1
                    if next_curr_index == len(curr_factorization):
                        break
                    if prev_factor == curr_factorization[next_curr_index]:
                        greatest_common_denominator *= prev_factor
                        next_curr_index += 1
                prev_LCM = prev_LCM * int(num / greatest_common_denominator)
            return prev_LCM

        def bar_charts_per_analysis(
                output_folder, comp_data, legend, element_name):

            if not(os.path.exists(output_folder)):
                os.makedirs(output_folder)

            for index in range(0, len(comp_data)):
                # Find composition sum per replicate
                comp_sum = [sum(comp_data[index][1][0]),
                            sum(comp_data[index][1][1]),
                            sum(comp_data[index][1][2])]

                filepath = (output_folder + "/" + comp_data[index][0] + ext)

                if sum(comp_sum) != 0:
                    height = least_common_multiple(comp_sum)
                else:
                    height = 1

                # Get length of each type for each replicate
                new_list = list()
                for row_num in range(len(legend)):
                    row = list()
                    for sample in range(3):
                        if comp_sum[sample] == 0:
                            row.append(0)
                        else:
                            proportion = int(height / comp_sum[sample])
                            row.append(
                                comp_data[index][1][sample][row_num]*proportion)
                    new_list.append(row)

                # Get color map ready
                if len(legend) == 4:
                    color_map = seaborn.color_palette("colorblind", n_colors=4)
                    multiplier = 1
                    first_index = 0
                elif len(legend) == 7:
                    color_map = seaborn.color_palette("colorblind", n_colors=13)
                    multiplier = -1
                    first_index = 12

                bottom_array = [0,0,0]
                pyplot.figure(figsize=[16,30])
                for row, c_index in zip(new_list, range(len(legend))):
                    curr_array = numpy.array(row)
                    pyplot.bar(['A','B','C'],
                                width=0.9,
                                height=curr_array,
                                bottom=bottom_array,
                                color=color_map[first_index + multiplier * c_index])
                    bottom_array = curr_array + bottom_array

                # pyplot.legend(legend,
                #               bbox_to_anchor=(0.5, -0.03),
                #               loc='lower center',
                #               ncol=len(legend),
                #               reverse=True,
                #               fontsize=20)
                # pyplot.title(comp_data[index][0] + " " + element_name, fontsize=20)
                pyplot.ylim(0,height)
                pyplot.xticks([])
                pyplot.yticks([])
                pyplot.tight_layout()
                pyplot.savefig(filepath)
                pyplot.close()

        if self.amr_analysis:
            comp_data = list()
            for c_r in self.amr_cr:
                comp_data.append([c_r, self.amr_cr[c_r].get_composition()])
            legend = ['Drugs', 'Metal', 'Biocide', 'Multi-Compound']

            bar_charts_per_analysis(
                arg_output_folder, comp_data, legend, 'ARG'
            )

        if self.mge_analysis:
            comp_data = list()
            for c_r in self.mge_cr:
                comp_data.append([c_r, self.mge_cr[c_r].get_composition()])
            legend = ['Plasmid', 'Prophage', 'TE', 'IS', 'ICE', 'Virus', 'Unclassified']

            bar_charts_per_analysis(
                mge_output_folder, comp_data, legend, 'MGE'
            )

    def write_richness(self, output_folder, richness_ext):
        if not(os.path.exists(output_folder)):
                os.makedirs(output_folder)
        if self.amr_analysis:
            with open(output_folder + 'ARG' + richness_ext, "w") as output_file:
                # Write header row
                output_file.write(
                    "Sample,Seq. platform,Enrichment platform,"
                    + "Probe type,ARG class richness Med (± sd),"
                    + "ARG mechanism richness Med (± sd),"
                    + "ARG group richness Med (± sd)\n")
                for c_r in self.amr_cr:
                    output_file.write(self.amr_cr[c_r].get_richness() + '\n')

        if self.mge_analysis:
            with open(output_folder + 'MGE' + richness_ext, "w") as output_file:
                # Write header row
                output_file.write(
                    "Sample,Seq. platform,Enrichment platform,"
                    + "Probe type,MGE type richness Med (± sd),"
                    + "MGE accession richness Med (± sd)\n")
                for c_r in self.mge_cr:
                    output_file.write(self.mge_cr[c_r].get_richness() + '\n')

