from tels_analysis.stacked_abundance_analysis.indiv_stacked_abundance import IndivStackedAbundance
from tels_analysis import get_sample_name_definition
from tels_analysis import get_mge_annot_dict
from tels_analysis import tels_file_path
from matplotlib import pyplot
import seaborn
import numpy
import os

class StackedAbundanceAnalyzer:
    
    def __init__(
            self, SOURCE_PREFIX, SOURCE_SUFFIX, AMR_DIV_EXT, 
            MGE_EXT, STATS_EXT, MGES_ANNOTATION):
        
        # Get tels output file names information
        self.source_prefix = SOURCE_PREFIX
        self.source_suffix = SOURCE_SUFFIX
        self.amr_reads_ext = AMR_DIV_EXT
        self.mge_reads_ext = MGE_EXT
        self.stats_ext = STATS_EXT

        # Get MGE types
        self.mges_annot = get_mge_annot_dict(MGES_ANNOTATION)

        # Estabilsh dictionaries for absolute abundance 
        self.abundance_dict_amr = dict()
        self.abundance_dict_mge = dict()

        # Make mutable category list of mge types and arg classes encountered
        self.amr_category_list = list()
        self.mge_category_list = list()

    def find_absolute_abundance(self, sample_name, amr_analysis, mge_analysis):

        amr_filepath = tels_file_path(self, sample_name, self.amr_reads_ext)
        mge_filepath = tels_file_path(self, sample_name, self.mge_reads_ext)
        stats_filepath = tels_file_path(self, sample_name, self.stats_ext)

        # Retrieve definition in tuple form Organism, Platform, Chemistry, Probe)
        # which must be used to declare the subtables (Organism + Chemistry) 
        # and the different groups in the legend (Platform + Probe)
        sample_name_definition = get_sample_name_definition(sample_name)
        subtable = sample_name_definition[0] + ' ' + sample_name_definition[2]
        legend = (sample_name_definition[1] if sample_name_definition[1] == 'PacBio'
                  else sample_name_definition[1] + ' ' + sample_name_definition[3])

        # If we are analyzing AMR
        if amr_analysis:

            # This is the first time we see this 
            # organism + chemistry combination for amr
            if subtable not in self.abundance_dict_amr:
                self.abundance_dict_amr[subtable] = IndivStackedAbundance(
                    True, self.amr_category_list)

            self.abundance_dict_amr[subtable].add_to_absolute(
                legend, amr_filepath, stats_filepath)

        # If we are analyzing MGEs
        if mge_analysis: 

        # This is the first time we see this 
        # organism + chemistry combination for mge
            if subtable not in self.abundance_dict_mge:
                self.abundance_dict_mge[subtable] = IndivStackedAbundance(
                    False, self.mge_category_list, self.mges_annot)
                
            self.abundance_dict_mge[subtable].add_to_absolute(
                legend, mge_filepath, stats_filepath)
                
    def make_stacked_barplot(
            self, output_folder, stacked_ext, amr_analysis, mge_analysis):
        
        def superplot_per_analysis(abundance_dict, category_list, element_name):

            # Go through absolute abundace information to make it relative
            for subtable in abundance_dict:
                abundance_dict[subtable].make_abundance_relative()

            # Set up main figure
            fig = pyplot.figure(figsize=(80, 27),layout='constrained')
            subfigs = fig.subfigures(
                nrows=2,
                height_ratios=[21,6]
            )

            axs = subfigs[0].subplots(
                nrows=2,
                ncols=8,
                gridspec_kw={'height_ratios': [20,1]})

            # Go through each subtable
            for index, subtable in enumerate(abundance_dict):

                # Stacked bar plot
                pyplot.sca(axs[0][index])
                color_list = seaborn.color_palette("colorblind", n_colors=4)
                sub_abundance = abundance_dict[subtable].get_abundance()
                bottom_array = None
                for l_index, legend in enumerate(sub_abundance):
                    current_array = numpy.array(list(sub_abundance[legend].values()))
                    x_coords = list(sub_abundance[legend].keys())
                    if bottom_array is None:
                        axs[0][index].bar(
                            x_coords, current_array, width=1.0,label=legend,
                            color=color_list[l_index], alpha=0.5)
                        bottom_array = current_array
                    else:
                        axs[0][index].bar(
                            x_coords, current_array, width=1.0, label=legend,
                            color=color_list[l_index], alpha=0.5, bottom=bottom_array)
                        bottom_array = bottom_array + current_array
                if index == 0:
                    axs[0][index].set_ylabel('Log Relative Abundance', size=60)
                    pyplot.yticks(fontsize=45)
                else:
                    pyplot.tick_params('y', labelleft=False)
                if index !=  0:
                    axs[0][index].sharey(axs[0][0])
                pyplot.xticks([])
                pyplot.margins(x=0)
                axs[0][index].set_title(subtable, size=70)
                axs[0][index].set_anchor('NE')
                if index == 7: 
                    axs[0][index].legend(prop={'size': 50})   # only show legend color at upper right subplot

                # Color-coded x-axis: 
                #   assign number to each gene based on category alphabetical sorting
                #   and use this number to determine color shown on x_axis
                pyplot.sca(axs[1][index])
                gene_to_category = abundance_dict[subtable].get_categories()
                x_matrix = list()
                for arg in gene_to_category:
                    x_matrix.append(category_list.index(gene_to_category[arg]))
                numpy_array = numpy.array([x_matrix])
                seaborn.heatmap(numpy_array, ax = axs[1][index], 
                                xticklabels=False, yticklabels=False, cbar=False, 
                                cmap='viridis')
                axs[1][index].set_anchor('NE')

            # Legend for color-coded x-axis
            legend_axs = subfigs[1].subplots(
                nrows=1,
                ncols=5,
                gridspec_kw={'width_ratios': [1,12,12,12,12]})
            label_matrix = [*range(len(category_list))]
            increment_val = (len(category_list) + 3) // 4
            for ax_index, first_label_index in enumerate(
                    range(0, len(category_list), increment_val),1):
                pyplot.sca(legend_axs[ax_index])
                legend_axs[ax_index].yaxis.tick_right()
                if ax_index < 4:
                    last_label_index = first_label_index+increment_val
                else:
                    last_label_index = len(category_list)
                numpy_array = numpy.array(label_matrix[first_label_index:last_label_index])
                numpy_array = numpy_array.reshape(last_label_index-first_label_index,1)
                seaborn.heatmap(numpy_array, ax = legend_axs[ax_index], 
                                square=True, xticklabels=False, cbar=False, 
                                cmap='viridis', linewidths=1, vmin=0, vmax=len(category_list)-1)
                pyplot.yticks(
                    ticks=pyplot.yticks()[0],
                    labels=category_list[first_label_index:last_label_index],
                    fontsize=50,
                    rotation=0)
                legend_axs[ax_index].set_anchor('NW')

            if not(os.path.exists(output_folder)):
                os.makedirs(output_folder)
            pyplot.gcf()
            #fig.tight_layout()
            fig.suptitle(
                'Relative Abundance & ' + element_name + ' Richness\n', fontsize=80)
            pyplot.savefig(output_folder + element_name + stacked_ext)
            pyplot.close()

        if amr_analysis:
            self.amr_category_list.sort()
            superplot_per_analysis(self.abundance_dict_amr, self.amr_category_list, 'ARG')
        if mge_analysis:
            self.mge_category_list.sort()
            superplot_per_analysis(self.abundance_dict_mge, self.mge_category_list, 'MGE')
