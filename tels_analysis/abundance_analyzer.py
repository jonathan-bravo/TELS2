from tels_analysis.abundance_analysis.indiv_abundance import IndivAbundance
from tels_analysis import get_sample_name_definition
from tels_analysis import get_mge_annot_dict
from tels_analysis import get_genes_length
from tels_analysis import tels_file_path
from matplotlib import pyplot
#pyplot.rcParams["font.family"] = "Times New Roman"
import seaborn
import pandas
import csv
import os

class AbundanceAnalyzer:
    
    def __init__(
            self, SOURCE_PREFIX, SOURCE_SUFFIX, FILE_SIZE_OUTPUT, AMR_DIV_EXT, 
            MGE_EXT, ACLAME_DB, ICEBERG_DB, PLASMID_DB, MEGARES_DB, MGES_ANNOTATION):
        
        # Retrieve size of dedup.fastq files
        self.all_file_sizes = dict()
        with open(FILE_SIZE_OUTPUT, "r") as file_of_sizes:
            csvreader = csv.reader(file_of_sizes)
            for row in csvreader:
                self.all_file_sizes.update({row[0]:row[1]})

        # Get MGE types
        self.mges_annot = get_mge_annot_dict(MGES_ANNOTATION)

        # Get tels output file names information
        self.source_prefix = SOURCE_PREFIX
        self.source_suffix = SOURCE_SUFFIX
        self.amr_reads_ext = AMR_DIV_EXT
        self.mge_reads_ext = MGE_EXT

        # Get gene lengths for mges and megares
        self.mge_genes_length = get_genes_length(ACLAME_DB)
        self.mge_genes_length.update(get_genes_length(ICEBERG_DB))
        self.mge_genes_length.update(get_genes_length(PLASMID_DB))

        self.megares_genes_length = get_genes_length(MEGARES_DB)

        # Estabilsh dictionaries for absolute abundance        
        self.abundance_dict_amr = {"BF":{}, "FMT":{}, "PPS":{}}
        self.abundance_dict_mge = {"BF":{}, "FMT":{}, "PPS":{}}
        self.initial_source_size = {"BF":{}, "FMT":{}, "PPS":{}}

    def find_absolute_abundance(self, sample_name, amr_analysis, mge_analysis):
        # Retrieve definition in tuple form Organism, Platform, Chemistry, Probe)
        # and save the organism (Bovine, Human, Soil, or Mock)
        # + probe type (PacBio, XT, or V2 and ARG, ARG-MGE, or MGE)
        sample_name_definition = get_sample_name_definition(sample_name, True)
        organism = sample_name_definition[0]
        if organism == "MOCK": return None
        if sample_name_definition[1] == 'PacBio':
            probe_type = sample_name_definition[1]
        else:
            # probe_type = (sample_name_definition[1] + ' ' 
            #             + sample_name_definition[2] + ' ' 
            #             + sample_name_definition[3])
            probe_type = (sample_name_definition[1] + ' ' 
                        + sample_name_definition[2] + ' ' 
                        + sample_name_definition[3] + ' '
                        + sample_name_definition[4])

        # If we are analyzing AMRs
        if amr_analysis:

            # This is the first time we see this probe for the organism for amr
            if probe_type not in list(self.abundance_dict_amr[organism].keys()):
                self.abundance_dict_amr[organism].update({probe_type:IndivAbundance(True)})

            self.abundance_dict_amr[organism][probe_type].add_to_absolute(
                tels_file_path(self, sample_name, self.amr_reads_ext))

        # If we are analyzing MGEs
        if mge_analysis:

            # This is the first time we see this probe for the organism for mge
            if probe_type not in list(self.abundance_dict_mge[organism].keys()):
                self.abundance_dict_mge[organism].update({probe_type:IndivAbundance(
                    False, self.mges_annot)})

            self.abundance_dict_mge[organism][probe_type].add_to_absolute(
                tels_file_path(self, sample_name, self.mge_reads_ext))
            
        # Overall, tThis is the first time we see this probe for the organism
        if probe_type not in list(self.initial_source_size[organism].keys()):
            self.initial_source_size[organism].update({probe_type:0}) 

        self.initial_source_size[organism][probe_type] += (
            float(self.all_file_sizes[sample_name]) / (10.0**9))

    def make_violin_plot(
            self, output_folder, violin_ext, amr_analysis, mge_analysis):

        def violin_plot_per_analysis(abundance_dict, genes_length, element_name):

            #print(abundance_dict)

            # Go through absolute abundace information to make it relative
            for organism in abundance_dict:
                #print(organism)
                for probe_type in abundance_dict[organism]:
                    #print(probe_type)
                    abundance_dict[organism][probe_type].make_abundance_relative(
                        self.initial_source_size[organism][probe_type], genes_length)
                    #print(f'{organism} - {probe_type}: {abundance_dict[organism][probe_type].get_abundance()}')
            

            # Set up main figure        
            fig, axs = pyplot.subplots(
                nrows=1,
                ncols=3,
                sharey="row",
                figsize=(60, 20),
                layout="constrained")
            fig.suptitle(element_name + " Relative Abundance\n", fontsize=80)
            fig.add_subplot(111, frame_on=False)
            pyplot.grid(False)
            pyplot.tick_params(labelcolor="none", bottom=False, left=False)
            for ax in axs.flat:
                ax.set_ylabel ("Log Relative Abundance", fontsize=60)
                ax.label_outer()
                ax.set_ylim(-4,6)

            # Set up individual violin plot
            legend_artists = list()

            for index, organism in enumerate(abundance_dict):
                #print(f'{index} {organism}')
                pyplot.sca(axs[index])
                x_axis_list = list()
                sample_abundance_dict = dict()
                for probe_type in abundance_dict[organism]:
                    #print(probe_type)
                    sample_abundance_dict[probe_type] = (
                        abundance_dict[organism][probe_type].get_abundance())
                    #print(sample_abundance_dict[probe_type])
                    x_axis_list.append(probe_type)
                df = pandas.DataFrame.from_dict(sample_abundance_dict)
                #print(df)
                seaborn.set_style("whitegrid")
                seaborn.set_context("paper")

                # custom_palette = {
                #     "TELSeq V2 RES": seaborn.color_palette("BuGn", 3)[0],
                #     "TELSeq V2 MOB": seaborn.color_palette("BuGn", 3)[1],
                #     "TELSeq V2 Combo": seaborn.color_palette("BuGn", 3)[2],
                #     "PacBio": 'grey',
                #     "TELSeq XT RES": seaborn.color_palette("Oranges", 3)[2],
                #     "TELSeq XT MOB": seaborn.color_palette("Oranges", 3)[1],
                #     "TELSeq XT Combo": seaborn.color_palette("Oranges", 3)[0]
                # }
                # custom_palette = {
                #     "TELSeq V2 RES A": seaborn.color_palette("BuGn", 3)[0],
                #     "TELSeq V2 MOB A": seaborn.color_palette("BuGn", 3)[1],
                #     "TELSeq V2 Combo A": seaborn.color_palette("BuGn", 3)[2],
                #     "TELSeq V2 RES B": seaborn.color_palette("BuGn", 3)[0],
                #     "TELSeq V2 MOB B": seaborn.color_palette("BuGn", 3)[1],
                #     "TELSeq V2 Combo B": seaborn.color_palette("BuGn", 3)[2],
                #     "TELSeq V2 RES C": seaborn.color_palette("BuGn", 3)[0],
                #     "TELSeq V2 MOB C": seaborn.color_palette("BuGn", 3)[1],
                #     "TELSeq V2 Combo C": seaborn.color_palette("BuGn", 3)[2],
                #     "PacBio": 'grey',
                #     "TELSeq XT RES A": seaborn.color_palette("Oranges", 3)[2],
                #     "TELSeq XT MOB A": seaborn.color_palette("Oranges", 3)[1],
                #     "TELSeq XT Combo A": seaborn.color_palette("Oranges", 3)[0],
                #     "TELSeq XT RES B": seaborn.color_palette("Oranges", 3)[2],
                #     "TELSeq XT MOB B": seaborn.color_palette("Oranges", 3)[1],
                #     "TELSeq XT Combo B": seaborn.color_palette("Oranges", 3)[0],
                #     "TELSeq XT RES C": seaborn.color_palette("Oranges", 3)[2],
                #     "TELSeq XT MOB C": seaborn.color_palette("Oranges", 3)[1],
                #     "TELSeq XT Combo C": seaborn.color_palette("Oranges", 3)[0]
                # }
                custom_palette = {
                    "TELSeq V2 RES A": '#F53224',
                    "TELSeq V2 RES B": '#F53224',
                    "TELSeq V2 RES C": '#F53224',
                    "TELSeq V2 MOB A": '#F8766D',
                    "TELSeq V2 MOB B": '#F8766D',
                    "TELSeq V2 MOB C": '#F8766D',
                    "TELSeq V2 Combo A": '#FCBBB6',
                    "TELSeq V2 Combo B": '#FCBBB6',
                    "TELSeq V2 Combo C": '#FCBBB6',
                    "PacBio": '#619CFF',
                    "TELSeq XT RES A": '#006E21',
                    "TELSeq XT RES B": '#006E21',
                    "TELSeq XT RES C": '#006E21',
                    "TELSeq XT MOB A": '#00BA38',
                    "TELSeq XT MOB B": '#00BA38',
                    "TELSeq XT MOB C": '#00BA38',
                    "TELSeq XT Combo A": '#07FF52',
                    "TELSeq XT Combo B": '#07FF52',
                    "TELSeq XT Combo C": '#07FF52'
                }
                axs[index] = seaborn.violinplot(
                        data=df, inner="box", palette=custom_palette)
                axs[index].grid()
                axs[index].set_title(organism, fontsize=70)

                if (index == 1) and (element_name == 'ARG'):
                    artists = axs[index-1].get_children()
                    #print(artists)
                    # legend_artists = [
                    #     artists[0], artists[4], artists[8],
                    #     artists[12], artists[16], artists[20],
                    #     artists[24]]
                    legend_artists = [
                        artists[0], artists[12], artists[24],
                        artists[32], artists[40], artists[52],
                        artists[64]]
                    #ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),ncol=3, fancybox=True, shadow=True)
                    # pyplot.legend(
                    #     handles=legend_artists,
                    #     labels = ['XT-HS2 RES', 'XT-HS2 MOB', 'XT-HS2 Combo', 'Non-enriched', 'XT RES', 'XT MOB', 'XT Combo'],
                    #     #labels=x_axis_list,
                    #     loc='center',
                    #     #bbox_to_anchor=(0, 0),
                    #     # fancybox=True,
                    #     # shadow=True,
                    #     frameon=False,
                    #     fontsize=45,
                    #     ncols=1)

                axs[index].set_xticklabels(labels=[])
                
                if index == 0:
                    pyplot.yticks(fontsize=50)
            

            if not(os.path.exists(output_folder)):
                os.makedirs(output_folder)
            pyplot.gcf()
            fig.legend(handles=legend_artists, labels = ['XT-HS2 RES', 'XT-HS2 MOB', 'XT-HS2 Combo', 'Non-enriched', 'XT RES', 'XT MOB', 'XT Combo'], loc='center', bbox_to_anchor=(.515,.05), frameon=False, fontsize=60, ncols=7)
            fig.tight_layout()
            fig.subplots_adjust(bottom=0.1)
            #pyplot.rcParams["font.family"] = "Times New Roman"
            pyplot.savefig(output_folder + element_name + violin_ext, transparent=True)
            pyplot.close()
            

        if amr_analysis:
            violin_plot_per_analysis(
                self.abundance_dict_amr, self.megares_genes_length, 'ARG')
        if mge_analysis:
            violin_plot_per_analysis(
                self.abundance_dict_mge, self.mge_genes_length, 'MGE')