from tels_analysis.stacked_abundance_analyzer import StackedAbundanceAnalyzer
from tels_analysis.compos_richness_analyzer import ComposRichnessAnalyzer
#from tels_analysis.colocalization_analyzer import colocalization_analyzer
from tels_analysis.linear_regression import megares_full_reference_length
from tels_analysis.linear_regression import mge_full_reference_length
from tels_analysis.statistical_analyzer import statistical_analyzer
from tels_analysis.abundance_analyzer import AbundanceAnalyzer
from tels_analysis.heatmap_analyzer import HeatmapAnalyzer
from tels_analysis.matrix_creator import matrix_creator
from tels_analysis.special_coloc import special_coloc
from tels_analysis.venn_analyzer import venn_analyzer
from tels_analysis.special import special
import multiprocessing
import configparser
import getopt
import shutil
import gzip
import sys
import os


file_list = ["BFV2AA",		"BFV2AB",		"BFV2AC",		"BFV2MA",		"BFV2MB",		"BFV2MC",
			"BFV2AMA",		"BFV2AMB",		"BFV2AMC",		"BFV2NEGA",		"BFV2NEGAM",	"BFV2NEGM",		
			"BFXTAA",		"BFXTAB",		"BFXTAC",		"BFXTMA",		"BFXTMB",		"BFXTMC",
			"BFXTAMA",		"BFXTAMB",		"BFXTAMC",		"BFXTNEGA",		"BFXTNEGAM",	"BFXTNEGM", 
			"HFV2AA",		"HFV2AB",		"HFV2AC",		"HFV2MA",		"HFV2MB",		"HFV2MC",
			"HFV2AMA",		"HFV2AMB",		"HFV2AMC", 		"HFV2NEGA",		"HFV2NEGAM",	"HFV2NEGM",
			"HFXTAA",		"HFXTAB",		"HFXTAC",		"HFXTMA",		"HFXTMB",		"HFXTMC",
			"HFXTAMA",		"HFXTAMB",		"HFXTAMC", 		"HFXTNEGA",		"HFXTNEGAM",	"HFXTNEGM",
			"MOV2AA",		"MOV2AB",		"MOV2AC",		"MOV2MA",		"MOV2MB",		"MOV2MC",
			"MOV2AMA",		"MOV2AMB",		"MOV2AMC",		"MOV2NEGA",		"MOV2NEGAM",	"MOV2NEGM",
			"MOXTAA",		"MOXTAB",		"MOXTAC",		"MOXTMA",		"MOXTMB",		"MOXTMC",
			"MOXTAMA",		"MOXTAMB",		"MOXTAMC",		"MOXTNEGA",		"MOXTNEGAM",	"MOXTNEGM",
			"SV2AA",		"SV2AB",		"SV2AC",		"SV2MA",		"SV2MB", 		"SV2MC",
			"SV2AMA",		"SV2AMB",		"SV2AMC",		"SV2NEGA",		"SV2NEGAM",		"SV2NEGM",
			"SXTAA",		"SXTAB", 		"SXTAC",		"SXTMA",		"SXTMB",		"SXTMC",
			"SXTAMA",		"SXTAMB",		"SXTAMC",		"SXTNEGA",		"SXTNEGAM",		"SXTNEGM"]

output_folder = "output/"
config_file = "config.ini"
try:
    options, args = getopt.getopt(sys.argv[1:], "hc:o:")
except getopt.GetoptError:
    print("tels_analysis.py -c <config_file> -o <output_folder>")
    sys.exit(-1)
for opt, arg in options:
    if opt == "-h":
        print("List of arguments:\n\n\n-c: config file\n-h: help\n-o: output folder\n")
        sys.exit()
    elif opt == "-c":
        config_file = arg
    elif opt == "-o":
        output_folder = arg

config = configparser.ConfigParser()
config.read(config_file)

if not os.path.exists(output_folder):
	os.makedirs(output_folder)

# Retrieving all constants from config

# Databases
ACLAME_DB = config.get("DATABASES", "ACLAME")
ICEBERG_DB = config.get("DATABASES", "ICEBERG")
PLASMID_DB = config.get("DATABASES", "PLASMID")
MEGARES_DB = config.get("DATABASES", "MEGARES")

# Annotation
MGES_ANNOTATION = config.get("ANNOTATIONS", "MGES")
MEGARES_ANNOTATION = config.get("ANNOTATIONS", "MEGARES")

# Tels Source
SOURCE_PREFIX = config.get("SOURCE_FILE", "SOURCE_PREFIX")
SOURCE_SUFFIX = config.get("SOURCE_FILE", "SOURCE_SUFFIX")
SOURCE_SUFFIX_DUP = config.get("SOURCE_FILE", "SOURCE_SUFFIX_DUP")

# Tels Source Rarefied
RAREFY_PREFIX = config.get("SOURCE_FILE", "RAREFY_PREFIX")
RAREFY_SUFFIX_ALL = config.get("SOURCE_FILE", "RAREFY_SUFFIX_ALL")
RAREFY_SUFFIX_CHEM = config.get("SOURCE_FILE", "RAREFY_SUFFIX_CHEM")
RAREFY_SUFFIX_ALL_DUP = config.get("SOURCE_FILE", "RAREFY_SUFFIX_ALL_DUP")
RAREFY_SUFFIX_CHEM_DUP = config.get("SOURCE_FILE", "RAREFY_SUFFIX_CHEM_DUP")

# Tels Source Extensions
SHORT_MGE_EXT = config.get("SOURCE_EXTENSION", "SHORT_MGE")
READS_LENGTH_EXT = config.get("SOURCE_EXTENSION", "READS_LENGTH")
SHORT_AMR_DIV_EXT = config.get("SOURCE_EXTENSION", "SHORT_AMR_DIV")

# Output Files
FILE_SIZE_OUTPUT = output_folder + config.get("SINGLE_OUTPUT_FILE", "FILE_SIZE")
AMR_MATRIX_OUTPUT = output_folder + config.get("SINGLE_OUTPUT_FILE", "AMR_MATRIX")
MGE_MATRIX_OUTPUT = output_folder + config.get("SINGLE_OUTPUT_FILE", "MGE_MATRIX")
STATISTICAL_ANALYSIS_OUTPUT = output_folder + config.get("SINGLE_OUTPUT_FILE", "STATISTICAL_ANALYSIS")

# Output Extension
VENN_EXT = config.get("OUTPUT_EXTENSION", "VENN")
STATS_EXT = config.get("OUTPUT_EXTENSION", "STATS")
VIOLIN_EXT = config.get("OUTPUT_EXTENSION", "VIOLIN")
STACKED_EXT = config.get("OUTPUT_EXTENSION", "STACKED")
HEATMAP_EXT = config.get("OUTPUT_EXTENSION", "HEATMAP")
INDIV_COMPOS_CHART_EXT = config.get("OUTPUT_EXTENSION", "INDIV_COMPOS_CHART")
COMPOS_RICHNESS_ANALYSIS_EXT = config.get("OUTPUT_EXTENSION", "COMPOS_RICHNESS_ANALYSIS")

# Parameters
AMR_ANALYSIS = config.getboolean("PARAMETERS", "AMR_ANALYSIS")
MGE_ANALYSIS = config.getboolean("PARAMETERS", "MGE_ANALYSIS")


if config.getboolean("STEPS", "STATS"):
	def calculate_stats(file_name, prefix, suffix, suffix_dup):
		stats_file_output = prefix + file_name + suffix + STATS_EXT
		if (not(os.path.exists(stats_file_output)) 
	 			and os.path.exists(prefix + file_name + suffix)):
			# Count dedup reads
			read_list_source = prefix + file_name + suffix
			input_file = open(read_list_source, "r")
			line_num = 0
			for line in input_file: line_num += 1
			input_file.close()

			# Count dup reads
			read_list_source = prefix + file_name + suffix_dup
			dup_file = open(read_list_source, "r")
			dupline_num = 0
			for line in dup_file: dupline_num += 1
			dup_file.close()

			output_file = open(stats_file_output, "w")
			output_file.write("DUPLICATED_STATS_NUM_OF_READS," + str(int(dupline_num/4)))
			output_file.write("\n")
			output_file.write("DEDUPLICATED_STATS_NUM_OF_READS," + str(int(line_num/4)))
			output_file.close()

	process = []
	for file_name in file_list:
		process.append(
			multiprocessing.Process(
				target=calculate_stats,
				args=(
					file_name, SOURCE_PREFIX, 
					SOURCE_SUFFIX, SOURCE_SUFFIX_DUP)))
		process.append(
			multiprocessing.Process(
				target=calculate_stats,
				args=(
					file_name, RAREFY_PREFIX, 
					RAREFY_SUFFIX_ALL, RAREFY_SUFFIX_ALL_DUP)))
		process.append(
			multiprocessing.Process(
				target=calculate_stats,
				args=(
					file_name, RAREFY_PREFIX, 
					RAREFY_SUFFIX_CHEM, RAREFY_SUFFIX_CHEM_DUP)))
	for i in range(0, 288, 8):
		for t in process[i:i+8]:
			t.start()
		for t in process[i:i+8]:
			t.join()

if config.getboolean("STEPS", "STATISTICAL_ANALYSIS"):
	stats_analyzer_object = statistical_analyzer(
		SOURCE_PREFIX, SOURCE_SUFFIX, SHORT_AMR_DIV_EXT,
		SHORT_MGE_EXT, STATS_EXT, READS_LENGTH_EXT)
	
	for file_name in file_list:
		stats_analyzer_object.analyzeFile(file_name)

	stats_analyzer_object.printAnalysis(STATISTICAL_ANALYSIS_OUTPUT)

if config.getboolean("STEPS", "COMPOS_RICHNESS_ANALYSIS"):
	cr_analyzer = ComposRichnessAnalyzer(
		SOURCE_PREFIX, SOURCE_SUFFIX, SHORT_AMR_DIV_EXT,
		SHORT_MGE_EXT, MGES_ANNOTATION, AMR_ANALYSIS, 
		MGE_ANALYSIS)
	
	for file_name in file_list:
		cr_analyzer.analyze_sample(file_name)
	cr_analyzer.write_richness(output_folder, COMPOS_RICHNESS_ANALYSIS_EXT)

	cr_analyzer.make_bar_charts(
		INDIV_COMPOS_CHART_EXT,
		output_folder + "/ARG_composition_new", 
		output_folder + "/MGE_composition_new")

if config.getboolean("STEPS", "HEATMAP"):
	heatmapAnalyzer = HeatmapAnalyzer(
		SOURCE_PREFIX, SOURCE_SUFFIX, SHORT_AMR_DIV_EXT,
		SHORT_MGE_EXT, MEGARES_ANNOTATION, MGES_ANNOTATION,
		AMR_ANALYSIS, MGE_ANALYSIS)
	
	for file_name in file_list:
		heatmapAnalyzer.add_to_maps(file_name)
	heatmapAnalyzer.make_maps(output_folder, HEATMAP_EXT)

if config.getboolean("STEPS", "FILE_SIZE"):
	fileOfSizes = open(FILE_SIZE_OUTPUT, "w")
	for file_name in file_list:
		size = str(os.stat(SOURCE_PREFIX + file_name + SOURCE_SUFFIX).st_size)
		fileOfSizes.write(file_name + "," + size + "\n")
	fileOfSizes.close()

if config.getboolean("STEPS", "VIOLIN"):
	abundance_analyzer = AbundanceAnalyzer(
		SOURCE_PREFIX, SOURCE_SUFFIX, FILE_SIZE_OUTPUT, 
		SHORT_AMR_DIV_EXT, SHORT_MGE_EXT, ACLAME_DB, 
		ICEBERG_DB, PLASMID_DB, MEGARES_DB, MGES_ANNOTATION)
	
	for file_name in file_list:
		abundance_analyzer.find_absolute_abundance(
			file_name, AMR_ANALYSIS, MGE_ANALYSIS)
		
	abundance_analyzer.make_violin_plot(
		output_folder, VIOLIN_EXT, AMR_ANALYSIS, MGE_ANALYSIS)

if config.getboolean("STEPS", "STACKED"):
	stackedAnalyzer = StackedAbundanceAnalyzer(
		SOURCE_PREFIX, SOURCE_SUFFIX, SHORT_AMR_DIV_EXT, 
		SHORT_MGE_EXT, STATS_EXT, MGES_ANNOTATION)
	
	for file_name in file_list:
		stackedAnalyzer.find_absolute_abundance(
			file_name, AMR_ANALYSIS, MGE_ANALYSIS)

	stackedAnalyzer.make_stacked_barplot(
		output_folder, STACKED_EXT, AMR_ANALYSIS, MGE_ANALYSIS)

if config.getboolean("STEPS", "VENN"):
	vennAnalyzer = venn_analyzer(SOURCE_PREFIX, 
								 SOURCE_SUFFIX,
								 SHORT_AMR_DIV_EXT, 
								 SHORT_MGE_EXT)
	for file_name in file_list:
		vennAnalyzer.addToCount(file_name)
	vennAnalyzer.makeVenn(output_folder, VENN_EXT)

if config.getboolean("STEPS", "MATRIX"):
	count_matrix = matrix_creator(RAREFY_PREFIX,
			       				  RAREFY_SUFFIX_ALL,
								  SHORT_AMR_DIV_EXT,
								  SHORT_MGE_EXT)
	
	for file_name in file_list:
		count_matrix.add_column(file_name)

	count_matrix.print_matrix(AMR_MATRIX_OUTPUT, MGE_MATRIX_OUTPUT)


# megares_full_reference_length(config['SOURCE_FILE']['MEGARES_FASTA'], output_folder + "/references_length.csv")
# mge_lock = multiprocessing.Lock()
# mp = list()
# for sample_name in file_list:
# 	mp.append(multiprocessing.Process(target=mge_full_reference_length, args=(sample_name, 
# 																			  config['SOURCE_FILE']['SOURCE_PREFIX'],
# 																			  config['SOURCE_FILE']['SOURCE_SUFFIX'],
# 																			  config['SOURCE_EXTENSION']['OVERLAP_OUTPUT'],
# 																			  config['SOURCE_FILE']['ACLAME_FASTA'],
# 																			  config['SOURCE_FILE']['ICEBERG_FASTA'],
# 																			  config['SOURCE_FILE']['PLASMID_FINDER_FASTA'],
# 																			  output_folder + "/references_length.csv",
# 																			  mge_lock)))
# for process in mp:
# 	process.start()
# for process in mp:
# 	process.join()

#if config.getboolean("STEPS", "COLOCALIZATION_ANALYSIS"):
#	for file_name in file_list:
#		colocalizationAnalyzer = colocalization_analyzer(file_name, SOURCE_PREFIX, 
#													   SOURCE_SUFFIX, 
#													   config.get("SOURCE_EXTENSION", "COLOCALIZATIONS"), 
#													   config.get("SOURCE_EXTENSION", "READS_LENGTH"))
#		colocalizationAnalyzer.makeChart( file_name, output_folder + "/" + 
#									  config.get("OUTPUT_FILE", "OUTPUT_PREFIX"), 
#									  config.get("OUTPUT_FILE", "OUTPUT_SUFFIX"), 
#									  config.get("OUTPUT_EXTENSION", "COLOCALIZATION_ANALYSIS"))

if config.getboolean("STEPS", "SPECIAL"):
	mgeInfo = special(SOURCE_PREFIX,
									"AlignedToMegares.csv",
									SOURCE_SUFFIX,
									SHORT_MGE_EXT,
									config.get("SOURCE_EXTENSION", "OVERLAP_OUTPUT"),
									MGES_ANNOTATION)
	colocInfo = special_coloc(SOURCE_PREFIX,
									"AlignedToMegares.csv",
									SOURCE_SUFFIX,
									config.get("SOURCE_EXTENSION", "COLOCALIZATIONS_RICHNESS"))
	for file_name in file_list:
		mgeInfo.addToMobilomeInfo(file_name)
		colocInfo.addColocInfo(file_name)
	mgeInfo.writeMobilomeInfo(output_folder + "/mobilome_info.csv")
	mgeInfo.addComparisonInfo()
	mgeInfo.writeComparisonMobilomeInfo(output_folder + "/comparison_mobilome_info.csv")
	mgeInfo.findInOverlapOutput(output_folder + "/comparison_with_other_samples.csv")
	colocInfo.writeColocInfo(output_folder + "/coloc_info.csv")