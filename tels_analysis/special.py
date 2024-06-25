from tels_analysis import get_mge_annot_dict
import itertools

class special:
    filename = lambda this, fileName, index : this.source_prefix + fileName + this.source_suffix + this.extension[index]

    def __init__(this, SOURCE_PREFIX, MGEAlignedToMegares, SOURCE_SUFFIX, SHORT_MGE, OVERLAP_OUTPUT, MGE_CLASSIFICATION):
        this.source_prefix = SOURCE_PREFIX
        this.mge_aligned_to_megares = []
        file = open(MGEAlignedToMegares, "r")
        for line in file:
            mge1, mge2 = tuple(line.split(","))
            mge2 = mge2[:-1]
            this.mge_aligned_to_megares.append(mge1)
            if mge2 != "":
                this.mge_aligned_to_megares.append(mge2)
        file.close()
        this.source_suffix = SOURCE_SUFFIX
        this.extension = [SHORT_MGE, OVERLAP_OUTPUT]
        this.mge_dict = get_mge_annot_dict(MGE_CLASSIFICATION)
        this.sample_list = []
        this.classified_arg = {}
        this.classified_aligned_arg = {}
        this.aligned_args_in_overlap = {}
        this.classified_args_in_overlap = {}
        this.classified_aligned_args_in_overlap = {}
        this.aligned_arg = {}
        this.richness = {}
        this.old_classified_arg = {}
        this.old_classified_aligned_arg = {}
        this.old_richness = {}
        this.new_mges = []

    def findInOverlapOutput(this, output):
        for sample in this.sample_list:
            with open(this.filename(sample, 1), "r") as file:
                for mge in file:
                    if mge[:-1] in list(this.aligned_arg.keys()):
                        if mge[:-1] not in list(this.aligned_args_in_overlap.keys()):
                            this.aligned_args_in_overlap[mge[:-1]] = list()
                        this.aligned_args_in_overlap[mge[:-1]].append(sample)
                    if mge[:-1] in list(this.classified_arg.keys()):
                        if mge[:-1] not in list(this.classified_args_in_overlap.keys()):
                            this.classified_args_in_overlap[mge[:-1]] = list()
                        this.classified_args_in_overlap[mge[:-1]].append(sample)
                    if mge[:-1] in list(this.classified_aligned_arg.keys()):
                        if mge[:-1] not in list(this.classified_aligned_args_in_overlap.keys()):
                            this.classified_aligned_args_in_overlap[mge[:-1]] = list()
                        this.classified_aligned_args_in_overlap[mge[:-1]].append(sample)
        with open(output, "w") as file:
            for sample in this.sample_list:
                file.write(',' + sample)
            for mge in this.classified_args_in_overlap:
                file.write('\n' + mge)
                for sample in this.sample_list:
                    if sample not in this.classified_args_in_overlap[mge]:
                        file.write(',')
                    else:
                        file.write(',Y')
            file.write('\n')
            for mge in this.classified_aligned_args_in_overlap:
                file.write('\n' + mge)
                for sample in this.sample_list:
                    if sample not in this.classified_aligned_args_in_overlap[mge]:
                        file.write(',')
                    else:
                        file.write(',Y')
            file.write('\n')
            for mge in this.aligned_args_in_overlap:
                file.write('\n' + mge)
                for sample in this.sample_list:
                    if sample not in this.aligned_args_in_overlap[mge]:
                        file.write(',')
                    else:
                        file.write(',Y')

    def writeComparisonMobilomeInfo(this, output):
        file = open(output, "w")
        for sample in this.sample_list:
            file.write("," + sample + ",")
        file.write("\n")
        for sample in this.sample_list:
            file.write(",new,old")
        file.write("\n")
        both = []
        for mge, info in this.classified_arg.items():
            file.write(mge)
            for sample in this.sample_list:
                if sample not in info:
                    file.write(",0")
                else:
                    file.write("," + str(info[sample]))
                if mge in this.old_classified_arg:
                    both.append(mge)
                    if sample in this.old_classified_arg[mge]:
                        file.write("," + str(this.old_classified_arg[mge][sample]))
                    else:
                        file.write(",0")
                else:
                    file.write(",0")
            file.write("\n")
        for mge, info in this.old_classified_arg.items():
            if mge in both:
                continue
            file.write(mge)
            for sample in this.sample_list:
                file.write(",0")
                if sample not in info:
                    file.write(",0")
                else:
                    file.write("," + str(info[sample]))
            file.write("\n")
        file.write("RICHNESS:")
        for sample in this.richness:
            file.write("," + str(this.richness[sample]["classified_arg"]) + "," + str(this.old_richness[sample]["classified_arg"]))
        file.write("\n")
        file.write("\n")
        both.clear()
        for mge, info in this.classified_aligned_arg.items():
            file.write(mge)
            for sample in this.sample_list:
                if sample not in info:
                    file.write(",0")
                else:
                    file.write("," + str(info[sample]))
                if mge in this.old_classified_aligned_arg:
                    both.append(mge)
                    if sample in this.old_classified_aligned_arg[mge]:
                        file.write("," + str(this.old_classified_aligned_arg[mge][sample]))
                    else:
                        file.write(",0")
                else:
                    file.write(",0")
            file.write("\n")
        for mge, info in this.old_classified_aligned_arg.items():
            if mge in both:
                continue
            file.write(mge)
            for sample in this.sample_list:
                file.write(",0")
                if sample not in info:
                    file.write(",0")
                else:
                    file.write("," + str(info[sample]))
            file.write("\n")
        file.write("RICHNESS:")
        for sample in this.richness:
            file.write("," + str(this.richness[sample]["classified_aligned_arg"]) + "," + str(this.old_richness[sample]["classified_aligned_arg"]))
        file.close()

    def addComparisonInfo(this):
        file = open("MGEs_Classification/TELS2_thershold_comparison/mobilome_info_edited.csv", "r")
        for sample in this.sample_list:
            this.old_richness.update({sample:{"classified_arg":0,"classified_aligned_arg":0}})
        sample = []
        for (line, lineNum) in zip(file, range(107)):
            if lineNum == 1:
                continue
            split_line = line.split(',')
            if lineNum == 0:
                for index in range(54):
                    sample.append(split_line[index*2+1])
                    this.old_richness.update({split_line[index*2+1]:{"classified_arg":0,"classified_aligned_arg":0}})
                continue
            if split_line[0] not in this.mge_aligned_to_megares:
                if split_line[0] not in this.old_classified_arg:
                    this.old_classified_arg.update({split_line[0]:dict()})
                for index in range(54):
                    this.old_classified_arg[split_line[0]].update({sample[index]:split_line[index*2+1]})
                    if split_line[index*2+1] != '0':
                        this.old_richness[sample[index]]["classified_arg"] += 1
            else:
                if split_line[0] not in this.old_classified_aligned_arg:
                    this.old_classified_aligned_arg.update({split_line[0]:dict()})
                    this.old_classified_arg.update({split_line[0]:dict()})
                for index in range(54):
                    this.old_classified_aligned_arg[split_line[0]].update({sample[index]:split_line[index*2+1]})
                    this.old_classified_arg[split_line[0]].update({sample[index]:split_line[index*2+1]})
                    if split_line[index*2+1] != '0':
                        this.old_richness[sample[index]]["classified_aligned_arg"] += 1
                        this.old_richness[sample[index]]["classified_arg"] += 1

    def writeMobilomeInfo(this, output):
        file = open(output, "w")
        for sample in this.sample_list:
            file.write("," + sample)
        file.write("\n")
        for mge, info in this.classified_arg.items():
            file.write(mge)
            for sample in this.sample_list:
                if sample not in info:
                    file.write(",0")
                else:
                    file.write("," + str(info[sample]))
            file.write("\n")
        file.write("RICHNESS:")
        for sample in this.richness:
            file.write("," + str(this.richness[sample]["classified_arg"]))
        file.write("\n")
        file.write("\n")
        for mge, info in this.classified_aligned_arg.items():
            file.write(mge)
            for sample in this.sample_list:
                if sample not in info:
                    file.write(",0")
                else:
                    file.write("," + str(info[sample]))
            file.write("\n")
        file.write("RICHNESS:")
        for sample in this.richness:
            file.write("," + str(this.richness[sample]["classified_aligned_arg"]))
        file.write("\n")
        file.write("\n")
        for mge, info in this.aligned_arg.items():
            file.write(mge)
            for sample in this.sample_list:
                if sample not in info:
                    file.write(",0")
                else:
                    file.write("," + str(info[sample]))
            file.write("\n")
        file.write("RICHNESS:")
        for sample in this.richness:
            file.write("," + str(this.richness[sample]["aligned_arg"]))
        file.close()
        print(this.new_mges)

    def addToMobilomeInfo(this, sample):
        this.richness.update({sample: {"classified_arg": 0, "classified_aligned_arg": 0, "aligned_arg": 0}})
        this.sample_list.append(sample)
        file = open(this.filename(sample, 0), "r")
        i = 0
        for line in file:
            i += 1
            if i < 21:
                continue
            mge = line.split(',')[0]
            count = line.split(',')[1][:-1]
            if (mge not in this.mge_dict) and (mge not in this.new_mges):
                this.new_mges.append(mge)
            annot = "UNKNOWN" if mge not in this.mge_dict else this.mge_dict[mge]
            if annot == "AMR":
                if mge not in this.mge_aligned_to_megares:
                    if mge not in this.classified_arg:
                        this.classified_arg.update({mge: dict()})
                    this.classified_arg[mge].update({sample: count})
                    this.richness[sample]["classified_arg"] += 1
                else:
                    if mge not in this.classified_aligned_arg:
                        this.classified_arg.update({mge: dict()})
                        this.classified_aligned_arg.update({mge: dict()})
                        this.aligned_arg.update({mge: dict()})
                    this.classified_arg[mge].update({sample: count})
                    this.classified_aligned_arg[mge].update({sample: count})
                    this.aligned_arg[mge].update({sample: count})
                    this.richness[sample]["classified_arg"] += 1
                    this.richness[sample]["classified_aligned_arg"] += 1
                    this.richness[sample]["aligned_arg"] += 1
            else:
                if mge in this.mge_aligned_to_megares:
                    if mge not in this.aligned_arg:
                        this.aligned_arg.update({mge: dict()})
                    this.aligned_arg[mge].update({sample: count})
                    this.richness[sample]["aligned_arg"] += 1
        file.close()
        