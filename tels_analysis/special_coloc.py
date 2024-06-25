class special_coloc:
    filename = lambda this, fileName: this.source_prefix + fileName + this.source_suffix + this.extension

    def __init__(this, SOURCE_PREFIX, MGEAlignedToMegares, SOURCE_SUFFIX, COLOCALIZATIONS_RICHNESS):
        this.source_prefix = SOURCE_PREFIX
        this.mge_classified_aligned_to_arg = []
        file = open(MGEAlignedToMegares, "r")
        for line in file:
            mge1, mge2 = tuple(line.split(","))
            mge2 = mge2[:-1]
            if mge2 == "":
                break
            this.mge_classified_aligned_to_arg.append(mge2)
        file.close()
        this.source_suffix = SOURCE_SUFFIX
        this.extension = COLOCALIZATIONS_RICHNESS
        this.sample_list = []
        this.coloc_info = {}

    def writeColocInfo(this, output):
        file = open(output, "w")
        for sample in this.sample_list:
            file.write(sample + "," + str(this.coloc_info[sample][0]) + "," + str(this.coloc_info[sample][1]) + "\n")
        file.close()
        

    def addColocInfo(this, sample):
        this.coloc_info.update({sample : [0, 0]})
        this.sample_list.append(sample)
        special_mge_list = []
        coloc_mge_list = []
        colocFile = open(this.filename(sample), "r")
        i = 0
        for line in colocFile:
            i += 1
            if i < 5:
                continue
            colocInfo = line.split(",")
            if colocInfo[1] not in coloc_mge_list:
                coloc_mge_list.append(colocInfo[1])
            if colocInfo[1] in this.mge_classified_aligned_to_arg:
                if colocInfo[1] not in special_mge_list:
                    special_mge_list.append(colocInfo[1])
        colocFile.close()
        this.coloc_info[sample] = [len(special_mge_list), len(coloc_mge_list)]
