from tels_analysis.venn_analysis.indiv_venn import indiv_venn
from tels_analysis import getSampleAndIndex

class venn_analyzer:

    filePath = lambda this, fileName, extension : this.source_prefix + fileName + this.source_suffix + extension

    def __init__(this, SOURCE_PREFIX, SOURCE_SUFFIX, SHORT_AMR_DIV, SHORT_MGE):
        this.source_prefix = SOURCE_PREFIX
        this.source_suffix = SOURCE_SUFFIX
        this.amr_reads = SHORT_AMR_DIV
        this.mge_reads = SHORT_MGE
        this.venn_dict = {}

    def addToCount(this, fileName):
        sample, index = getSampleAndIndex(fileName)
        if sample not in this.venn_dict:
            this.venn_dict[sample] = indiv_venn(sample)
        this.venn_dict[sample].addToCount(this.filePath(fileName, this.amr_reads), index)

    def makeVenn(this, outputFolder, VENN):
        for sample in this.venn_dict:
            this.venn_dict[sample].findFinalCount()
            this.venn_dict[sample].makeFigure(outputFolder + "/ARG_venn_diagram", VENN)