from tels_analysis.statistical_analysis.indiv_stats import indiv_stats
import os

class statistical_analyzer:
    filePath = lambda this, fileName, extension : this.source_prefix + fileName + this.source_suffix + extension

    def __init__(
            this, SOURCE_PREFIX, SOURCE_SUFFIX, SHORT_AMR_DIV, 
            SHORT_MGE, STATS, READS_LENGTH):
        this.source_prefix = SOURCE_PREFIX
        this.source_suffix = SOURCE_SUFFIX
        this.reads_length = READS_LENGTH
        this.amr_reads = SHORT_AMR_DIV
        this.mge_reads = SHORT_MGE
        this.stats = STATS
        this.statList = []

    def analyzeFile(this, fileName):
        fileStats = indiv_stats(fileName)
        if os.path.exists(this.filePath(fileName, this.stats)):
            fileStats.findAllStats(
                this.filePath(fileName, this.stats), this.filePath(fileName, this.reads_length), 
                this.filePath(fileName, this.amr_reads), this.filePath(fileName, this.mge_reads))
            this.statList.append(fileStats.getStats())

    def printAnalysis(this, STATISTICAL_ANALYSIS):
        analysis = open(STATISTICAL_ANALYSIS, "w")
        analysis.write("Organism,Sequencing platform,Chemistry,Probe,Replicate,Yield,Raw reads,De-duplicated reads,Duplication (%),ARG On-target (%),MGE On-target (%)\n")
        for stat in this.statList:
            analysis.write(stat + "\n")
        analysis.close()