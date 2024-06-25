from tels_analysis import get_sample_name_definition
import json

class indiv_stats:
    def __init__ (this, fileName):
        this.sample_name_def = get_sample_name_definition(fileName)
        if this.sample_name_def[1] == "PacBio":
            if fileName[-2:] == 'AM': 
                this.replicate = 'B'
            elif fileName[-1] == 'A':
                this.replicate = 'A'
            else:
                this.replicate = 'C'
            
            this.sample_name_def = (
                this.sample_name_def[0], this.sample_name_def[1],
                this.sample_name_def[2], "__"
            )
        else:
            this.replicate = fileName[-1]
        this.raw_reads = 0
        this.deduplicated_reads = 0
        this.total_reads_length = 0
        this.duplication = 0
        this.ARG_on_target = 0
        this.MGE_on_target = 0

    def getStats(this):
        toReturn = ','.join(this.sample_name_def)
        toReturn = toReturn + "," + this.replicate
        toReturn = toReturn + "," + str(this.total_reads_length / 1000000.0)
        toReturn = toReturn + "," + str(this.raw_reads)
        toReturn = toReturn + "," + str(this.deduplicated_reads)
        toReturn = toReturn + "," + str(this.duplication)
        toReturn = toReturn + "," + str(this.ARG_on_target)
        toReturn = toReturn + "," + str(this.MGE_on_target)
        return toReturn

    def findAllStats(
            this, filePath_stats, filepath_reads_length, filePath_ARG, filePath_MGE):
        this.findReadStats(filePath_stats)
        this.findARGStats(filePath_ARG)
        this.findMGEStats(filePath_MGE)
        this.find_total_reads_length(filepath_reads_length)

    def find_total_reads_length(self, filepath):
        json_file = open(filepath)
        json_obj = json.load(json_file)
        self.total_reads_length = sum(list(json_obj.values()))

    def findReadStats(this, filePath):
        statFile = open(filePath, "r")
        this.raw_reads = int(statFile.readline().split(',')[1])
        this.deduplicated_reads = int(statFile.readline().split(',')[1])
        statFile.close()
        if this.raw_reads == 0:
            this.raw_reads = "__"
            this.duplication = "__"
        else:
            this.duplication = round(100 - ((this.deduplicated_reads/this.raw_reads) * 100), 1)
    def findARGStats(this, filePath):
        ARGstatFile = open(filePath, "r")
        ARGstatFile.readline()
        this.ARG_on_target = round((int(ARGstatFile.readline().split(',')[1]) / this.deduplicated_reads) * 100,1)
        ARGstatFile.close()

    def findMGEStats(this, filePath):
        MGEstatFile = open(filePath, "r")
        MGEstatFile.readline()
        this.MGE_on_target = round((int(MGEstatFile.readline().split(',')[1]) / this.deduplicated_reads) * 100, 1)
        MGEstatFile.close()