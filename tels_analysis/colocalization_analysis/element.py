class element:
    def __init__(this, element_name, element_start, element_end, element_type):
        this.name = element_name
        this.start = element_start
        this.end = element_end
        this.type = element_type
        this.color = -1
        if this.type == "ARG":
            if (this.name.find("Metals") != -1) or (this.name.find("Biocides") != -1):
                this.color = 3
            else:
                this.color = 2
        elif this.type == "MGE":
            this.color = 1
        elif this.type == "KEGG":
            this.color = 0
        
    def getStart(this):
        return this.start

    def getElementInfo(this):
        return (this.name, this.start, this.end, this.color)