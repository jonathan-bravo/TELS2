from . import element

class colocalization:
    def __init__(this, read_name, read_length, ARG, ARG_start, ARG_end, MGE_list, MGE_start_list, MGE_end_list, KEGG_list, KEGG_start_list, KEGG_end_list):
        this.elements = []
        this.length = read_length
        this.name = read_name
        this.elements.append(element.element(ARG, ARG_start, ARG_end, "ARG"))
        for i in range(0, len(MGE_list)):
            this.elements.append(element.element(MGE_list[i], MGE_start_list[i], MGE_end_list[i], "MGE"))
        for i in range(0, len(KEGG_list)):
            if KEGG_list[i] != '':
                this.elements.append(element.element(KEGG_list[i], KEGG_start_list[i], KEGG_end_list[i], "KEGG"))
        this.orderElements()

    def orderElements(this):
        def sort(eleList):
            if len(eleList) == 1:
                return eleList
            else:
                listA = sort(eleList[0:int(len(eleList)/2)])
                listB = sort(eleList[int(len(eleList)/2):])
                toReturn = []
                for i in range(0, len(eleList)):
                    if len(listA) == 0:
                        toReturn.append(listB[0])
                        listB.pop(0)
                    elif len(listB) == 0:
                        toReturn.append(listA[0])
                        listA.pop(0)
                    elif listA[0].getStart() < listB[0].getStart():
                        toReturn.append(listA[0])
                        listA.pop(0)
                    else:
                        toReturn.append(listB[0])
                        listB.pop(0)
                return toReturn
        this.elements = sort(this.elements)

    def getColocInfo(this):
        toReturn = []
        for ele in this.elements:
            temp = (this.name, this.length)
            toAdd = ele.getElementInfo() + temp
            toReturn.append(toAdd)
        return toReturn