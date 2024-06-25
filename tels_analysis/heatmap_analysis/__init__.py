def megares_analyzer(megaresFile):
    def sortList(mechanismsList):
        if len(mechanismsList) == 1:
            return mechanismsList
        else:
            listA = sortList(mechanismsList[0:int(len(mechanismsList)/2)])
            listB = sortList(mechanismsList[int(len(mechanismsList)/2):])
            toReturn = []
            for i in range(0,len(mechanismsList)):
                if len(listA) == 0:
                    toReturn.append(listB.pop(0))
                elif len(listB) == 0:
                    toReturn.append(listA.pop(0))
                elif listA[0][0] < listB[0][0]:
                    toReturn.append(listA.pop(0))
                else:
                    toReturn.append(listB.pop(0))
            return toReturn
    drugList = []
    otherList = []
    megares = open(megaresFile, "r")
    megares.readline()
    for line in megares:
        splitLine = line.split('|')
        tempTuple = (splitLine[2], splitLine[3])
        if splitLine[2] == "betalactams": tempTuple = ("Betalactams", splitLine[3])
        if splitLine[1] == "Drugs":
            if drugList.count(tempTuple) == 0:
                drugList.append(tempTuple)
        else:
            if otherList.count(tempTuple) == 0:
                otherList.append(tempTuple)
    megares.close()
    drugList = sortList(drugList)
    otherList = sortList(otherList)
    return (drugList, otherList)
