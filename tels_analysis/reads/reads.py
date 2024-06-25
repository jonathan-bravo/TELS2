class reads:
    def __init__(this, fileName, SOURCE_PREFIX, SOURCE_SUFFIX, READS_LENGTH):
        this.readsDict = {}
        Reads_length_file = open(SOURCE_PREFIX + fileName + SOURCE_SUFFIX + READS_LENGTH, "r")
        line = Reads_length_file.readline()
        readNameBool = False
        lengthBool = False
        name = ""
        length = ""
        for c in line:
            if c == '\"':
                readNameBool = not(readNameBool)
            elif (c == ',') or (c == '}'):
                lengthBool = False
                this.readsDict.update({name:int(length)})
                name = ""
                length = ""
            elif c == ':':
                lengthBool = True
            elif c == ' ':
                continue
            elif readNameBool:
                name = name + c
            elif lengthBool:
                length = length + c
    def getLength(this, name):
        return this.readsDict[name]