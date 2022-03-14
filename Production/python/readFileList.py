def addList(fileList, inputFiles, fileNamePrefix=None):
    """read intput file list from a another list"""
    if fileNamePrefix is None:
        for name in inputFiles:
            if len(name) > 0 and name[0] != '#':
                fileList.append(name)
    else:
        for name in inputFiles:
            if len(name) > 0 and name[0] != '#':
                print(fileNamePrefix + name)
                fileList.append(fileNamePrefix + name)
