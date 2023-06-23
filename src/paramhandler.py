import pickle

def test():
    param_dict = {"config1": "glyco", "config2": 12345678}
    with open("config.pickle", "wb") as f:
        pickle.dump(param_dict, f)

def test2():
    param_dict = {"config1": "glyco", "config2": 12345678}
    filesaver = open("config.pickle", "wb")
    pickle.dump(param_dict, filesaver)
    filesaver.close()


def loadtest():
    loadfile = open("config.pickle", "rb")
    unpack = pickle.load(loadfile)
    loadfile.close()
    print(unpack)


def addparam(a,b):
    print("add param a and b when calling it")


def deleteconfig(filename):
    print("List the config files and let user select the file to be deleted")
