#testfunc

def test1():
    print('test function')
    userinput = input("Key anything for test")
    return userinput

def platformtest():
    import os
    a = os.name
    if a == 'nt':
        #print ('Windows')
        returnvalue= 'Win'
    elif a == 'posix':
        #print ('Mac')
        returnvalue= 'Mac'
    else:
        print("Unknown")
    #print(a)
    print("You are using", returnvalue)
