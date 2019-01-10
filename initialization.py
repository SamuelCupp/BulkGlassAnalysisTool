import math
import re
from os import path
#import subprocess
import random
from error import *
from functions import *
from structs import *

def inputdata(fragments):
    "This reads in initial data from the input files."

    vanderwaals = []

    with open("vanderwaals.txt", "r") as f:
        for line in f:
            #This matches the (element) and the (radius) for each element.
            if re.match("([A-Z][a-z]?)\s+([0-9]\.[0-9]+)",line):
                m = re.match("([A-Z][a-z]?)\s+([0-9]\.[0-9]+)",line)
                vanderwaals.append([m.group(1),float(m.group(2))])
    
    with open("input.txt", "r") as f:
        for line in f:
            if re.match("package = (.*)",line):
                m = re.match("package = (.*)",line)
                package = m.group(1)
                if(not path.isdir("packages/"+package)):
                    msg = "\nThe package '"+package+"' is either not supported "
                    msg = msg+"or has been misspelled.\nPlease check for valid"
                    msg = msg+"packages in the packages/ subdirectory."
                    raise InputError(msg)
            elif re.match("threads = (.*)",line):
                m = re.match("threads = (.*)",line)
                threads = int(m.group(1))
                if(threads < 1):
                    msg = "The 'threads' value in input.txt is invalid."
                    raise InputError(msg)
            elif re.match("\s*package location\s*",line):
                m = re.match("(.*)",next(f))
                callpath = m.group(1)
                # There should be a way to test the exe with subprocess,
                # but I haven't figured out how to apply it correctly.
                # mopac needs an enter to finish, but I don't see how to pass one.
#                try:
#                    subprocess.run(callpath,timeout=0.1,check=True)
#                except subprocess.CalledProcessError:
#                    print('oh no')
#                except OSError:
#                    print('oh no2')
            elif re.match("opt = (.*)",line):
                m = re.match("opt = (.*)",line)
                flag = m.group(1)
                if (flag=='yes' or flag == '1'):
                    flag = 1
                elif (flag=='no' or flag == '0'):
                    flag = 0
                else:
                    msg = "The 'opt' value in input.txt is invalid."
                    raise InputError(msg)
                opt = flag
            elif re.match("fragment_count = (.*)",line):
                m = re.match("fragment_count = (.*)",line)
                fragNumber = int(m.group(1))
                if(fragNumber <= 0):
                    msg = "The 'fragment_count' value in input.txt is invalid."
                    raise InputError(msg)
    
    with open("fragments.txt", "r") as f:
        for line in f:
            if re.match("(\w+)\n",line):
                m = re.match("(\w+)\n",line)
                fragments.append(fragment(m.group(1),[0,0,0]))
                for line in f:
                    #This matches the (element) and the (x), (y), (z) position for each atom in the fragment.
                    if re.match("([A-Z][a-z]?)\s+((?:-|)[0-9]+\.[0-9]+)\s+((?:-|)[0-9]+\.[0-9]+)\s+((?:-|)[0-9]+\.[0-9]+)",line):
                        m = re.match("([A-Z][a-z]?)\s+((?:-|)[0-9]+\.[0-9]+)\s+((?:-|)[0-9]+\.[0-9]+)\s+((?:-|)[0-9]+\.[0-9]+)",line)
                        ierr = 1
                        for tup in vanderwaals:
                            if (tup[0] == m.group(1)):
                                index = vanderwaals.index(tup)
                                fragments[-1].addAtom(m.group(1),[float(m.group(2)),float(m.group(3)),float(m.group(4))],float(vanderwaals[index][1]))
                                ierr = 0
                        if(ierr == 1):
                            msg = '\nThe element '+m.group(1)+' was missing from vanderwaals.txt.\n'
                            msg = msg + 'Please ensure that all elements used in fragments.txt\n'
                            msg = msg + 'are also present in vanderwaals.txt.'
                            raise InputError(msg)
                    elif re.match("\s*\n",line):
                        break
                    else:
                        msg = '\nThe following atom in fragments.txt has an incorrect input formatting:\n\n'
                        msg = msg + "    " + line
                        raise InputError(msg)
        return (threads, opt,fragNumber,package,callpath)

def initial_network(network, fragNumber, fragments):
    import copy
    network.clear()
    network.append(copy.deepcopy(fragments[random.randrange(0,len(fragments))]))
    randomRotate(network[0])
    rejections = 0
    for i in range(1,int(fragNumber)):
        ierr = 1
        while(ierr == 1):
            network.append(copy.deepcopy(fragments[random.randrange(0,len(fragments))]))
            #This randomly translates fragments [1,fragNumber] angstroms. The size of the
            #random sphere should be considered in more detail for real runs.
            randomTranslate(network[i],fragNumber)
            randomRotate(network[i])
            ierr = posCheck(network)
            if(ierr):
                rejections += 1
                network.pop()
    print("Network construction finished with", rejections,"fragment placements failed.")
    names = []
    for frag in network:
        names.append(frag.name)
    from collections import Counter
    print(Counter(names))
    return rejections
