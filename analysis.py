from structs import *
from math import sqrt

def pairDistance(network,bins,low,high,binSize):
    for frag in network:
        for atom in frag.atomList:
            for frag2 in network:
                for atom2 in frag2.atomList:
                    if(frag == frag2 and atom == atom2):
                        continue
                    dist = [abs(atom.loc[i]-atom2.loc[i]) for i in range(0,2)]
                    length = 0
                    for elem in dist:
                       length += elem*elem
                    length = sqrt(length)
                    if (atom.elem,atom2.elem) in bins:
                        step = 0
                        while(low + step*binSize < high):
                            if(length > low + step*binSize and length < low + (step+1)*binSize):
                                bins[(atom.elem,atom2.elem)][step] += 1
                                break
                            step += 1
                    else:
                        bins[(atom.elem,atom2.elem)] = [0 for i in range(0,int((high-low)/binSize))]
                        step = 0
                        while(low + step*binSize <= high):
                            if(length > low + step*binSize and length < low + (step+1)*binSize):
                                bins[(atom.elem,atom2.elem)][step] += 1
                                unmatched = 0
                                break
                            step += 1
