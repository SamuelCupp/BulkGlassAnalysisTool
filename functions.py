import math as m
import numpy as np
import random as r

def randomRotate(fragment):
    # The uniformly random rotation matrix used here is from 
    # Arvo, James (1991). “Fast Random Rotation Matrices”.
    # Graphics Gems III (David Kirk, ed.), pg. 117–120.
    # Academic Press, Boston.
    # Using the variables defined below, the rotation matrix is given by
    # | -a*b  |  a*d  | sqm*c |
    # |  b*c  | -c*d  | sqm*a |
    # | sqm*d | sqm*b |  -x3  |

    x1 = r.random()
    x2 = r.random()
    x3 = r.random()

    sq = m.sqrt(x3)
    sqm = m.sqrt(1.0-x3)
    temp1 = 2.0 * m.pi * x2*sq
    temp2 = 2.0 * m.pi * (x1 + x2*sq)

    a = m.sin(temp1)
    b = m.sin(temp2)
    c = m.cos(temp1)
    d = m.cos(temp2)

    Rot = [[-a*b,a*d,sqm*c],[b*c,-c*d,sqm*a],[sqm*d,sqm*b,-x3]]

    for atom in fragment.atomList:
        atom.loc = np.matmul(Rot,atom.loc)

def randomTranslate(fragment,maxsep):
    # Random angles for translation generated using first method from
    # http://corysimon.github.io/articles/uniformdistn-on-sphere/
    # "Generating uniformly distributed numbers on a sphere"
    # Accessed Sept-01-2018
    minsep = 1.0
    radius = minsep + (m.sqrt(maxsep)-minsep)*r.random()
    theta = 2 * m.pi *r.random()
    phi = m.acos(1 - 2 * r.random())
    fragment.offset = [radius*m.sin(phi) * m.cos(theta), radius*m.sin(phi)*m.sin(theta), radius*m.cos(phi)]

def posCheck(network):
    temp = 0
    for i in network[:-1]:
        temp += len(i.atomList)
    for atom in network[-1].atomList:
        p1 = [atom.loc[i] + network[-1].offset[i] for i in range(0,3)]
        test = 0
        for frag in network[:-1]:
            for atom2 in frag.atomList:
                test += 1
                p2 = [atom2.loc[i] + frag.offset[i] for i in range(0,3)]
                diff = [p1[i]-p2[i] for i in range(0,3)]
                if(abs(np.linalg.norm(diff)) < 1.0):
                    return 1
    return 0
