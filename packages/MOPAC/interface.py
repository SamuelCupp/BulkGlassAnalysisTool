import subprocess
from error import PackageError

def commOut(mynum,threads,opt,network):
    with open('temp-'+str(mynum)+'.mop', 'w') as out:
        out.write ("CHARGE=0 PM7 CYCLES=10000 threads="+str(threads) + (" OPT" if opt == 1 else " NOOPT") + "\n\n")
        for frag in network:
            for atom in frag.atomList:
                loc = "  "+ str(atom.loc[0]+frag.offset[0]) +" "+str(opt)+" "+ str(atom.loc[1]+frag.offset[1]) +" "+str(opt)+" "+ str(atom.loc[2]+frag.offset[2]) +" "+str(opt)
                out.write("  " + atom.elem + loc + '\n')
    return

def commRun(callpath,mynum,threads):
    subprocess.run([callpath,'temp-'+str(mynum)+'.mop'])

def commIn(mynum,opt,network,energy):
    import re
    fragnum = 0
    atomnum = 0
    with open('temp-'+str(mynum)+'.out', 'r') as f:
        for line in f:
            if re.match('\s+(TOTAL ENERGY)\s+=\s+((?:-|)[0-9]+\.[0-9]+)(.*)',line):
                m = re.match('\s+(TOTAL ENERGY)\s+=\s+((?:-|)[0-9]+\.[0-9]+)(.*)',line)
                energy += float(m.group(2))
                if(opt == 1):
                    for line in f:
                        if re.match('\s+CARTESIAN COORDINATES\s*',line):
                            next(f)
                            for line in f:
                                if re.match("\s*[0-9]+\s+([A-Z][a-z]?)\s+((?:-|)[0-9]+\.[0-9]+)\s+((?:-|)[0-9]+\.[0-9]+)\s+((?:-|)[0-9]+\.[0-9]+)",line):
                                    m = re.match("\s*[0-9]+\s+([A-Z][a-z]?)\s+((?:-|)[0-9]+\.[0-9]+)\s+((?:-|)[0-9]+\.[0-9]+)\s+((?:-|)[0-9]+\.[0-9]+)",line)
                                    if(atomnum < len(network[fragnum].atomList)):
                                        frag = network[fragnum]
                                        if(frag.atomList[atomnum].elem != m.group(1)):
                                            msg = 'Something has gone wrong with the MOPAC interface! The atoms are out of order!\n'
                                            msg = msg + "I don't know which atoms are which, so I can't store them in fragments!"
                                            PackageError(msg)
                                        frag.atomList[atomnum].loc = [float(m.group(2))-frag.offset[0],float(m.group(3))-frag.offset[1],float(m.group(4))-frag.offset[2]]
                                        atomnum += 1
                                    else:
                                        fragnum += 1
                                        atomnum = 0
                                        if(fragnum >= len(network)):
                                            msg = 'Something has gone wrong with the MOPAC interface! There seem to be more\n'
                                            msg = msg + "atoms than I started with, as I have run out of fragments in the network!"
                                            PackageError(msg)
                                elif re.match('\s*\n',line):
                                    return
                                else:
                                    msg = '\nThe following line from temp-'+str(mynum)+'.out has unexpected formatting.\n'
                                    msg = line+'\n'
                                    msg = 'Expected lines are the CARTESIAN COORDINATES for each atom followed by a blank line.'
                                    raise PackageError(msg)
