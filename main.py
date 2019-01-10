import importlib
import re
from structs import *
from initialization import *
from error import *
from analysis import *

# Below are the mpi calls from Fortran code.
#call MPI_INIT(ierr)
#if (ierr.ne.0) then
#        write(*,*) 'init err'
#        stop
#call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierr)
#if (ierr.ne.0) then
#        write(*,*) 'size err'
#        stop
#call MPI_COMM_RANK(MPI_COMM_WORLD,mynum,ierr)
#if (ierr.ne.0) then
#        write(*,*) 'rank err'
#        stop

mynum = 0 #placeholder for processor number
#beginning at line 51 of main.f90 of master branch

avgTotal = 1 #number of iterations/thread
lowBin = 0.5
highBin = 6.0
binSize = 0.1

pairBins = {}

opt = 1        #default averaging method

rejection=0 #rejection count for the initial conditions.

Eavg = 0 #average network energy
fragments = []
network = []

threads, opt, fragNumber, package, callpath = inputdata(fragments)
commIn = importlib.import_module("packages."+package+".interface").commIn
commOut = importlib.import_module("packages."+package+".interface").commOut
commRun = importlib.import_module("packages."+package+".interface").commRun

if opt:
    for i in range(0,avgTotal):
        rejection += initial_network(network, fragNumber, fragments)
#        commOut(mynum,threads,opt,network)
#        commRun(callpath,mynum,threads)
        commIn(mynum,opt,network,Eavg)
        pairDistance(network,pairBins,lowBin,highBin,binSize)
else:
    print("The NOOPT feature is non-functional at this time. Please run with opt = 1")

Eavg = Eavg/avgTotal

avgDist = {}
for key in pairBins:
    subtotal = 0
    pairCount = 0
    for value in pairBins[key]:
        subtotal += value * (lowBin + (pairBins[key].index(value)+0.5)*binSize)
        pairCount += value
    avgDist[key] = subtotal/pairCount
print(avgDist)

# Below are the mpi calls from Fortran code.
#CALL MPI_ALLREDUCE(Eavg,Elast,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
#        MPI_COMM_WORLD,ierr)
#if(ierr.ne.0) then
#   write(12,*) 'reduction err for energy'
#endif
#
#allocate(temp2(1))
#CALL MPI_ALLREDUCE(rejection,temp2(1),1,MPI_INTEGER,MPI_SUM,&
#        MPI_COMM_WORLD,ierr)
#rejection = temp2(1)
#
#CALL MPI_ALLREDUCE(failrate,temp2(1),1,MPI_INTEGER,MPI_SUM,&
#        MPI_COMM_WORLD,ierr)
#failrate = temp2(1)
#
#allocate(temp2(elemNumber*(elemNumber+1)/2))
#do i=1,binNumber
#   call MPI_ALLREDUCE(bins(i,:),temp2,((elemNumber*(elemNumber+1))/2),MPI_INTEGER,MPI_SUM,&
#           MPI_COMM_WORLD,ierr)
#   if(ierr.ne.0) then
#      write(12,*) 'reduction err for bin number ',i
#   endif
#   bins(i,:) = temp2
