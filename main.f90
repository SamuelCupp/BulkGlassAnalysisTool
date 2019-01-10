program optimization
use mpi
use structures
use interfaces
use initMopac
use randomSeed
use vandcheck
use tally
implicit none
integer,dimension(:,:),allocatable  ::  bins
integer,dimension(:),allocatable  ::  temp2
type(fragment),dimension(:),allocatable  ::  fragArray,oldPos,newPos
type(atom),dimension(:),allocatable  ::  vanderwaals
character(len=2),dimension(:),allocatable  ::  elemArray, elems
character(len=100)  ::  line
character(len=30)  ::  filename,tempname
character(len=16)  ::  endchar
character(len=10)  ::  iteration
double precision  ::  Elast, Eavg, E
double precision  ::  lowbin, binSize
double precision  ::  t1, t2
integer  ::  n, i, j, m, fragNumber, elemNumber, avgloop, avgtotal
integer  ::  opt,failrate,iflag
integer  ::  atoms, rejection, binNumber
integer  ::  ierr, base, shift, mynum, numproc

call MPI_INIT(ierr)
if (ierr.ne.0) then
        write(*,*) 'init err'
        stop
endif
call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierr)
if (ierr.ne.0) then
        write(*,*) 'size err'
        stop
endif
call MPI_COMM_RANK(MPI_COMM_WORLD,mynum,ierr)
if (ierr.ne.0) then
        write(*,*) 'rank err'
        stop
endif

! iteration variables:  i,j,m

write(endchar,*) mynum
open(unit=12,file=trim(adjustL(endchar))//'-notes.dat')

call cpu_time ( t1 )
write(12,*) 'processor ',mynum+1,' of ',numproc,' online!'

avgtotal = 1            ! Number of iterations/thread

lowbin = 0.5d0
binNumber = 60
binSize = 0.1d0

endchar = '==End File=='

elemNumber = 0
failrate=0  !fail count for mopac convergence.
rejection=0 !rejection count for the initial conditions.

call init_random_seed (mynum)

filename = 'init-'
tempname = 'temp-'
base=6
shift=base
if (mynum.ge.1000) then
   shift=base+3
   write (filename(base:shift),'(i4)') mynum
   write (tempname(base:shift),'(i4)') mynum
else if (mynum.ge.100) then
   shift=base+2
   write (filename(base:shift),'(i3)') mynum
   write (tempname(base:shift),'(i3)') mynum
else if (mynum.ge.10) then
   shift=base+1
   write (filename(base:shift),'(i2)') mynum
   write (tempname(base:shift),'(i2)') mynum
else
   write (filename(base:shift),'(i1)') mynum
   write (tempname(base:shift),'(i1)') mynum
end if

open(unit=15,file=trim(filename)//'.err')
iflag = 0

open(unit=1,file='fragments.txt',action='read')

read(1,*)
read(1,*)
read(1,*)
read(1,*)
read(1,*)
read(1,*) j

allocate(fragArray(j))

do i=1,size(fragArray)
   read(1,*)
   read(1,*) fragArray(i)%fragName
   read(1,*) fragArray(i)%natoms
   allocate(fragArray(i)%atomSet(fragArray(i)%natoms))
   do j=1,fragArray(i)%natoms
      read(1,*) fragArray(i)%atomSet(j)%atomType,fragArray(i)%atomSet(j)%location
      fragArray(i)%offset = fragArray(i)%offset + fragArray(i)%atomSet(j)%location
   enddo
   fragArray(i)%offset = fragArray(i)%offset/dble(fragArray(i)%natoms)
   do j=1,fragArray(i)%natoms
      fragArray(i)%atomSet(j)%location=fragArray(i)%atomSet(j)%location-fragArray(i)%offset
   enddo
enddo

close(1)

open(unit=1,file='input.txt',action='read')

read(1,*)
read(1,*)
read(1,*)
read(1,*)
read(1,*) opt

if ((opt.ne.0) .and. (opt.ne.1)) then
   write(12,*) "Optimization was not properly formatted. Continuing with NOOPT flag."
   opt = 0
endif

read(1,*)
read(1,*)
read(1,*) fragNumber

read(1,*)
read(1,*)
read(1,*) elemNumber

allocate(elems(elemNumber),bins(binNumber,(elemNumber*(elemNumber+1))/2))

close(1)
allocate(vanderwaals(20))

open(unit=2,file='vanderwaals.txt',action='read')

do i=1,size(vanderwaals)
   read(2,'(A)') line
   if (line(1:16)==endchar) exit
   read(line,*) vanderwaals(i)%atomType,vanderwaals(i)%location(1)
enddo
close(2)

bins = 0
Eavg = 0.0d0
E = 0.0d0

allocate(oldPos(fragNumber),newPos(fragNumber))

write(iteration,'(i3)') mynum
iteration = adjustl(iteration)
!open(unit=100,file='energy-'//trim(iteration)//'.dat')

!AVERAGING LOOP STARTS HERE
do avgloop=1,avgtotal

  call init_pos(filename,opt,newPos,fragArray,vanderwaals,fragNumber,iflag,rejection,i)
  
  if (iflag==1) then
    write(12,*) 'No match was found for atoms in fragment ',i,'.'
    write(12,*) 'Check that all atoms used are contained in the vanderwaals file.'
    call exit
  elseif (iflag==2) then
    write(12,*) 'No match was found for atoms in fragment ',i,'.'
    write(12,*) 'Check that all atoms used are contained in the file.'
    call exit
  endif
  
  write(12,*) 'failures before acceptable network:',i
  
  elems(1) = newPos(1)%atomSet(1)%atomType
  m=2
  do i=1,fragNumber
    do j=1,newPos(i)%natoms
      if(m-1>elemNumber) then
        write(*,*) m,elemNumber
        write(12,*) 'More elements were found in the data than was provided in the&
                     input file. Please recheck input and fragment files.'
        call exit
      endif
      if (all(newPos(i)%atomSet(j)%atomType .ne. elems)) then
        elems(m) = newPos(i)%atomSet(j)%atomType
        m=m+1
      endif
    enddo
  enddo 
  
  call system('/opt/mopac/MOPAC2016.exe '//trim(filename)//'.mop')
  
  !Finds energy (and positions if optimization is used)
  call readmopac(filename,newPos,fragNumber,E,opt,iflag)
  if (iflag==1) then
    write(12,*) "Warning! No energy was found in processor",mynum,"'s mopac .out file."
    call exit
  elseif(iflag==2) then
    failrate = failrate+1
    avgtotal = avgtotal-1
    cycle
  endif
  !Former command for iterative use of MOPAC.
  !Obsolete in this version of the code, but legacy still
  !present in file systems and updated for new 
  !call writemopac(tempname,newPos,fragNumber,opt)
  
  Eavg = Eavg + E
  
  call bin(newPos,fragNumber,bins,elems,lowbin,binSize,binNumber,iflag)
  
  if (iflag==1) then
    write(12,*) 'Atoms have been lost on processor ',mynum,'! An unknown error'&
                //'has occurred in the bin subroutine, found in tally.f90.'
    call exit
  endif
  
  !WRITE(100,*) E
  
enddo !STOP AVERAGING LOOP

!close(100)

Eavg = Eavg/avgtotal

CALL MPI_ALLREDUCE(Eavg,Elast,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
        MPI_COMM_WORLD,ierr)

if(ierr.ne.0) then
   write(12,*) 'reduction err for energy'
endif

allocate(temp2(1))

CALL MPI_ALLREDUCE(rejection,temp2(1),1,MPI_INTEGER,MPI_SUM,&
        MPI_COMM_WORLD,ierr)

rejection = temp2(1)

CALL MPI_ALLREDUCE(failrate,temp2(1),1,MPI_INTEGER,MPI_SUM,&
        MPI_COMM_WORLD,ierr)

failrate = temp2(1)

deallocate(temp2)

Eavg=Elast

allocate(temp2(elemNumber*(elemNumber+1)/2))
do i=1,binNumber
   call MPI_ALLREDUCE(bins(i,:),temp2,((elemNumber*(elemNumber+1))/2),MPI_INTEGER,MPI_SUM,&
           MPI_COMM_WORLD,ierr)
   if(ierr.ne.0) then
      write(12,*) 'reduction err for bin number ',i
   endif
   bins(i,:) = temp2
enddo
deallocate(temp2)

if(mynum==0) then

!Eavg=Eavg/numproc
!print*,Eavg

m=0
do i=1,elemNumber
   do j=i,elemNumber
      open(unit=50,file=trim(elems(i))//trim(elems(j))//'.dat')
write(50,*) sum(bins(:,m+j))
      do n=1,binNumber
         write(50,*) (lowbin+binSize*(n-1)),bins(n,m+j)
      enddo
      close(50)
   enddo
   m = m + elemNumber - i
enddo

endif
close(12)

do i=1,fragNumber
   deallocate(newPos(i)%atomSet)
enddo
deallocate(newPos,bins)

CALL MPI_FINALIZE(ierr)

end
