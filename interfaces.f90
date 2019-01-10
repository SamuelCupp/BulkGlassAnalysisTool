module interfaces
contains

subroutine readmopac(filename,newPos,fragNumber,E,opt,iflag)
use structures
integer  ::  n, i, j, fragNumber
integer  ::  opt,iflag
double precision  :: E 
character(len=30)  ::  filename
character(len=100)  :: line
character(len=65)  :: testchar,evstr,geom1,geom2,geom3,endstr
character(len=30)  ::  energychar
character(len=41)  ::  excess
integer  ::  finder, start, finish, atoms
logical :: ok,engfound
type(fragment),dimension(:),allocatable  ::  newPos

iflag = 0

testchar = 'TOTAL ENERGY            ='
endstr= ' == MOPAC DONE =='
geom1='ATOM   CHEMICAL          X               Y               Z'
geom2='NUMBER   SYMBOL      (ANGSTROMS)     (ANGSTROMS)     (ANGSTROMS)'
geom3='           Empirical Formula:'
evstr='EV'
excess='     EXCESS NUMBER OF OPTIMIZATION CYCLES'


engfound=.false.
open(unit=3,file=trim(filename)//'.out',action='read')

do
  read(3,'(A)') line
  if (line(1:17).eq.endstr) exit
  if (line(1:).eq.excess) then
    iflag = 2
    return
  endif
  finder = index(line,testchar)
  if (line(11:35).eq.testchar) then
    finder=index(line,'=')
    start=finder+1
    finish=index(line,'EV')-1
    energychar = line(start:finish)
    read(line(start:finish),*) E
    engfound=.true.
  endif
  if ((line(4:61).eq.geom1).and.(engfound).and.(opt==1)) then
    read(3,'(A)') line
    if(line(3:66).eq.geom2) then
      read(3,'(A)') line
      ok=.true.
      atoms=0
      i = 1
      do while(ok)
        read(3,'(A)') line
        finder=index(line,'*')
        atoms = atoms+1
        if (finder.gt.0) then
          if(atoms>newPos(i)%natoms) then
            atoms = atoms - newPos(i)%natoms
            i = i+1
          endif
          read(line(23:35),*) newPos(i)%atomSet(atoms)%location(1)
          read(line(39:51),*) newPos(i)%atomSet(atoms)%location(2)
          read(line(55:67),*) newPos(i)%atomSet(atoms)%location(3)
        end if
        if (finder.eq.0) exit
      enddo
    endif
  endif
enddo

close(3)

if(.not. engfound) then
  iflag = 1
  return
endif

do i=1,fragNumber
  newPos(i)%offset = 0.0d0
  do j=1,newPos(i)%natoms
     newPos(i)%offset = newPos(i)%offset + newPos(i)%atomSet(j)%location
  enddo
  newPos(i)%offset = newPos(i)%offset/dble(newPos(i)%natoms)
enddo

do i=1,fragNumber
  do j=1,newPos(i)%natoms
    newPos(i)%atomSet(j)%location = newPos(i)%atomSet(j)%location&
                                   -newPos(i)%offset
  enddo
enddo

end subroutine

subroutine writemopac(filename,Pos,fragNumber,opt)
use structures
implicit none
character(len=30)  ::  filename
character(len=60)  ::  parameters
type(fragment),dimension(:)  ::  Pos
integer  ::  fragNumber, i, j, opt, error

parameters = "CHARGE=0 SINGLET PM7 CYCLES=10000"

open(unit=1,file=trim(filename)//'.mop')
if (opt==1) then
  write(1,*) trim(parameters)//" OPT"
else
  write(1,*) trim(parameters)//" NOOPT"
endif
!Title & comments
write(1,*) ""
write(1,*) ""
do i=1,fragNumber
  do j=1,Pos(i)%natoms
    write(1,'("  ", a2, " ", f15.8, "  ", I1, " ", f15.8, "  ", I1, " ", f15.8, " ", I1)',iostat=error)&
      Pos(i)%atomSet(j)%atomType,Pos(i)%atomSet(j)%location(1)+Pos(i)%offset(1),&
      opt,Pos(i)%atomSet(j)%location(2)+Pos(i)%offset(2),opt,&
      Pos(i)%atomSet(j)%location(3)+Pos(i)%offset(3),opt
    if (error.ne.0) then
      print*,"Writing the .mop file has resulted in error code ",error
      call exit
    endif
  enddo
enddo

close(1)

end subroutine

end module
