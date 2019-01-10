module vandcheck
use structures
contains

subroutine rejectPos(newPos,fragNumber,iflag)
implicit none
type(fragment),dimension(:),intent(in)  ::  newPos
double precision,dimension(:,:),allocatable  ::  location
double precision,dimension(:),allocatable  ::  rads
integer,intent(in)  ::  fragNumber
integer  ::  i, j, iflag, atomNumber

iflag = 0
atomNumber = 0
do i=1,fragNumber
   atomNumber = atomNumber + newPos(i)%natoms
enddo

!allocate(location(atomNumber,4),rads(atomNumber))
allocate(location(atomNumber,3),rads(atomNumber))

atomNumber = 0

do i=1,fragNumber
   do j=1,newPos(i)%natoms
!      location(j+atomNumber,1:3) = newPos(i)%atomSet(j)%location&
      location(j+atomNumber,:) = newPos(i)%atomSet(j)%location&
                            + newPos(i)%offset
!      location(j+atomNumber,4) = newPos(i)%atomSet(j)%vandRad
   enddo
   atomNumber = atomNumber + newPos(i)%natoms
enddo


do i=1,atomNumber-1
   do j=i+1,atomNumber
      rads(i) = sqrt(dot_product(location(i,1:3)-location(j,1:3),location(i,1:3)-location(j,1:3)))
!      if (rads(i)<(location(i,4) + location(j,4))/2.0d0) then
      if (rads(i)<1.0d0) then
         iflag = 1
         return
      endif
   enddo
enddo

deallocate(location,rads)

end subroutine


subroutine distcheck(newPos,fragNumber,iflag)
implicit none
type(fragment),dimension(:),intent(in)  ::  newPos
double precision,dimension(:,:),allocatable  ::  location
double precision,dimension(:),allocatable  ::  rads
integer,intent(in)  ::  fragNumber
integer  ::  i, j, iflag, atomNumber

iflag = 1
atomNumber = 0
do i=1,fragNumber
   atomNumber = atomNumber + newPos(i)%natoms
enddo

allocate(location(atomNumber,4),rads(atomNumber))

atomNumber = 0

do i=1,fragNumber
   do j=1,newPos(i)%natoms
      location(j+atomNumber,1:3) = newPos(i)%atomSet(j)%location&
                            + newPos(i)%offset
      location(j+atomNumber,4) = newPos(i)%atomSet(j)%vandRad
   enddo
   atomNumber = atomNumber + newPos(i)%natoms
enddo


do i=1,atomNumber-1
   do j=i+1,atomNumber
      rads(i) = sqrt(dot_product(location(i,1:3)-location(j,1:3),location(i,1:3)-location(j,1:3)))
      if (rads(i)<(location(i,4) + location(j,4))/2.0d0) then
         iflag = 0
         return
      endif
   enddo
enddo

deallocate(location,rads)

end subroutine
end module
