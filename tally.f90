module tally
use structures
contains

subroutine bin(pos,fragNumber,bins,elemTypes,lowest,binSize,binNumber,iflag)
implicit none
type(fragment),dimension(:),intent(in)  ::  pos
character(len=2),dimension(:),intent(in)  ::  elemTypes
double precision,intent(in)  ::  lowest,binSize
integer,intent(in)  ::  binNumber, fragNumber
integer,dimension(:,:)  ::  bins

double precision,dimension(:,:,:),allocatable  ::  elemPosition
integer,dimension(:),allocatable :: elemCount

double precision  ::  lowbin, highbin, r
integer,dimension(:),allocatable  ::  atomCount
integer  ::  i, j, k, m, n, iflag, counter,elemNumber

iflag = 0

elemNumber = size(elemTypes)

n = 0
do i=1,fragNumber
   n = n + pos(i)%natoms
enddo

allocate(elemCount(elemNumber),elemPosition(elemNumber,n,3))

elemCount = 0
elemPosition = 0.0d0
do k=1,elemNumber
   do i=1,fragNumber
      do j=1,pos(i)%natoms
         if(pos(i)%atomSet(j)%atomType==elemTypes(k)) then
            elemCount(k) = elemCount(k)+1
            elemPosition(k,elemCount(k),:) = pos(i)%atomSet(j)%location&
                +pos(i)%offset
         endif
      enddo
   enddo
enddo

if(sum(elemCount) .ne. n) then
   iflag = 1
   return
endif

do i=1,binNumber
   lowbin = lowest + binSize*dble(i)
   highbin = lowbin + binSize
   counter = 0
   do j=1,elemNumber
      do k=1,elemCount(j)
         do m=k+1,elemCount(j)
            r = sqrt(dot_product(elemPosition(j,k,:)-elemPosition(j,m,:),&
                                 elemPosition(j,k,:)-elemPosition(j,m,:)))
            if ((r .ge. lowbin) .and. (r .lt. highbin)) then
               bins(i,1+counter) = bins(i,1+counter) + 1
            endif
         enddo !end of self-self loop (m)
         do m=j+1,elemNumber
            do n=1,elemCount(m)
               r = sqrt(dot_product(elemPosition(j,k,:)-elemPosition(j,m,:),&
                                    elemPosition(j,k,:)-elemPosition(j,m,:)))
               if ((r .ge. lowbin) .and. (r .lt. highbin)) then
                  bins(i,m+counter+1-j) = bins(i,m+counter+1-j) + 1
               endif
               enddo !loop over atoms in secondary element loop (n)
         enddo !loop over other elements (m)
      enddo !loop over atoms in primary loop (k)
      counter = counter + elemNumber + 1 - j
   enddo !outer loop over elements (j)
enddo !loop over bins (i)
deallocate(elemPosition,elemCount)
end subroutine
end module
