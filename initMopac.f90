module initMopac
contains

subroutine init_pos(filename,opt,initPos,fragArray,vanderwaals,fragNumber,iflag,rejectCounter,failureValue)
use structures
use interfaces
use vandcheck
implicit none
type(fragment),dimension(:)  ::  initPos, fragArray
type(atom),dimension(:),intent(in)  ::  vanderwaals
integer, intent(in)  ::  opt,fragNumber
character(len=30),intent(in)  ::  filename
double precision,dimension(3)  ::  random
integer  ::  rejectCounter, i, j, n, iflag, t1,t2, failureValue
integer  ::  repeater

iflag = 0
rejectCounter = 0
!do same as in loop for one, but leave offset at zero. This keeps the network
!close to the center
call random_number(random(1))
n =floor(dble(size(fragArray))*random(1)+1.0d0)
initPos(1) = fragArray(n)
call rotate(initPos(1))
call vandAssign(initPos(1),iflag,vanderwaals)
if (iflag==1) then
  failureValue = n
  return
endif

do j=2,fragNumber
  n = floor((size(fragArray))*random(1)+1.0d0)
  initPos(j) = fragArray(n)
  call vandAssign(initPos(j),iflag,vanderwaals)
  if (iflag==2) then
    failureValue = n
    return
  endif

  repeater = 1
  do while(repeater==1)
    call random_number(random(1))
    call rotate(initPos(j))
    call translate(initPos(j),dble(fragNumber))
    call rejectPos(initPos(1:j),j,iflag)
    if (iflag==1) then
      rejectCounter = rejectCounter + 1
    elseif(iflag==0) then
      repeater=0
    endif
  enddo
enddo

call writemopac(filename,initPos,fragNumber,opt)

end subroutine init_pos

subroutine vandAssign(frag,iflag,vanderwaals)
use structures
type(fragment)  ::  frag
integer  ::  iflag,j,n
type(atom),dimension(:),intent(in)  :: vanderwaals
iflag = 0
do j=1,frag%natoms
   if (iflag==0) exit
   iflag = 0
   do n=1,size(vanderwaals)
      if(vanderwaals(n)%atomType==frag%atomSet(j)%atomType) then
         frag%atomSet(j)%vandRad = vanderwaals(n)%location(1)
         iflag = 1
         exit
      endif
   enddo
enddo
end subroutine vandAssign

subroutine rotate(frag)
use structures
implicit none
type(fragment)  ::  frag,rotatedfrag
double precision, dimension(3,3)  ::  rotation
double precision  ::  x, y, z, pi
integer  ::  i

pi = 4.0d0*datan(1.0d0)
call random_number(x)
x = x*2.0d0*pi

call random_number(y)
y = y*2.0d0*pi

call random_number(z)
z = z*2.0d0*pi

!Rotation matrix generated via Rx.Ry.Rz using matrices from Wikipedia, accessed
!3/3/2015.

rotation(1,1) = dcos(y)*dcos(z)
rotation(1,2) =dcos(z)*dsin(x)*dsin(y) + dcos(x)*dsin(z)
rotation(1,3) = dsin(x)*dsin(z) - dcos(x)*dcos(z)*dsin(y)
rotation(2,1) = -dcos(y)*dsin(z)
rotation(2,2) = dcos(x)*dcos(z) - dsin(x)*dsin(y)*dsin(z)
rotation(2,3) = dcos(z)*dsin(x) + dcos(x)*dsin(y)*dsin(z)
rotation(3,1) = dsin(y)
rotation(3,2) = -dcos(y)*dsin(x)
rotation(3,3) = dcos(x)*dcos(y)

do i=1,frag%natoms
   frag%atomSet(i)%location = matmul(rotation,frag%atomSet(i)%location)
enddo
end subroutine

subroutine translate(frag,b)
use structures
implicit none
type(fragment)  ::  frag
double precision,dimension(3)  ::  rvec
double precision  ::  randNum, a, b, pi

pi = 4.0d0*datan(1.0d0)
a = 4.0d0

call random_number(randNum)
rvec(1) = randNum*(b-a)+a
call random_number(randNum)
rvec(2) = randNum*2.0d0*pi
call random_number(randNum)
rvec(3) = randNum*pi
frag%offset(1) = rvec(1)*dcos(rvec(2))*dsin(rvec(3))
frag%offset(2) = rvec(1)*dsin(rvec(2))*dsin(rvec(3))
frag%offset(3) = rvec(1)*dcos(rvec(2))

end subroutine

end module initMopac
