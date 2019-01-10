module structures

type  ::  atom
   character(len=2)  :: atomType
   double precision,dimension(3)  :: location
   double precision  ::  vandRad
end type atom

type :: fragment
   character(len=20)  ::  fragName
   integer  ::  natoms
   type(atom),dimension(:),allocatable  ::  atomSet
   double precision,dimension(3)  ::  offset = 0.0d0
end type fragment

interface assignment(=)
   module procedure atom_equality
   module procedure fragment_equality
end interface

contains
subroutine fragment_equality(frag1,frag2)
type(fragment),intent(out)  ::  frag1
type(fragment),intent(in)  ::  frag2

frag1%natoms = frag2%natoms
frag1%fragName = frag2%fragName
if(allocated(frag1%atomSet)) deallocate(frag1%atomSet)
allocate(frag1%atomSet(frag1%natoms))
frag1%atomSet = frag2%atomSet
end subroutine

subroutine atom_equality(atom1,atom2)
type(atom),intent(out)  ::  atom1
type(atom),intent(in)  ::  atom2

atom1%atomType = atom2%atomType
atom1%location = atom2%location
atom1%vandRad = atom2%vandRad
end subroutine

end module
