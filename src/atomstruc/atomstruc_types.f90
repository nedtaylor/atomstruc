module atomstruc__types
  !! Module to store, read and write geometry files
  !!
  !! This module contains the procedures to read and write geometry files.
  !! It also contains the derived types used to store the geometry data.
  use coreutils, only: real32, pi, strip_null, stop_program, inverse_3x3
  use atomstruc__elements, only: get_element_properties
  implicit none


  private

  public :: basis_type, species_type



  type :: species_type
     !! Derived type to store information about a species/element.
     integer, allocatable, dimension(:) :: atom_idx
     !! The indices of the atoms of this species in the basis.
     !! For ASE compatibility.
     logical, dimension(:), allocatable :: atom_mask
     !! The mask of the atoms of this species.
     real(real32), allocatable ,dimension(:,:) :: atom
     !! The atomic positions of the species (axis, atom).
     real(real32), allocatable, dimension(:,:) :: force
     !! The forces on the atoms of the species (axis, atom).
     real(real32) :: mass
     !! The mass of the species.
     real(real32) :: charge
     !! The charge of the species.
     real(real32) :: radius
     !! The radius of the species.
     character(len=3) :: name
     !! The name of the species.
     integer :: num
     !! The number of atoms of this species.
  end type species_type
  type :: basis_type
     !! Derived type to store information about a basis.
     type(species_type), allocatable, dimension(:) :: spec
     !! Information about each species in the basis.
     integer :: nspec = 0
     !! The number of species in the basis.
     integer :: natom = 0
     !! The number of atoms in the basis.
     real(real32) :: energy = 0._real32
     !! The energy of the basis.
     real(real32) :: lat(3,3) = 0._real32
     !! The lattice vectors of the basis.
     logical :: lcart = .false.
     !! Boolean whether the basis is in cartesian coordinates.
     logical, dimension(3) :: pbc = .true.
     !! Boolean whether the basis has periodic boundary conditions.
     character(len=128) :: sysname = "default"
     !! The name of the system.
   contains
     procedure, pass(this) :: allocate_species
     !! Procedure to allocate the species in the basis.
     procedure, pass(this) :: convert
     !! Procedure to convert the basis to cartesian coordinates.
     procedure, pass(this) :: copy
     !! Procedure to copy the basis.
     procedure, pass(this) :: set_atom_mask
     !! Procedure to set the atom mask of the basis.
     procedure, pass(this) :: get_lattice_constants
     !! Procedure to get the lattice constants of the basis.
     procedure, pass(this) :: add_atom
     !! Procedure to add an atom to the basis.
     procedure, pass(this) :: remove_atom
     !! Procedure to remove an atom from the basis.
     procedure, pass(this) :: remove_atoms
     !! Procedure to remove atoms from the basis.
     procedure, pass(this) :: set_element_properties_to_default
     !! Procedure to set the element properties to default values.
  end type basis_type


  interface basis_type
     module function init_basis_type(basis) result(output)
       !! Initialise the basis type.
       type(basis_type), intent(in), optional :: basis
       !! Optional. Basis to copy.
       type(basis_type) :: output
       !! The basis to initialise.
     end function init_basis_type
  end interface basis_type



contains

!###############################################################################
  module function init_basis_type(basis) result(output)
    !! Initialise the basis type.
    implicit none

    ! Arguments
    type(basis_type), intent(in), optional :: basis
    !! Optional. Basis to copy.
    type(basis_type) :: output
    !! The basis to initialise.

    if(present(basis)) call output%copy(basis)

  end function init_basis_type
!###############################################################################


!###############################################################################
  subroutine allocate_species( &
       this, num_species, &
       species_symbols, species_count, atoms, atom_idx_list )
    !! Allocate the species in the basis.
    implicit none

    ! Arguments
    class(basis_type), intent(inout) :: this
    !! Parent. The basis to allocate the species in.
    integer, intent(in), optional :: num_species
    !! Optional. The number of species in the basis.
    character(3), dimension(:), intent(in), optional :: species_symbols
    !! Optional. The symbols of the species.
    integer, dimension(:), intent(in), optional :: species_count
    !! Optional. The number of atoms of each species.
    real(real32), dimension(:,:), intent(in), optional :: atoms
    !! Optional. The atomic positions of the species.
    integer, dimension(:), intent(in), optional :: atom_idx_list
    !! Optional. The indices of the atoms of the species.

    ! Local variables
    integer :: i, j, istart, iend
    !! Loop index.

    if(present(num_species)) this%nspec = num_species

    if(allocated(this%spec)) deallocate(this%spec)
    allocate(this%spec(this%nspec))

    species_check: if(present(species_symbols))then
       if(size(species_symbols).ne.this%nspec) exit species_check
       this%spec(:)%name = species_symbols
    end if species_check

    natom_check: if(present(species_count))then
       if(size(species_count).ne.this%nspec) exit natom_check
       this%spec(:)%num = species_count
       istart = 1
       do i = 1, this%nspec
          iend = istart + this%spec(i)%num - 1
          allocate(this%spec(i)%atom_mask(this%spec(i)%num), source = .true.)
          allocate(this%spec(i)%atom_idx(this%spec(i)%num))
          allocate(this%spec(i)%atom(3,this%spec(i)%num))
          if(present(atoms))then
             this%spec(i)%atom = atoms(:3,istart:iend)
          end if
          if(present(atom_idx_list))then
             this%spec(i)%atom_idx = atom_idx_list(istart:iend)
          else
             this%spec(i)%atom_idx = [ ( j, j = istart, iend, 1 ) ]
          end if
          istart = iend + 1
       end do
    end if natom_check

    do i = 1, this%nspec
       call get_element_properties( &
            this%spec(i)%name, &
            mass = this%spec(i)%mass, &
            charge = this%spec(i)%charge, &
            radius = this%spec(i)%radius )
    end do

  end subroutine allocate_species
!###############################################################################



!###############################################################################
  subroutine convert(this)
    !! Convert the basis between direct and cartesian coordinates.
    implicit none

    ! Arguments
    class(basis_type), intent(inout) :: this
    !! Parent. The basis to convert.

    ! Local variables
    integer :: is, ia
    !! Loop index.
    real(real32), dimension(3,3) :: lattice
    !! The reciprocal lattice vectors.


    if(this%lcart)then
       lattice = inverse_3x3( this%lat )
    else
       lattice = this%lat
    end if

    this%lcart = .not.this%lcart
    do is = 1, this%nspec
       do ia = 1, this%spec(is)%num
          this%spec(is)%atom(1:3,ia) = &
               matmul( this%spec(is)%atom(1:3,ia), lattice )
       end do
    end do

  end subroutine convert
!###############################################################################


!###############################################################################
  function get_lattice_constants(this, radians) result(output)
    !! Convert the lattice from lattice matrix to abc and αβγ.
    implicit none

    ! Arguments
    class(basis_type), intent(in) :: this
    !! Parent. The basis.
    logical, intent(in), optional :: radians
    !! Optional. Boolean whether to return angles in radians.
    real(real32), dimension(2,3) :: output
    !! The lattice constants and angles.

    ! Local variables
    logical :: radians_
    !! Boolean whether to return angles in radians.


    radians_ = .true.
    if(present(radians)) radians_ = radians

    output = convert_lat_to_abc(this%lat, radians_)

  end function get_lattice_constants
!-------------------------------------------------------------------------------
  function convert_lat_to_abc(lattice, radians) result(abc_angle)
    !! Convert the lattice from lattice matrix to abc and αβγ.
    implicit none

    ! Arguments
    real(real32), dimension(3,3), intent(in) :: lattice
    !! The lattice matrix.
    logical, intent(in), optional :: radians
    !! Optional. Boolean whether to return angles in radians.
    real(real32), dimension(2,3) :: abc_angle
    !! The lattice constants and angles.

    ! Local variables
    integer :: i
    !! Loop index.


    do i = 1, 3
       abc_angle(1,i)=norm2(lattice(i,:))
    end do
    do i = 1, 3
    end do
    abc_angle(2,1)=acos(dot_product(lattice(2,:),lattice(3,:))/&
         (abc_angle(1,2)*abc_angle(1,3)))
    abc_angle(2,3)=acos(dot_product(lattice(1,:),lattice(3,:))/&
         (abc_angle(1,1)*abc_angle(1,3)))
    abc_angle(2,3)=acos(dot_product(lattice(1,:),lattice(2,:))/&
         (abc_angle(1,1)*abc_angle(1,2)))

    if(present(radians))then
       if(.not.radians) abc_angle(2,:)=abc_angle(2,:)*180._real32/pi
    end if

  end function convert_lat_to_abc
!###############################################################################


!###############################################################################
  subroutine copy(this, basis, length)
    !! Copy the basis.
    implicit none

    ! Arguments
    class(basis_type), intent(inout) :: this
    !! Parent. The basis to copy into.
    class(basis_type), intent(in) :: basis
    !! The basis to copy from.
    integer, intent(in), optional :: length
    !! The dimension of the basis atom positions.


    ! Local variables
    integer :: i, j
    !! Loop indices.
    integer :: length_, length_input
    !! The dimension of the basis atom positions.


    !---------------------------------------------------------------------------
    ! determines whether user wants output basis extra translational dimension
    !---------------------------------------------------------------------------
    length_input = size(basis%spec(1)%atom,dim=1)
    if(present(length))then
       length_ = length
    else
       length_ = length_input
    end if


    !---------------------------------------------------------------------------
    ! if already allocated, deallocates output basis
    !---------------------------------------------------------------------------
    if(allocated(this%spec))then
       do i = 1, this%nspec
          if(allocated(this%spec(i)%atom_mask)) &
               deallocate(this%spec(i)%atom_mask)
          if(allocated(this%spec(i)%atom_idx)) deallocate(this%spec(i)%atom_idx)
          if(allocated(this%spec(i)%atom)) deallocate(this%spec(i)%atom)
          if(allocated(this%spec(i)%force)) deallocate(this%spec(i)%force)
       end do
       deallocate(this%spec)
    end if


    !---------------------------------------------------------------------------
    ! allocates output basis and clones data from input basis to output basis
    !---------------------------------------------------------------------------
    allocate(this%spec(basis%nspec))
    do i = 1, basis%nspec
       allocate(this%spec(i)%atom_mask(basis%spec(i)%num), source = .true.)
       allocate(this%spec(i)%atom_idx(basis%spec(i)%num))
       allocate(this%spec(i)%atom(length_,basis%spec(i)%num))
       allocate(this%spec(i)%force(3,basis%spec(i)%num))

       if(allocated(basis%spec(i)%atom_mask)) &
            this%spec(i)%atom_mask = basis%spec(i)%atom_mask

       if(allocated(basis%spec(i)%atom_idx))then
          this%spec(i)%atom_idx = basis%spec(i)%atom_idx
       else
          this%spec(i)%atom_idx = [ ( j, j = sum(basis%spec(1:i-1:1)%num) + 1, &
               sum(basis%spec(1:i)%num) ) ]
       end if
       if(length_input.eq.length_)then
          this%spec(i)%atom(:length_,:) = basis%spec(i)%atom(:length_,:)
       elseif(length_input.gt.length_)then
          this%spec(i)%atom(:3,:) = basis%spec(i)%atom(:3,:)
       elseif(length_input.lt.length_)then
          this%spec(i)%atom(:3,:) = basis%spec(i)%atom(:3,:)
          this%spec(i)%atom(4,:) = 1._real32
       end if
       if(allocated(basis%spec(i)%force)) &
            this%spec(i)%force(:,:) = basis%spec(i)%force(:,:)
       this%spec(i)%num = basis%spec(i)%num
       this%spec(i)%name = strip_null(basis%spec(i)%name)

       this%spec(i)%mass = basis%spec(i)%mass
       this%spec(i)%charge = basis%spec(i)%charge
       this%spec(i)%radius = basis%spec(i)%radius
    end do
    this%nspec = basis%nspec
    this%natom = basis%natom
    this%lcart = basis%lcart
    this%sysname = basis%sysname
    this%energy = basis%energy
    this%lat = basis%lat
    this%pbc = basis%pbc

  end subroutine copy
!###############################################################################


!###############################################################################
  subroutine set_atom_mask(this, index_list)
    !! Set the mask for the atoms in the basis.
    implicit none

    ! Arguments
    class(basis_type), intent(inout) :: this
    !! Parent. The basis.
    integer, dimension(:,:), intent(in), optional :: index_list
    !! The list of indices to set the mask for.

    ! Local variables
    integer :: i
    !! Loop index.


    do i = 1, this%nspec
       if(.not.allocated(this%spec(i)%atom_mask))then
          allocate( &
               this%spec(i)%atom_mask(this%spec(i)%num), source = .true. &
          )
       end if
    end do

    if(present(index_list))then
       do i = 1, size(index_list,2)
          this%spec(index_list(1,i))%atom_mask(index_list(2,i)) = &
               .not.this%spec(index_list(1,i))%atom_mask(index_list(2,i))
       end do
    end if

  end subroutine set_atom_mask
!###############################################################################


!###############################################################################
  subroutine add_atom(this, species, position, is_cartesian, mask)
    !! Add an atom to the basis.
    implicit none

    ! Arguments
    class(basis_type), intent(inout) :: this
    !! Parent. The basis.
    character(len=3), intent(in) :: species
    !! The species of the atom to add.
    real(real32), dimension(3), intent(in) :: position
    !! The position of the atom to add.
    logical, intent(in), optional :: is_cartesian
    !! Optional. Boolean whether the position is in cartesian coordinates.
    !! NOT YET IMPLEMENTED.
    logical, intent(in), optional :: mask
    !! Optional. Boolean whether to add a mask for the atom.

    ! Local variables
    integer :: j
    !! Loop index.
    integer :: idx
    !! The index of the species in the basis.
    integer :: length
    !! The dimension of the basis atom positions.
    logical :: mask_
    !! Boolean mask for the atom.
    integer, dimension(:), allocatable :: atom_idx
    !! Temporary array.
    logical, dimension(:), allocatable :: atom_mask
    !! Temporary array.
    real(real32), dimension(:,:), allocatable :: positions
    !! Temporary array.
    type(species_type), dimension(:), allocatable :: species_list
    !! Temporary array.


    mask_ = .true.
    if(present(mask)) mask_ = mask

    this%natom = this%natom + 1
    length = size(this%spec(1)%atom,dim=1)
    idx = findloc(this%spec(:)%name, strip_null(species), dim=1)
    if(idx.eq.0)then
       this%nspec = this%nspec + 1
       allocate(species_list(this%nspec))
       species_list(1:this%nspec-1) = this%spec(1:this%nspec-1)
       deallocate(this%spec)
       species_list(this%nspec)%name = strip_null(species)
       species_list(this%nspec)%num = 1
       call get_element_properties(species_list(this%nspec)%name, &
            species_list(this%nspec)%mass, &
            species_list(this%nspec)%charge, &
            species_list(this%nspec)%radius &
       )
       allocate(species_list(this%nspec)%atom_idx(1))
       allocate(species_list(this%nspec)%atom_mask(1), source = mask_)
       species_list(this%nspec)%atom_idx(1) = this%natom
       allocate(species_list(this%nspec)%atom(length,1))
       species_list(this%nspec)%atom(:,1) = 0._real32
       species_list(this%nspec)%atom(:3,1) = position
       this%spec = species_list
       deallocate(species_list)
    else
       allocate(atom_mask(this%spec(idx)%num+1), source = .true.)
       if(allocated(this%spec(idx)%atom_mask))then
          atom_mask(1:this%spec(idx)%num) = this%spec(idx)%atom_mask
       end if
       atom_mask(this%spec(idx)%num+1) = mask_
       allocate(atom_idx(this%spec(idx)%num+1))
       if(allocated(this%spec(idx)%atom_idx))then
          atom_idx(1:this%spec(idx)%num) = this%spec(idx)%atom_idx
       else
          atom_idx(1:this%spec(idx)%num) = [ ( j, j = 1, this%spec(idx)%num ) ]
       end if
       atom_idx(this%spec(idx)%num+1) = this%natom
       allocate(positions(length,this%spec(idx)%num+1))
       positions = 0._real32
       positions(:,1:this%spec(idx)%num) = this%spec(idx)%atom
       positions(:3,this%spec(idx)%num+1) = position
       this%spec(idx)%num = this%spec(idx)%num + 1
       this%spec(idx)%atom_mask = atom_mask
       this%spec(idx)%atom_idx = atom_idx
       this%spec(idx)%atom = positions
       deallocate(atom_mask)
       deallocate(atom_idx)
       deallocate(positions)
    end if

  end subroutine add_atom
!###############################################################################


!###############################################################################
  subroutine remove_atom(this, ispec, iatom)
    !! Remove an atom from the basis.
    implicit none

    ! Arguments
    class(basis_type), intent(inout) :: this
    !! Parent. The basis.
    integer, intent(in) :: ispec, iatom
    !! The species and atom to remove.

    ! Local variables
    integer :: i
    !! Loop index.
    integer :: remove_idx
    !! The index associated with the atom to remove.
    integer, dimension(:), allocatable :: atom_idx
    !! Temporary array to store the atomic indices.
    logical, dimension(:), allocatable :: atom_mask
    !! Temporary array to store the atomic masks.
    real(real32), dimension(:,:), allocatable :: atom
    !! Temporary array to store the atomic positions.


    !---------------------------------------------------------------------------
    ! remove atom from basis
    !---------------------------------------------------------------------------
    remove_idx = this%spec(ispec)%atom_idx(iatom)
    do i = 1, this%nspec
       if(i.eq.ispec)then
          if(iatom.gt.this%spec(i)%num)then
             call stop_program("Atom to remove does not exist")
             return
          end if
          allocate(atom_mask(this%spec(i)%num-1), source = .true.)
          allocate(atom_idx(this%spec(i)%num-1))
          allocate(atom(size(this%spec(i)%atom,1),this%spec(i)%num-1))
          if(iatom.eq.1)then
             atom_mask(1:this%spec(i)%num-1) = &
                  this%spec(i)%atom_mask(2:this%spec(i)%num:1)
             atom_idx(1:this%spec(i)%num-1) = &
                  this%spec(i)%atom_idx(2:this%spec(i)%num:1)
             atom(:,1:this%spec(i)%num-1:1) = &
                  this%spec(i)%atom(:,2:this%spec(i)%num:1)
          elseif(iatom.eq.this%spec(i)%num)then
             atom_mask(1:this%spec(i)%num-1) = &
                  this%spec(i)%atom_mask(1:this%spec(i)%num-1:1)
             atom_idx(1:this%spec(i)%num-1) = &
                  this%spec(i)%atom_idx(1:this%spec(i)%num-1:1)
             atom(:,1:this%spec(i)%num-1:1) = &
                  this%spec(i)%atom(:,1:this%spec(i)%num-1:1)
          else
             atom_mask(1:iatom-1:1) = this%spec(i)%atom_mask(1:iatom-1:1)
             atom_idx(1:iatom-1:1) = this%spec(i)%atom_idx(1:iatom-1:1)
             atom_idx(iatom:this%spec(i)%num-1:1) = &
                  this%spec(i)%atom_idx(iatom+1:this%spec(i)%num:1)
             atom(:,1:iatom-1:1) = this%spec(i)%atom(:,1:iatom-1:1)
             atom(:,iatom:this%spec(i)%num-1:1) = &
                  this%spec(i)%atom(:,iatom+1:this%spec(i)%num:1)
          end if
          where(atom_idx(1:this%spec(i)%num-1:1).gt.remove_idx)
             atom_idx(1:this%spec(i)%num-1:1) = &
                  atom_idx(1:this%spec(i)%num-1:1) - 1
          end where
          this%spec(i)%atom_mask = atom_mask
          this%spec(i)%atom_idx = atom_idx
          this%spec(i)%atom = atom
          deallocate(atom_mask)
          deallocate(atom_idx)
          deallocate(atom)
          this%spec(i)%num = this%spec(i)%num - 1
          this%natom = this%natom - 1
          if(this%spec(i)%num.eq.0)then
             deallocate(this%spec(i)%atom)
             if(this%nspec.eq.0)then
                deallocate(this%spec)
                this%lcart = .true.
                this%sysname = ""
                this%energy = 0._real32
                this%lat = 0._real32
                this%pbc = .true.
             end if
          end if
       end if
    end do

  end subroutine remove_atom
!###############################################################################


!###############################################################################
  subroutine remove_atoms(this, atoms)
    !! Remove atoms from the basis.
    use coreutils, only: swap
    implicit none

    ! Arguments
    class(basis_type), intent(inout) :: this
    !! Parent. The basis.
    integer, dimension(:,:), intent(in) :: atoms
    !! The atoms to remove (2, number of atoms to remove)
    !! 1st value of 1st dimension is the species number
    !! 2nd value of 1st dimension is the atom number
    !! 2nd dimension is the number of atoms to remove

    ! Local variables
    integer :: is, ia, i
    !! Loop index.
    integer :: n, m, start_idx, end_idx, loc
    !! Index variables.
    integer :: num_species
    !! The number of species.
    integer, dimension(:,:), allocatable :: atoms_ordered
    !! The atoms to remove ordered by species and atom


    !---------------------------------------------------------------------------
    ! reorder atoms to remove
    !---------------------------------------------------------------------------
    allocate(atoms_ordered, source=atoms)
    n = size(atoms_ordered, 1)
    m = size(atoms_ordered, 2)

    do i = 1, m
       loc = maxloc(atoms_ordered(1, i:n), dim=1) + i - 1
       if (loc .ne. i) then
          call swap(atoms_ordered(1, i), atoms_ordered(1, loc))
          call swap(atoms_ordered(2, i), atoms_ordered(2, loc))
       end if
    end do
    num_species = this%nspec
    do is = 1, num_species
       start_idx = findloc(atoms_ordered(1, :), is, dim=1)
       end_idx   = findloc(atoms_ordered(1, :), is, dim=1, back=.true.)
       if (start_idx .eq. 0) cycle
       do ia = start_idx, end_idx, 1
          loc = maxloc( &
               atoms_ordered(2, ia:end_idx), &
               dim=1 &
          ) + ia - 1
          if (loc .ne. ia) then
             call swap(atoms_ordered(1, ia), atoms_ordered(1, loc))
             call swap(atoms_ordered(2, ia), atoms_ordered(2, loc))
          end if
       end do
    end do


    !---------------------------------------------------------------------------
    ! remove atoms from basis
    !---------------------------------------------------------------------------
    do i = 1, size(atoms_ordered, 2)
       call this%remove_atom(atoms_ordered(1, i), atoms_ordered(2, i))
    end do

    do is = 1, this%nspec
       if (this%spec(is)%num .eq. 0) then
          this%spec = [ this%spec(1:is-1), this%spec(is+1:) ]
          this%nspec = this%nspec - 1
       end if
    end do

  end subroutine remove_atoms
!###############################################################################


!###############################################################################
  subroutine set_element_properties_to_default(this)
    !! Set the element properties to default values.
    implicit none

    ! Arguments
    class(basis_type), intent(inout) :: this
    !! Parent. The basis.

    ! Local variables
    integer :: i
    !! Loop index.

    do i = 1, this%nspec
       call get_element_properties( &
            element = this%spec(i)%name, &
            charge = this%spec(i)%charge, &
            mass = this%spec(i)%mass, &
            radius = this%spec(i)%radius &
       )
    end do

  end subroutine set_element_properties_to_default
!###############################################################################

end module atomstruc__types
