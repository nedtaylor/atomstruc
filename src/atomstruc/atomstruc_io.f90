module atomstruc__io
  !! Module to store, read and write geometry files
  !!
  !! This module contains the procedures to read and write geometry files.
  !! It also contains the derived types used to store the geometry data.
  use coreutils, only: real32, pi, stop_program, print_warning, &
       to_upper, to_lower, jump, icount
  use atomstruc__types, only: basis_type, species_type
  use atomstruc__elements, only: get_element_properties
  implicit none


  private

  public :: igeom_input, igeom_output
  public :: geom_read, geom_write


  integer :: igeom_input = 1
  !! geometry input file format
  !! 1 = VASP
  !! 2 = CASTEP
  !! 3 = Quantum Espresso
  !! 4 = CRYSTAL
  !! 5 = XYZ
  !! 6 = extended XYZ
  integer :: igeom_output = 1
  !! geometry output file format




contains

!###############################################################################
  subroutine geom_read(UNIT, basis, length, iostat)
    !! Read geometry from a file.
    implicit none

    ! Arguments
    integer, intent(in) :: UNIT
    !! The unit number of the file to read from.
    type(basis_type), intent(out) :: basis
    !! The basis to read the geometry into.
    integer, optional, intent(in) :: length
    !! Optional. The dimension of the basis atom positions.
    integer, optional, intent(out) :: iostat
    !! Optional. The I/O status of the read.

    ! Local variables
    integer :: i
    !! Loop index.
    integer :: length_
    !! The dimension of the basis atom positions.
    integer :: iostat_
    !! The I/O status of the read.


    length_ = 3
    iostat_ = 0
    if(present(length)) length_=length

    select case(igeom_input)
    case(1)
       call VASP_geom_read(UNIT, basis, length_, iostat_)
    case(2)
       call CASTEP_geom_read(UNIT, basis, length_)
    case(3)
       call QE_geom_read(UNIT, basis, length_)
    case(4)
       call stop_program("Not yet set up for CRYSTAL")
       return
    case(5)
       call XYZ_geom_read(UNIT, basis, length_, iostat_)
       call print_warning("XYZ file format does not contain lattice data")
    case(6)
       call extXYZ_geom_read(UNIT, basis, length_, iostat_)
    end select
    if(iostat_.ne.0) then
       if(present(iostat)) iostat = iostat_
       return
    else
       if(present(iostat)) iostat = 0
    end if
    if(length_.eq.4)then
       do i = 1, basis%nspec
          basis%spec(i)%atom(4,:) = 1._real32
       end do
    end if
    do i = 1, basis%nspec
       call get_element_properties( &
            basis%spec(i)%name, &
            mass = basis%spec(i)%mass, &
            charge = basis%spec(i)%charge, &
            radius = basis%spec(i)%radius )
    end do

  end subroutine geom_read
!###############################################################################


!###############################################################################
  subroutine geom_write(UNIT, basis)
    !! Write geometry to a file.
    implicit none

    ! Arguments
    integer, intent(in) :: UNIT
    !! The unit number of the file to write to.
    class(basis_type), intent(in) :: basis
    !! The basis to write the geometry from.

    ! MAKE IT CHANGE HERE IF USER SPECIFIES LCART OR NOT
    ! AND GIVE IT THE CASTEP AND QE OPTION OF LABC !

    select case(igeom_output)
    case(1)
       call VASP_geom_write(UNIT,basis)
    case(2)
       call CASTEP_geom_write(UNIT,basis)
    case(3)
       call QE_geom_write(UNIT,basis)
    case(4)
       call stop_program("ERROR: Not yet set up for CRYSTAL")
       return
    case(5)
       call XYZ_geom_write(UNIT,basis)
    case(6)
       call extXYZ_geom_write(UNIT,basis)
    end select

  end subroutine geom_write
!###############################################################################


!###############################################################################
  subroutine VASP_geom_read(UNIT, basis, length, iostat)
    !! Read the structure in vasp poscar style format.
    implicit none

    ! Arguments
    integer, intent(in) :: UNIT
    !! The unit number of the file to read from.
    type(basis_type), intent(out) :: basis
    !! The basis to read the geometry into.
    integer, intent(in), optional :: length
    !! Optional. The dimension of the basis atom positions.
    integer, intent(out), optional :: iostat
    !! Optional. The I/O status of the read.

    integer :: Reason
    !! The I/O status of the read.
    integer :: pos, count, natom
    !! Temporary integer variables.
    real(real32) :: scal
    !! The scaling factor of the lattice.
    character(len=100) :: lspec
    !! The species names and number of each atomic species.
    character(len=1024) :: buffer
    !! Temporary character variable.
    integer :: i, j, k
    !! Loop index.
    integer :: length_
    !! The dimension of the basis atom positions.
    integer :: iostat_
    !! The I/O status of the read.


    length_ = 3
    iostat_ = 0
    !---------------------------------------------------------------------------
    ! determine dimension of basis (include translation dimension for symmetry?)
    !---------------------------------------------------------------------------
    if(present(length)) length_ = length


    !---------------------------------------------------------------------------
    ! read system name
    !---------------------------------------------------------------------------
    read(UNIT,'(A)',iostat=Reason) basis%sysname
    if(Reason.ne.0)then
       write(0,'("ERROR: The file is not in POSCAR format.")')
       write(0,*) "Expected system name, got: ",trim(basis%sysname)
       iostat_ = 1
       if(present(iostat)) iostat = iostat_
       return
    end if
    read(UNIT,*) scal


    !---------------------------------------------------------------------------
    ! read lattice
    !---------------------------------------------------------------------------
    do i = 1, 3
       read(UNIT,*) (basis%lat(i,j),j=1,3)
    end do
    basis%lat=scal*basis%lat


    !---------------------------------------------------------------------------
    ! read species names and number of each atomic species
    !---------------------------------------------------------------------------
    read(UNIT,'(A)') lspec
    basis%nspec = icount(lspec)
    allocate(basis%spec(basis%nspec))
    if(verify(lspec,' 0123456789').ne.0) then
       count=0;pos=1
       speccount: do
          i=verify(lspec(pos:), ' ')
          if (i.eq.0) exit speccount
          count=count+1
          pos=i+pos-1
          i=scan(lspec(pos:), ' ')
          if (i.eq.0) exit speccount
          basis%spec(count)%name=lspec(pos:pos+i-1)
          pos=i+pos-1
       end do speccount

       read(UNIT,*) (basis%spec(j)%num,j=1,basis%nspec)
    else !only numbers
       do count = 1, basis%nspec
          write(basis%spec(count)%name,'(I0)') count
       end do
       read(lspec,*) (basis%spec(j)%num,j=1,basis%nspec)
    end if


    !---------------------------------------------------------------------------
    ! determines whether input basis is in direct or cartesian coordinates
    !---------------------------------------------------------------------------
    basis%lcart=.false.
    read(UNIT,'(A)') buffer
    buffer = to_lower(buffer)
    if(verify(trim(buffer),'direct').eq.0) basis%lcart=.false.
    if(verify(trim(buffer),'cartesian').eq.0) basis%lcart=.true.


    !---------------------------------------------------------------------------
    ! read basis
    !---------------------------------------------------------------------------
    natom = 0
    do i = 1, basis%nspec
       allocate(basis%spec(i)%atom_idx(basis%spec(i)%num))
       allocate(basis%spec(i)%atom_mask(basis%spec(i)%num), source = .true.)
       allocate(basis%spec(i)%atom(length_,basis%spec(i)%num))
       basis%spec(i)%atom(:,:) = 0._real32
       do j = 1, basis%spec(i)%num
          natom = natom + 1
          basis%spec(i)%atom_idx(j) = natom
          read(UNIT,*) (basis%spec(i)%atom(k,j),k=1,3)
       end do
    end do


    !---------------------------------------------------------------------------
    ! convert basis if in cartesian coordinates
    !---------------------------------------------------------------------------
    if(basis%lcart) call basis%convert()


    !---------------------------------------------------------------------------
    ! normalise basis to between 0 and 1 in direct coordinates
    !---------------------------------------------------------------------------
    do i = 1, basis%nspec
       do j = 1, basis%spec(i)%num
          do k = 1, 3
             basis%spec(i)%atom(k,j) = &
                  basis%spec(i)%atom(k,j) - floor(basis%spec(i)%atom(k,j))
          end do
       end do
    end do
    basis%natom = sum(basis%spec(:)%num)

    if(present(iostat)) iostat = iostat_

  end subroutine VASP_geom_read
!###############################################################################


!###############################################################################
  subroutine VASP_geom_write(UNIT, basis, cartesian)
    !! Write the structure in vasp poscar style format.
    implicit none

    ! Arguments
    integer, intent(in) :: UNIT
    !! The unit number of the file to write to.
    class(basis_type), intent(in) :: basis
    !! The basis to write the geometry from.
    logical, intent(in), optional :: cartesian
    !! Optional. Whether to write the basis in cartesian coordinates.

    ! Local variables
    integer :: i,j
    !! Loop index.
    character(100) :: fmt
    !! Format string.
    character(10) :: string
    !! String to determine whether to write in direct or cartesian coordinates.


    string="Direct"
    if(present(cartesian))then
       if(cartesian) string="Cartesian"
    end if

    write(UNIT,'(A)') trim(adjustl(basis%sysname))
    write(UNIT,'(F15.9)') 1._real32
    do i = 1, 3
       write(UNIT,'(3(F15.9))') basis%lat(i,:)
    end do
    write(fmt,'("(",I0,"(A,1X))")') basis%nspec
    write(UNIT,trim(adjustl(fmt))) (adjustl(basis%spec(j)%name),j=1,basis%nspec)
    write(fmt,'("(",I0,"(I0,5X))")') basis%nspec
    write(UNIT,trim(adjustl(fmt))) (basis%spec(j)%num,j=1,basis%nspec)
    write(UNIT,'(A)') trim(adjustl(string))
    do i = 1, basis%nspec
       do j = 1, basis%spec(i)%num
          write(UNIT,'(3(F15.9))') basis%spec(i)%atom(1:3,j)
       end do
    end do

  end subroutine VASP_geom_write
!###############################################################################


!###############################################################################
  subroutine QE_geom_read(UNIT,basis,length)
    !! Read the structure in Quantum Espresso style format.
    implicit none

    ! Arguments
    integer, intent(in) :: UNIT
    !! The unit number of the file to read from.
    type(basis_type), intent(out) :: basis
    !! The basis to read the geometry into.
    integer, intent(in), optional :: length
    !! Optional. The dimension of the basis atom positions.

    ! Local variables
    integer :: Reason
    !! The I/O status of the read.
    integer :: i, j, k, iline
    !! Loop index.
    integer :: length_ = 3
    !! The dimension of the basis atom positions.
    integer, dimension(1000) :: tmp_natom
    !! Temporary array to store the number of atoms of each species.
    real(real32), dimension(3) :: tmpvec
    !! Temporary array to store the atomic positions.
    character(len=3) :: ctmp
    !! Temporary character variable.
    character(256) :: stop_msg
    !! Error message.
    character(len=3), dimension(1000) :: tmp_spec
    !! Temporary array to store the species names.
    character(len=1024) :: buffer, buffer2
    !! Temporary character variables.


    !---------------------------------------------------------------------------
    ! determine dimension of basis (include translation dimension for symmetry?)
    !---------------------------------------------------------------------------
    if(present(length)) length_ = length
    basis%lcart = .false.
    basis%sysname = "Converted_from_geom_file"


    !---------------------------------------------------------------------------
    ! read lattice
    !---------------------------------------------------------------------------
    rewind UNIT
    cellparam: do
       read(UNIT,'(A)',iostat=Reason) buffer
       if(Reason.ne.0)then
          call stop_program( &
               "An issue with the QE input file format has been encountered." &
          )
          return
       end if
       if(index(trim(buffer),"ibrav").ne.0)then
          write(stop_msg,*) &
               "Internal error in QE_geom_read" // &
               achar(13) // achar(10) // &
               "  Subroutine not yet set up to read IBRAV lattices"
          call stop_program(stop_msg)
          return
       end if
       if(verify("CELL_PARAMETERS",buffer).eq.0) then
          exit cellparam
       end if
    end do cellparam
    do i = 1, 3
       read(UNIT,*) (basis%lat(i,j),j=1,3)
    end do


    !---------------------------------------------------------------------------
    ! determines whether input basis is in direct or cartesian coordinates
    !---------------------------------------------------------------------------
    iline=0
    rewind UNIT
    basfind: do
       read(UNIT,'(A)',iostat=Reason) buffer
       iline=iline+1
       if(verify("ATOMIC_POSITIONS",buffer).eq.0)then
          backspace(UNIT)
          read(UNIT,*) buffer,buffer2
          if(verify("crystal",buffer2).eq.0) basis%lcart = .false.
          if(verify("angstrom",buffer2).eq.0) basis%lcart = .true.
          exit basfind
       end if
    end do basfind


    !---------------------------------------------------------------------------
    ! read basis
    !---------------------------------------------------------------------------
    basis%natom = 0
    basis%nspec = 0
    tmp_natom   = 1
    basread: do
       read(UNIT,'(A)',iostat=Reason) buffer
       read(buffer,*) ctmp
       if(Reason.ne.0) exit
       if(trim(ctmp).eq.'') exit
       if(verify(buffer,' 0123456789').eq.0) exit
       basis%natom = basis%natom + 1
       if(.not.any(tmp_spec(1:basis%nspec).eq.ctmp))then
          basis%nspec = basis%nspec + 1
          tmp_spec(basis%nspec) = ctmp
       else
          where(tmp_spec(1:basis%nspec).eq.ctmp)
             tmp_natom(1:basis%nspec) = tmp_natom(1:basis%nspec) + 1
          end where
       end if
    end do basread

    allocate(basis%spec(basis%nspec))
    basis%spec(1:basis%nspec)%name = tmp_spec(1:basis%nspec)
    do i = 1, basis%nspec
       basis%spec(i)%num = 0
       allocate(basis%spec(i)%atom_idx(tmp_natom(i)))
       allocate(basis%spec(i)%atom_mask(tmp_natom(i)), source = .true.)
       allocate(basis%spec(i)%atom(tmp_natom(i),length_))
    end do

    call jump(UNIT,iline)
    basread2: do i = 1, basis%natom
       read(UNIT,*,iostat=Reason) ctmp,tmpvec(1:3)
       do j = 1, basis%nspec
          if(basis%spec(j)%name.eq.ctmp)then
             basis%spec(j)%num = basis%spec(j)%num + 1
             basis%spec(j)%atom_idx(basis%spec(j)%num) = i
             basis%spec(j)%atom(1:3,basis%spec(j)%num) = tmpvec(1:3)
             exit
          end if
       end do
    end do basread2


    !---------------------------------------------------------------------------
    ! convert basis if in cartesian coordinates
    !---------------------------------------------------------------------------
    if(basis%lcart) call basis%convert()


    !---------------------------------------------------------------------------
    ! normalise basis to between 0 and 1 in direct coordinates
    !---------------------------------------------------------------------------
    do i = 1, basis%nspec
       do j = 1, basis%spec(i)%num
          do k = 1, 3
             basis%spec(i)%atom(k,j) = &
                  basis%spec(i)%atom(k,j) - floor( basis%spec(i)%atom(k,j) )
          end do
       end do
    end do
    basis%natom = sum(basis%spec(:)%num)

  end subroutine QE_geom_read
!###############################################################################


!###############################################################################
  subroutine QE_geom_write(UNIT, basis, cartesian)
    !! Write the structure in Quantum Espresso style format.
    implicit none

    ! Arguments
    integer, intent(in) :: UNIT
    !! The unit number of the file to write to.
    class(basis_type), intent(in) :: basis
    !! The basis to write the geometry from.
    logical, intent(in), optional :: cartesian
    !! Optional. Whether to write the basis in cartesian coordinates.

    ! Local variables
    integer :: i,j
    !! Loop index.
    character(10) :: string
    !! String to determine whether to write in crystal or angstrom coordinates.


    string="crystal"
    if(present(cartesian))then
       if(cartesian) string="angstrom"
    end if


    write(UNIT,'("CELL_PARAMETERS angstrom")')
    do i = 1, 3
       write(UNIT,'(3(F15.9))') basis%lat(i,:)
    end do
    write(UNIT,'("ATOMIC_SPECIES")')
    do i = 1, basis%nspec
       write(UNIT,'(A)') trim(adjustl(basis%spec(i)%name))
    end do
    write(UNIT,'("ATOMIC_POSITIONS",1X,A)') trim(adjustl(string))
    do i = 1, basis%nspec
       do j = 1, basis%spec(i)%num
          write(UNIT,'(A5,1X,3(F15.9))') &
               basis%spec(i)%name,basis%spec(i)%atom(1:3,j)
       end do
    end do

  end subroutine QE_geom_write
!###############################################################################


!###############################################################################
  subroutine CASTEP_geom_read(UNIT, basis, length)
    !! Read the structure in CASTEP style format.
    implicit none

    ! Arguments
    integer, intent(in) :: UNIT
    !! The unit number of the file to read from.
    type(basis_type), intent(out) :: basis
    !! The basis to read the geometry into.
    integer, intent(in), optional :: length
    !! Optional. The dimension of the basis atom positions.

    ! Local variables
    integer :: Reason
    !! The I/O status of the read.
    integer :: i, j, k, iline
    !! Loop index.
    integer :: length_ = 3
    !! The dimension of the basis atom positions.
    integer :: itmp1
    !! Temporary integer variable.
    character(len=3) :: ctmp
    !! Temporary character variable.
    character(len=20) :: units
    !! Units of the lattice vectors.
    character(len=200) :: buffer, store
    !! Temporary character variables.
    logical :: labc
    !! Logical variable to determine whether the lattice is in abc or
    !! cartesian coordinates.
    integer, dimension(1000) :: tmp_natom
    !! Temporary array to store the number of atoms of each species.
    real(real32), dimension(3) :: abc, angle, dvtmp1
    !! Temporary arrays to store the lattice vectors.
    character(len=3), dimension(1000) :: tmp_spec
    !! Temporary array to store the species names.


    !---------------------------------------------------------------------------
    ! determine dimension of basis (include translation dimension for symmetry?)
    !---------------------------------------------------------------------------
    if(present(length)) length_ = length


    !---------------------------------------------------------------------------
    ! reading loop of file
    !---------------------------------------------------------------------------
    tmp_spec = ""
    tmp_natom = 0
    iline = 0
    labc = .true.
    basis%sysname = "from CASTEP"
    rewind(UNIT)
    readloop: do
       iline=iline+1
       read(UNIT,'(A)',iostat=Reason) buffer
       if(Reason.ne.0) exit
       buffer=to_upper(buffer)
       if(scan(trim(adjustl(buffer)),'%').ne.1) cycle readloop
       if(index(trim(adjustl(buffer)),'%END').eq.1) cycle readloop
       read(buffer,*) store, buffer
       if(trim(buffer).eq.'') cycle readloop
       !------------------------------------------------------------------------
       ! read lattice
       !------------------------------------------------------------------------
       lattice_if: if(index(trim(buffer),"LATTICE").eq.1)then
          if(index(trim(buffer),"ABC").ne.0) labc = .true.
          if(index(trim(buffer),"CART").ne.0) labc = .false.
          store = ""
          itmp1 = 0
          lattice_loop: do
             itmp1 = itmp1 + 1
             read(UNIT,'(A)',iostat=Reason) buffer
             if(Reason.ne.0) exit lattice_loop
             if(scan(trim(adjustl(buffer)),'%').eq.1) exit lattice_loop
             if(itmp1.eq.5)then
                call stop_program( &
                     "Too many lines in LATTICE block of structure file" &
                )
                return
             end if
             store=trim(store)//" "//trim(buffer)
          end do lattice_loop
          iline=iline+itmp1

          if(labc)then
             read(store,*) units,(abc(i),i=1,3), (angle(j),j=1,3)
             basis%lat = convert_abc_to_lat(abc,angle,.false.)
          else
             read(store,*) units,(basis%lat(i,:),i=1,3)
          end if
          cycle readloop
       end if lattice_if

       !------------------------------------------------------------------------
       ! read basis
       !------------------------------------------------------------------------
       basis_if: if(index(trim(buffer),"POSITIONS").eq.1) then
          if(index(trim(buffer),"ABS").ne.0) basis%lcart=.true.
          if(index(trim(buffer),"FRAC").ne.0) basis%lcart=.false.
          itmp1 = 0
          basis_loop1: do
             read(UNIT,'(A)',iostat=Reason) buffer
             if(Reason.ne.0) exit basis_loop1
             if(scan(trim(adjustl(buffer)),'%').eq.1) exit basis_loop1
             read(buffer,*) ctmp
             if(trim(ctmp).eq.'') exit
             if(verify(buffer,' 0123456789').eq.0) exit
             basis%natom = basis%natom + 1
             if(.not.any(tmp_spec(1:basis%nspec).eq.ctmp))then
                basis%nspec = basis%nspec+1
                tmp_natom(basis%nspec) = 1
                tmp_spec(basis%nspec)  = ctmp
             else
                where(tmp_spec(1:basis%nspec).eq.ctmp)
                   tmp_natom(1:basis%nspec) = tmp_natom(1:basis%nspec) + 1
                end where
             end if
          end do basis_loop1

          allocate(basis%spec(basis%nspec))
          basis%spec(1:basis%nspec)%name = tmp_spec(1:basis%nspec)
          do i = 1, basis%nspec
             basis%spec(i)%num = 0
             allocate(basis%spec(i)%atom(length_,tmp_natom(i)))
          end do

          call jump(UNIT,iline)
          basis_loop2: do i = 1, basis%natom
             read(UNIT,'(A)',iostat=Reason) buffer
             if(Reason.ne.0)then
                call stop_program("Internal error in assigning the basis")
                return
             end if
             read(buffer,*) ctmp,dvtmp1(1:3)
             species_loop: do j = 1, basis%nspec
                if(basis%spec(j)%name.eq.ctmp)then
                   basis%spec(j)%num = basis%spec(j)%num + 1
                   basis%spec(j)%atom(1:3,basis%spec(j)%num) = dvtmp1(1:3)
                   exit species_loop
                end if
             end do species_loop
          end do basis_loop2

       end if basis_if
    end do readloop


    !---------------------------------------------------------------------------
    ! convert basis if in cartesian coordinates
    !---------------------------------------------------------------------------
    if(basis%lcart) call basis%convert()


    !---------------------------------------------------------------------------
    ! normalise basis to between 0 and 1 in direct coordinates
    !---------------------------------------------------------------------------
    do i = 1, basis%nspec
       do j = 1, basis%spec(i)%num
          do k = 1, 3
             basis%spec(i)%atom(k,j) = &
                  basis%spec(i)%atom(k,j) - floor( basis%spec(i)%atom(k,j) )
          end do
       end do
    end do
    basis%natom = sum(basis%spec(:)%num)

  end subroutine CASTEP_geom_read
!###############################################################################


!###############################################################################
  subroutine CASTEP_geom_write(UNIT, basis, labc, cartesian)
    !! Write the structure in CASTEP style format.
    implicit none

    ! Arguments
    integer :: UNIT
    !! The unit number of the file to write to.
    class(basis_type), intent(in) :: basis
    !! The basis to write the geometry from.
    logical, intent(in), optional :: labc
    !! Optional. Boolean whether to write the lattice in abc format.
    logical, intent(in), optional :: cartesian
    !! Optional. Boolean whether to write basis in cartesian coordinates.

    ! Local variables
    integer :: i, j
    !! Loop index.
    real(real32), dimension(2,3) :: abc_angle
    !! Temporary arrays to store the lattice vectors.
    character(4) :: string_lat, string_bas
    !! Strings specifying lattice and basis format
    character(len=256) :: stop_msg
    !! Error message.


    string_lat="CART"
    if(present(labc))then
       if(labc) string_lat="ABC"
    end if

    string_bas="FRAC"
    if(present(cartesian))then
       if(cartesian)then
          string_bas="ABS"
          write(stop_msg,*) &
               "Internal error in CASTEP_geom_write" // &
               achar(13) // achar(10) // &
               "  Subroutine not yet set up to output cartesian coordinates"
          call stop_program(stop_msg)
          return
       end if
    end if

    write(UNIT,'("%block LATTICE_",A)') trim(string_lat)
    write(UNIT,'("ang")')
    if(present(labc))then
       if(labc)then
          abc_angle = convert_lat_to_abc(basis%lat)
          write(UNIT,'(3(F15.9))') abc_angle(1,:)
          write(UNIT,'(3(F15.9))') abc_angle(2,:)
          goto 10
       end if
    end if
    do i = 1, 3
       write(UNIT,'(3(F15.9))') basis%lat(i,:)
    end do

10  write(UNIT,'("%endblock LATTICE_",A)') trim(string_lat)

    write(UNIT,*)
    write(UNIT,'("%block POSITIONS_",A)') trim(string_bas)
    do i = 1, basis%nspec
       do j = 1, basis%spec(i)%num
          write(UNIT,'(A5,1X,3(F15.9))') &
               basis%spec(i)%name,basis%spec(i)%atom(1:3,j)
       end do
    end do
    write(UNIT,'("%endblock POSITIONS_",A)') trim(string_bas)

  end subroutine CASTEP_geom_write
!###############################################################################


!###############################################################################
  subroutine XYZ_geom_read(UNIT, basis, length, iostat)
    !! Read the structure in xyz style format.
    implicit none

    ! Arguments
    integer, intent(in) :: UNIT
    !! The unit number of the file to read from.
    type(basis_type), intent(out) :: basis
    !! The basis to read the geometry into.
    integer, intent(in), optional :: length
    !! Optional. The dimension of the basis atom positions.
    integer, intent(out), optional :: iostat
    !! Optional. The I/O status of the read.

    ! Local variables
    integer :: Reason
    !! The I/O status of the read.
    integer :: i, j
    !! Loop index.
    integer, allocatable, dimension(:) :: tmp_num
    !! Temporary array to store the number of atoms of each species.
    real(real32), dimension(3) :: vec
    !! Temporary array to store the atomic positions.
    real(real32), allocatable, dimension(:,:,:) :: tmp_bas
    !! Temporary array to store the atomic positions.
    character(len=3) :: ctmp
    !! Temporary character variable.
    character(len=3), allocatable, dimension(:) :: tmp_spec
    !! Temporary array to store the species names.
    integer :: length_
    !! The dimension of the basis atom positions.
    integer :: iostat_
    !! The I/O status of the read.


    length_ = 3
    iostat_ = 0
    if(present(length)) length_ = length


    read(UNIT,*,iostat=Reason) basis%natom
    if(Reason.ne.0)then
       write(0,'("ERROR: The file is not in xyz format.")')
       iostat_ = 1
       if(present(iostat)) iostat = iostat_
       return
    end if
    read(UNIT,'(A)',iostat=Reason) basis%sysname


    !---------------------------------------------------------------------------
    ! read basis
    !---------------------------------------------------------------------------
    allocate(tmp_spec(basis%natom))
    allocate(tmp_num(basis%natom))
    allocate(tmp_bas(length_,basis%natom,basis%natom))
    tmp_num(:) = 0
    tmp_spec = ""
    tmp_bas = 0
    basis%nspec = 0
    do i = 1, basis%natom
       read(UNIT,*,iostat=Reason) ctmp, vec(1:3)
       if(.not.any(tmp_spec(1:basis%nspec).eq.ctmp))then
          basis%nspec = basis%nspec + 1
          tmp_spec(basis%nspec) = ctmp
          tmp_bas(1:3,1,basis%nspec) = vec(1:3)
          tmp_num(basis%nspec) = 1
       else
          checkspec: do j = 1, basis%nspec
             if(tmp_spec(j).eq.ctmp)then
                tmp_num(j) = tmp_num(j)+1
                tmp_bas(1:3,tmp_num(j),j) = vec(1:3)
                exit checkspec
             end if
          end do checkspec
       end if
    end do


    !---------------------------------------------------------------------------
    ! move basis from temporary basis to main basis.
    ! done to allow for correct allocation of number of and per species
    !---------------------------------------------------------------------------
    allocate(basis%spec(basis%nspec))
    do i = 1, basis%nspec
       basis%spec(i)%name = tmp_spec(i)
       basis%spec(i)%num  = tmp_num(i)
       allocate(basis%spec(i)%atom(length_,tmp_num(i)))
       basis%spec(i)%atom(:,:) = 0
       basis%spec(i)%atom(1:3,1:tmp_num(i)) = tmp_bas(1:3,1:tmp_num(i),i)
    end do

    if(present(iostat)) iostat = iostat_

  end subroutine XYZ_geom_read
!###############################################################################


!###############################################################################
  subroutine XYZ_geom_write(UNIT,basis)
    !! Write the structure in xyz style format.
    implicit none

    ! Arguments
    integer, intent(in) :: UNIT
    !! The unit number of the file to write to.
    class(basis_type), intent(in) :: basis
    !! The basis to write the geometry from.

    ! Local variables
    integer :: i, j
    !! Loop index.


    write(UNIT,'("I0")') basis%natom
    write(UNIT,'("A")') basis%sysname
    do i = 1, basis%nspec
       do j = 1, basis%spec(i)%num
          write(UNIT,'(A5,1X,3(F15.9))') &
               basis%spec(i)%name,basis%spec(i)%atom(1:3,j)
       end do
    end do

  end subroutine XYZ_geom_write
!###############################################################################


!###############################################################################
  subroutine extXYZ_geom_read(UNIT, basis, length, iostat)
    !! Read the structure in extended xyz style format.
    implicit none

    ! Arguments
    integer, intent(in) :: UNIT
    !! The unit number of the file to read from.
    type(basis_type), intent(out) :: basis
    !! The basis to read the geometry into.
    integer, intent(in), optional :: length
    !! Optional. The dimension of the basis atom positions.
    integer, intent(out), optional :: iostat
    !! Optional. The I/O status of the read.

    ! Local variables
    integer :: Reason
    !! The I/O status of the read.
    integer :: i, j
    !! Loop index.
    integer :: index1, index2
    !! Index variables.
    integer, allocatable, dimension(:) :: tmp_num
    !! Temporary array to store the number of atoms of each species.
    real(real32), dimension(3) :: vec
    !! Temporary array to store the atomic positions.
    real(real32), allocatable, dimension(:,:,:) :: tmp_bas
    !! Temporary array to store the atomic positions.
    character(len=3) :: ctmp
    !! Temporary character variable.
    character(len=3), allocatable, dimension(:) :: tmp_spec
    !! Temporary array to store the species names.
    character(len=1024) :: buffer
    !! Temporary character variable.
    integer :: length_ = 3
    !! The dimension of the basis atom positions.
    integer :: iostat_ = 0
    !! The I/O status of the read.


    basis%lcart=.true.
    if(present(length)) length_ = length


    !---------------------------------------------------------------------------
    ! read system information
    !---------------------------------------------------------------------------
    read(UNIT,*,iostat=Reason) basis%natom
    if(Reason.ne.0)then
       write(0,'("ERROR: The file is not in xyz format.")')
       iostat_ = 1
       if(present(iostat)) iostat = iostat_
       return
    end if
    read(UNIT,'(A)',iostat=Reason) buffer
    if(Reason.ne.0)then
       write(0,'("ERROR: The file is not in xyz format.")')
       iostat_ = 1
       if(present(iostat)) iostat = iostat_
       return
    end if
    index1 = index(buffer,'Lattice="') + 9
    index2 = index(buffer(index1:),'"') + index1 - 2
    read(buffer(index1:index2),*) ( ( basis%lat(i,j), j = 1, 3), i = 1, 3)

    index1 = index(buffer,'free_energy=') + 12
    read(buffer(index1:),*) basis%energy


    !---------------------------------------------------------------------------
    ! read basis
    !---------------------------------------------------------------------------
    allocate(tmp_spec(basis%natom))
    allocate(tmp_num(basis%natom))
    allocate(tmp_bas(length_,basis%natom,basis%natom))
    tmp_num(:) = 0
    tmp_spec = ""
    tmp_bas = 0
    basis%nspec = 0
    do i = 1, basis%natom
       read(UNIT,*,iostat=Reason) ctmp, vec(1:3)
       if(.not.any(tmp_spec(1:basis%nspec).eq.ctmp))then
          basis%nspec=basis%nspec+1
          tmp_spec(basis%nspec) = trim(adjustl(ctmp))
          tmp_bas(1:3,1,basis%nspec) = vec(1:3)
          tmp_num(basis%nspec) = 1
       else
          checkspec: do j = 1, basis%nspec
             if(tmp_spec(j).eq.ctmp)then
                tmp_num(j) = tmp_num(j) + 1
                tmp_bas(1:3,tmp_num(j),j) = vec(1:3)
                exit checkspec
             end if
          end do checkspec
       end if
    end do


    !---------------------------------------------------------------------------
    ! move basis from temporary basis to main basis.
    ! done to allow for correct allocation of number of and per species
    !---------------------------------------------------------------------------
    allocate(basis%spec(basis%nspec))
    basis%sysname = ""
    do i = 1, basis%nspec
       basis%spec(i)%name = tmp_spec(i)
       basis%spec(i)%num = tmp_num(i)
       allocate(basis%spec(i)%atom(length_,tmp_num(i)))
       basis%spec(i)%atom(:,:) = 0
       basis%spec(i)%atom(1:3,1:tmp_num(i)) = tmp_bas(1:3,1:tmp_num(i),i)
       write(buffer,'(I0,A)') basis%spec(i)%num,trim(basis%spec(i)%name)
       basis%sysname = basis%sysname//trim(buffer)
       if(i.lt.basis%nspec) basis%sysname = trim(adjustl(basis%sysname))//"_"
    end do

    if(present(iostat)) iostat = iostat_

  end subroutine extXYZ_geom_read
!###############################################################################


!###############################################################################
  subroutine extXYZ_geom_write(UNIT, basis)
    !! Write the structure in extended xyz style format.
    implicit none

    ! Arguments
    integer, intent(in) :: UNIT
    !! The unit number of the file to write to.
    class(basis_type), intent(in) :: basis
    !! The basis to write the geometry from.

    ! Local variables
    integer :: i, j
    !! Loop index.


    write(UNIT,'(I0)') basis%natom
    write(UNIT,'(A,8(F0.8,1X),F0.8,A)', advance="no") &
         'Lattice="',((basis%lat(i,j),j=1,3),i=1,3),'"'
    write(UNIT,'(A,F0.8)', advance="no") ' free_energy=',basis%energy
    write(UNIT,'(A)', advance="no") ' pbc="T T T"'
    if(basis%lcart)then
       do i = 1, basis%nspec
          do j = 1, basis%spec(i)%num
             write(UNIT,'(A8,3(1X, F16.8))') &
                  basis%spec(i)%name,basis%spec(i)%atom(1:3,j)
          end do
       end do
    else
       do i = 1, basis%nspec
          do j = 1, basis%spec(i)%num
             write(UNIT,'(A8,3(1X, F16.8))') basis%spec(i)%name, &
                  matmul(basis%spec(i)%atom(1:3,j),basis%lat)
          end do
       end do
    end if

  end subroutine extXYZ_geom_write
!###############################################################################



!#############################################################################
  function convert_abc_to_lat(abc,angle,radians) result(lattice)
    !! Convert the lattice from abc and αβγ to lattice matrix.
    implicit none

    ! Arguments
    real(real32), dimension(3), intent(in) :: abc, angle
    !! lattice constants
    logical, intent(in), optional :: radians
    !! Optional. Boolean whether angles are in radians.
    real(real32), dimension(3,3) :: lattice
    !! The lattice matrix.

    ! Local variables
    real(real32), dimension(3) :: in_angle
    !! The lattice angles in radians.



    in_angle = angle
    if(present(radians))then
       if(.not.radians) in_angle = angle*pi/180._real32
    end if

    lattice=0._real32

    lattice(1,1)=abc(1)
    lattice(2,:2)=(/abc(2)*cos(in_angle(3)),abc(2)*sin(in_angle(3))/)

    lattice(3,1) = abc(3)*cos(in_angle(2))
    lattice(3,2) = abc(3)*(cos(in_angle(1)) - cos(in_angle(2))*&
         cos(in_angle(3)))/sin(in_angle(3))
    lattice(3,3) = sqrt(abc(3)**2._real32 - &
         lattice(3,1)**2._real32 - &
         lattice(3,2)**2._real32)

  end function convert_abc_to_lat
!###############################################################################


!###############################################################################
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
!###############################################################################

end module atomstruc__io
