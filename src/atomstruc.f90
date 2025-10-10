module atomstruc
  !! This is the top-level module for the atomstruc Fortran library.
  use atomstruc__types, only: basis_type, species_type
  use atomstruc__io, only: geom_read, geom_write, igeom_input, igeom_output
  use atomstruc__elements, only: get_element_properties
  use atomstruc__utils, only: basis_merge
  use atomstruc__extd, only: extended_basis_type

  implicit none

  public

end module atomstruc
