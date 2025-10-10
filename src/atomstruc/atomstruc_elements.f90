module atomstruc__elements
  !! Module containing atomic element data.
  use coreutils, only: real32
  implicit none


  private

  public :: get_element_properties



contains

!###############################################################################
  subroutine get_element_properties(element, charge, mass, radius)
    !! Set the mass and charge of the element
    implicit none

    ! Arguments
    character(len=3), intent(in) :: element
    !! Element name.
    real(real32), intent(out), optional :: charge
    !! Charge of the element.
    real(real32), intent(out), optional :: mass
    !! Mass of the element.
    real(real32), intent(out), optional :: radius
    !! Radius of the element.

    ! Local variables
    real(real32) :: mass_, charge_, radius_
    !! Mass, charge and radius of the element.

    select case(element)
    case('H')
       mass_ = 1.00784_real32
       charge_ = 1.0_real32
       radius_ = 0.31_real32
    case('He')
       mass_ = 4.0026_real32
       charge_ = 2.0_real32
       radius_ = 0.28_real32
    case('Li')
       mass_ = 6.94_real32
       charge_ = 3.0_real32
       radius_ = 1.28_real32
    case('Be')
       mass_ = 9.0122_real32
       charge_ = 4.0_real32
       radius_ = 0.96_real32
    case('B')
       mass_ = 10.81_real32
       charge_ = 5.0_real32
       radius_ = 0.84_real32
    case('C')
       mass_ = 12.011_real32
       charge_ = 6.0_real32
       radius_ = 0.76_real32
    case('N')
       mass_ = 14.007_real32
       charge_ = 7.0_real32
       radius_ = 0.71_real32
    case('O')
       mass_ = 15.999_real32
       charge_ = 8.0_real32
       radius_ = 0.66_real32
    case('F')
       mass_ = 18.998_real32
       charge_ = 9.0_real32
       radius_ = 0.57_real32
    case('Ne')
       mass_ = 20.180_real32
       charge_ = 10.0_real32
       radius_ = 0.58_real32
    case('Na')
       mass_ = 22.989_real32
       charge_ = 11.0_real32
       radius_ = 1.66_real32
    case('Mg')
       mass_ = 24.305_real32
       charge_ = 12.0_real32
       radius_ = 1.41_real32
    case('Al')
       mass_ = 26.982_real32
       charge_ = 13.0_real32
       radius_ = 1.21_real32
    case('Si')
       mass_ = 28.085_real32
       charge_ = 14.0_real32
       radius_ = 1.11_real32
    case('P')
       mass_ = 30.974_real32
       charge_ = 15.0_real32
       radius_ = 1.07_real32
    case('S')
       mass_ = 32.06_real32
       charge_ = 16.0_real32
       radius_ = 1.05_real32
    case('Cl')
       mass_ = 35.453_real32
       charge_ = 17.0_real32
       radius_ = 1.02_real32
    case('Ar')
       mass_ = 39.948_real32
       charge_ = 18.0_real32
       radius_ = 1.06_real32
    case('K')
       mass_ = 39.098_real32
       charge_ = 19.0_real32
       radius_ = 2.03_real32
    case('Ca')
       mass_ = 40.078_real32
       charge_ = 20.0_real32
       radius_ = 1.74_real32
    case('Sc')
       mass_ = 44.956_real32
       charge_ = 21.0_real32
       radius_ = 1.44_real32
    case('Ti')
       mass_ = 47.867_real32
       charge_ = 22.0_real32
       radius_ = 1.32_real32
    case('V')
       mass_ = 50.942_real32
       charge_ = 23.0_real32
       radius_ = 1.22_real32
    case('Cr')
       mass_ = 51.996_real32
       charge_ = 24.0_real32
       radius_ = 1.18_real32
    case('Mn')
       mass_ = 54.938_real32
       charge_ = 25.0_real32
       radius_ = 1.17_real32
    case('Fe')
       mass_ = 55.845_real32
       charge_ = 26.0_real32
       radius_ = 1.17_real32
    case('Co')
       mass_ = 58.933_real32
       charge_ = 27.0_real32
       radius_ = 1.16_real32
    case('Ni')
       mass_ = 58.693_real32
       charge_ = 28.0_real32
       radius_ = 1.15_real32
    case('Cu')
       mass_ = 63.546_real32
       charge_ = 29.0_real32
       radius_ = 1.17_real32
    case('Zn')
       mass_ = 65.38_real32
       charge_ = 30.0_real32
       radius_ = 1.25_real32
    case('Ga')
       mass_ = 69.723_real32
       charge_ = 31.0_real32
       radius_ = 1.26_real32
    case('Ge')
       mass_ = 72.63_real32
       charge_ = 32.0_real32
       radius_ = 1.22_real32
    case('As')
       mass_ = 74.922_real32
       charge_ = 33.0_real32
       radius_ = 1.19_real32
    case('Se')
       mass_ = 78.971_real32
       charge_ = 34.0_real32
       radius_ = 1.16_real32
    case('Br')
       mass_ = 79.904_real32
       charge_ = 35.0_real32
       radius_ = 1.14_real32
    case('Kr')
       mass_ = 83.798_real32
       charge_ = 36.0_real32
       radius_ = 1.12_real32
    case('Rb')
       mass_ = 85.468_real32
       charge_ = 37.0_real32
       radius_ = 2.16_real32
    case('Sr')
       mass_ = 87.62_real32
       charge_ = 38.0_real32
       radius_ = 1.91_real32
    case('Y')
       mass_ = 88.906_real32
       charge_ = 39.0_real32
       radius_ = 1.62_real32
    case('Zr')
       mass_ = 91.224_real32
       charge_ = 40.0_real32
       radius_ = 1.45_real32
    case('Nb')
       mass_ = 92.906_real32
       charge_ = 41.0_real32
       radius_ = 1.34_real32
    case('Mo')
       mass_ = 95.95_real32
       charge_ = 42.0_real32
       radius_ = 1.3_real32
    case('Tc')
       mass_ = 98.0_real32
       charge_ = 43.0_real32
       radius_ = 1.27_real32
    case('Ru')
       mass_ = 101.07_real32
       charge_ = 44.0_real32
       radius_ = 1.25_real32
    case('Rh')
       mass_ = 102.91_real32
       charge_ = 45.0_real32
       radius_ = 1.25_real32
    case('Pd')
       mass_ = 106.42_real32
       charge_ = 46.0_real32
       radius_ = 1.28_real32
    case('Ag')
       mass_ = 107.87_real32
       charge_ = 47.0_real32
       radius_ = 1.34_real32
    case('Cd')
       mass_ = 112.41_real32
       charge_ = 48.0_real32
       radius_ = 1.48_real32
    case('In')
       mass_ = 114.82_real32
       charge_ = 49.0_real32
       radius_ = 1.44_real32
    case('Sn')
       mass_ = 118.71_real32
       charge_ = 50.0_real32
       radius_ = 1.41_real32
    case('Sb')
       mass_ = 121.76_real32
       charge_ = 51.0_real32
       radius_ = 1.38_real32
    case('Te')
       mass_ = 127.6_real32
       charge_ = 52.0_real32
       radius_ = 1.35_real32
    case('I')
       mass_ = 126.9_real32
       charge_ = 53.0_real32
       radius_ = 1.33_real32
    case('Xe')
       mass_ = 131.29_real32
       charge_ = 54.0_real32
       radius_ = 1.31_real32
    case('Cs')
       mass_ = 132.91_real32
       charge_ = 55.0_real32
       radius_ = 2.35_real32
    case('Ba')
       mass_ = 137.33_real32
       charge_ = 56.0_real32
       radius_ = 1.98_real32
    case('La')
       mass_ = 138.91_real32
       charge_ = 57.0_real32
       radius_ = 1.69_real32
    case('Ce')
       mass_ = 140.12_real32
       charge_ = 58.0_real32
       radius_ = 1.65_real32
    case('Pr')
       mass_ = 140.91_real32
       charge_ = 59.0_real32
       radius_ = 1.65_real32
    case('Nd')
       mass_ = 144.24_real32
       charge_ = 60.0_real32
       radius_ = 1.64_real32
    case('Pm')
       mass_ = 145.0_real32
       charge_ = 61.0_real32
       radius_ = 1.63_real32
    case('Sm')
       mass_ = 150.36_real32
       charge_ = 62.0_real32
       radius_ = 1.62_real32
    case('Eu')
       mass_ = 152.0_real32
       charge_ = 63.0_real32
       radius_ = 1.85_real32
    case('Gd')
       mass_ = 157.25_real32
       charge_ = 64.0_real32
       radius_ = 1.61_real32
    case('Tb')
       mass_ = 158.93_real32
       charge_ = 65.0_real32
       radius_ = 1.59_real32
    case('Dy')
       mass_ = 162.5_real32
       charge_ = 66.0_real32
       radius_ = 1.59_real32
    case('Ho')
       mass_ = 164.93_real32
       charge_ = 67.0_real32
       radius_ = 1.58_real32
    case('Er')
       mass_ = 167.26_real32
       charge_ = 68.0_real32
       radius_ = 1.57_real32
    case('Tm')
       mass_ = 168.93_real32
       charge_ = 69.0_real32
       radius_ = 1.56_real32
    case('Yb')
       mass_ = 173.05_real32
       charge_ = 70.0_real32
       radius_ = 1.74_real32
    case('Lu')
       mass_ = 174.97_real32
       charge_ = 71.0_real32
       radius_ = 1.56_real32
    case('Hf')
       mass_ = 178.49_real32
       charge_ = 72.0_real32
       radius_ = 1.44_real32
    case('Ta')
       mass_ = 180.95_real32
       charge_ = 73.0_real32
       radius_ = 1.34_real32
    case('W')
       mass_ = 183.84_real32
       charge_ = 74.0_real32
       radius_ = 1.3_real32
    case('Re')
       mass_ = 186.21_real32
       charge_ = 75.0_real32
       radius_ = 1.28_real32
    case('Os')
       mass_ = 190.23_real32
       charge_ = 76.0_real32
       radius_ = 1.26_real32
    case('Ir')
       mass_ = 192.22_real32
       charge_ = 77.0_real32
       radius_ = 1.27_real32
    case('Pt')
       mass_ = 195.08_real32
       charge_ = 78.0_real32
       radius_ = 1.3_real32
    case('Au')
       mass_ = 196.97_real32
       charge_ = 79.0_real32
       radius_ = 1.34_real32
    case('Hg')
       mass_ = 200.59_real32
       charge_ = 80.0_real32
       radius_ = 1.49_real32
    case('Tl')
       mass_ = 204.38_real32
       charge_ = 81.0_real32
       radius_ = 1.48_real32
    case('Pb')
       mass_ = 207.2_real32
       charge_ = 82.0_real32
       radius_ = 1.47_real32
    case('Bi')
       mass_ = 208.98_real32
       charge_ = 83.0_real32
       radius_ = 1.46_real32
    case('Po')
       mass_ = 209.0_real32
       charge_ = 84.0_real32
       radius_ = 1.45_real32
    case('At')
       mass_ = 210.0_real32
       charge_ = 85.0_real32
       radius_ = 1.44_real32
    case('Rn')
       mass_ = 222.0_real32
       charge_ = 86.0_real32
       radius_ = 1.43_real32
    case('Fr')
       mass_ = 223.0_real32
       charge_ = 87.0_real32
       radius_ = 2.6_real32
    case('Ra')
       mass_ = 226.0_real32
       charge_ = 88.0_real32
       radius_ = 2.21_real32
    case('Ac')
       mass_ = 227.0_real32
       charge_ = 89.0_real32
       radius_ = 1.86_real32
    case('Th')
       mass_ = 232.04_real32
       charge_ = 90.0_real32
       radius_ = 1.75_real32
    case('Pa')
       mass_ = 231.04_real32
       charge_ = 91.0_real32
       radius_ = 1.61_real32
    case('U')
       mass_ = 238.03_real32
       charge_ = 92.0_real32
       radius_ = 1.58_real32
    case('Np')
       mass_ = 237.0_real32
       charge_ = 93.0_real32
       radius_ = 1.55_real32
    case('Pu')
       mass_ = 244.0_real32
       charge_ = 94.0_real32
       radius_ = 1.53_real32
    case('Am')
       mass_ = 243.0_real32
       charge_ = 95.0_real32
       radius_ = 1.51_real32
    case('Cm')
       mass_ = 247.0_real32
       charge_ = 96.0_real32
       radius_ = 1.69_real32
    case('Bk')
       mass_ = 247.0_real32
       charge_ = 97.0_real32
       radius_ = 1.48_real32
    case('Cf')
       mass_ = 251.0_real32
       charge_ = 98.0_real32
       radius_ = 1.47_real32
    case('Es')
       mass_ = 252.0_real32
       charge_ = 99.0_real32
       radius_ = 1.46_real32
    case('Fm')
       mass_ = 257.0_real32
       charge_ = 100.0_real32
       radius_ = 1.45_real32
    case('Md')
       mass_ = 258.0_real32
       charge_ = 101.0_real32
       radius_ = 1.44_real32
    case('No')
       mass_ = 259.0_real32
       charge_ = 102.0_real32
       radius_ = 1.43_real32
    case('Lr')
       mass_ = 262.0_real32
       charge_ = 103.0_real32
       radius_ = 1.62_real32
    case('Rf')
       mass_ = 267.0_real32
       charge_ = 104.0_real32
       radius_ = 1.57_real32
    case('Db')
       mass_ = 270.0_real32
       charge_ = 105.0_real32
       radius_ = 1.49_real32
    case('Sg')
       mass_ = 271.0_real32
       charge_ = 106.0_real32
       radius_ = 1.43_real32
    case('Bh')
       mass_ = 270.0_real32
       charge_ = 107.0_real32
       radius_ = 1.41_real32
    case('Hs')
       mass_ = 277.0_real32
       charge_ = 108.0_real32
       radius_ = 1.34_real32
    case('Mt')
       mass_ = 276.0_real32
       charge_ = 109.0_real32
       radius_ = 1.29_real32
    case('Ds')
       mass_ = 281.0_real32
       charge_ = 110.0_real32
       radius_ = 1.28_real32
    case('Rg')
       mass_ = 280.0_real32
       charge_ = 111.0_real32
       radius_ = 1.21_real32
    case('Cn')
       mass_ = 285.0_real32
       charge_ = 112.0_real32
       radius_ = 1.22_real32
    case('Nh')
       mass_ = 284.0_real32
       charge_ = 113.0_real32
       radius_ = 1.21_real32
    case('Fl')
       mass_ = 289.0_real32
       charge_ = 114.0_real32
       radius_ = 1.21_real32
    case('Mc')
       mass_ = 288.0_real32
       charge_ = 115.0_real32
       radius_ = 1.21_real32
    case('Lv')
       mass_ = 293.0_real32
       charge_ = 116.0_real32
       radius_ = 1.21_real32
    case('Ts')
       mass_ = 294.0_real32
       charge_ = 117.0_real32
       radius_ = 1.21_real32
    case('Og')
       mass_ = 294.0_real32
       charge_ = 118.0_real32
       radius_ = 1.21_real32
    case default
       ! handle unknown element
       mass_ = 0.0_real32
       charge_ = 0.0_real32
       radius_ = 0.0_real32
    end select

    !---------------------------------------------------------------------------
    ! Return the values
    !---------------------------------------------------------------------------
    if(present(mass)) mass = mass_
    if(present(charge)) charge = charge_
    if(present(radius)) radius = radius_

  end subroutine get_element_properties
!###############################################################################

end module atomstruc__elements
