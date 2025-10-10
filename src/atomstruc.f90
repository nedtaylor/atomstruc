module atomstruc
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, atomstruc!"
  end subroutine say_hello
end module atomstruc
