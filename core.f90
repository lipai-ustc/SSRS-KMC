module core
   
  !use teio,    only: LINELEN,&
  !                   SPCSLEN

  !use input,   only: InputData

  !use io,      only: io_adjustl

  implicit none
  private 
  save

  public  :: att_core,  &
             det_core

  !private :: 

  integer  :: edge

contains

  !----------------------------------------------------------------!

  subroutine att_core(id,x,y)
    implicit none
    integer,       intent(in)      :: id
    integer,       intent(in)      :: x,y
    integer                        :: i,j,k
  end subroutine att_core
  !----------------------------------------------------------------!

  subroutine det_core(id)
    implicit none
    integer,       intent(in)      :: id
    integer                        :: i,j,k
  end subroutine det_core
 
end module core
