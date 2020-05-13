module mesh

  use input,   only: InputData,&
                     print_err

  use io,      only: io_adjustl

  use statistics,    only: SpecI

  implicit none
  private
  save

  public :: get_neighb,&
            get_nb_direct,&
            init_mesh,&
            save_grid

!  private ::

  !----------------------------------------!
  !          grid data type                !
  !----------------------------------------!

  type, public  :: MeshGrid
    integer                 :: occ
    type(SpecI), pointer    :: spec
  end type MeshGrid
  type(MeshGrid), dimension(:,:), allocatable,public  :: Grid

  type, public :: CoorType
    integer                                    :: NBNum
    ! NBNum for calculate neighboring sites
    integer, dimension(:,:,:),allocatable      :: NBSite
    integer, dimension(:,:),  allocatable      :: NB
    integer, dimension(2)                      :: Drct
  end type CoorType

  type(CoorType),                             public  :: Coor
  integer,                                    public  :: MeshType
  integer, dimension(2),                      public  :: MeshSize
  integer, dimension(2),                      public  :: CoreSize
  integer,                                    public  :: SiteNum  ! Num of all available sites

  !---------------------------------------------------!

contains
  
  !---------------------------------------------------!
  
  subroutine init_mesh(inp)
    implicit none
    type(InputData), intent(in)   :: inp
    integer                       :: i,j

    ! update meshgrid
    MeshType=inp%MeshType
    MeshSize=inp%MeshSize
    allocate(Grid(MeshSize(1),MeshSize(2)))
    do i=1,MeshSize(1)
      do j=1,MeshSize(2)
        Grid(i,j)%occ=0    !no occupacy
        Grid(i,j)%spec=>null()    !no occupacy
      end do
    end do

    if(inp%FlagCore) then
        CoreSize=inp%CoreSize
        ! Build Core
        ! some points that cannot be reached in C++; check here!
        do i=MeshSize(1)/2-CoreSize(1)/2,MeshSize(1)/2+CoreSize(1)/2
          do j=MeshSize(2)/2-CoreSize(2)/2,MeshSize(2)/2+CoreSize(2)/2
            Grid(i,j)%occ=-1    !occupied by core
          end do
        end do
    else
        CoreSize=0
    end if

    SiteNum=MeshSize(1)*MeshSize(2)-CoreSize(1)*CoreSize(2)

    if(MeshType==1) then 
      Coor%NBNum=3    ! three neighoring sites for tri mesh
      allocate(Coor%NBSite(2,3,2),Coor%NB(3,2))
      ! NBSite(x even/odd; 3 directions; x/y)
      Coor%NBSite(1,1,1)=-1
      Coor%NBSite(1,1,2)=1
      Coor%NBSite(1,2,1)=-1
      Coor%NBSite(1,2,2)=0
      Coor%NBSite(1,3,1)=1
      Coor%NBSite(1,3,2)=0
      ! x odd
      Coor%NBSite(2,1,1)=-1
      Coor%NBSite(2,1,2)=0
      Coor%NBSite(2,2,1)=1
      Coor%NBSite(2,2,2)=0
      Coor%NBSite(2,3,1)=1
      Coor%NBSite(2,3,2)=-1
    else 
      call print_err("MeshType /=1 is not ready now")
    end if

  end subroutine init_mesh

  !---------------------------------------------------!

  subroutine get_neighb(x,y)
    implicit none
    integer, intent(in)                :: x,y
    integer                            :: i
    if(MeshType==1) then
      if(mod(x,2)==0) then  ! for x even
        do i=1, Coor%NBNum
          Coor%NB(i,1)=x+Coor%NBSite(1,i,1) 
          Coor%NB(i,2)=y+Coor%NBSite(1,i,2)
        end do
      else                       ! for x odd
        do i=1, Coor%NBNum
          Coor%NB(i,1)=x+Coor%NBSite(2,i,1)
          Coor%NB(i,2)=y+Coor%NBSite(2,i,2)
        end do
      endif
      do i=1, Coor%NBNum
        if(Coor%NB(i,1)>MeshSize(1)) Coor%NB(i,1)=Coor%NB(i,1)-MeshSize(1) ! for x
        if(Coor%NB(i,1)<1)           Coor%NB(i,1)=Coor%NB(i,1)+MeshSize(1)
        if(Coor%NB(i,2)>MeshSize(2)) Coor%NB(i,2)=Coor%NB(i,2)-MeshSize(2) ! for y
        if(Coor%NB(i,2)<1)           Coor%NB(i,2)=Coor%NB(i,2)+MeshSize(2)
      end do
    else 
      call print_err("MeshType /=1 is not ready now")
    endif
  end subroutine get_neighb

  !---------------------------------------------------!

  subroutine get_nb_direct(x,y,drct)
    implicit none
    integer, intent(in)               :: x,y
    integer, intent(in)               :: drct
    if(MeshType==1) then
      if(mod(x,2)==0) then  ! for x even
        Coor%Drct(1)=x+Coor%NBSite(1,drct,1)
        Coor%Drct(2)=y+Coor%NBSite(1,drct,2)
      else
        Coor%Drct(1)=x+Coor%NBSite(2,drct,1)
        Coor%Drct(2)=y+Coor%NBSite(2,drct,2)
      end if
    else 
      call print_err("MeshType /=1 is not ready now")
    end if

    !periodic condition
    if(Coor%Drct(1)>MeshSize(1))  Coor%Drct(1)=Coor%Drct(1)-MeshSize(1)
    if(Coor%Drct(1)<1)            Coor%Drct(1)=Coor%Drct(1)+MeshSize(1)
    if(Coor%Drct(2)>MeshSize(2))  Coor%Drct(2)=Coor%Drct(2)-MeshSize(2)
    if(Coor%Drct(2)<1)            Coor%Drct(2)=Coor%Drct(2)+MeshSize(2)
  end subroutine

  !---------------------------------------------------!

  subroutine save_grid(u)
    implicit none
    integer, intent(in)               :: u
    integer                           :: i,j
    integer                           :: n

    n=0
    write(u,'(A)') "GridOcc"
    do i=1, MeshSize(1)
        do j=1, MeshSize(2)
            if(Grid(i,j)%occ/=0.and.Grid(i,j)%occ/=-1) then
                n=n+1
                write(u,'(5A)') trim(io_adjustl(i)),' ',trim(io_adjustl(j)),&
                           ' ',trim(io_adjustl(Grid(i,j)%occ))
            end if
        end do
    end do
    write(u,'(2A)') "GridOccNum ",trim(io_adjustl(n))

  end subroutine

end module
