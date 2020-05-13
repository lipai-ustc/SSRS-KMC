module statistics
   
  use teio,    only: LINELEN,&
                     SPCSLEN

  use input,   only: InputData

  use io,      only: io_adjustl

  implicit none
  private 
  save

  public  :: init_statistics,  &
             update_event_num, &
             print_event_count

  !private :: 

  type, public :: SpecI
      integer                        :: SID
      integer                        :: Coor(2)
      type(SpecI), pointer           :: prev
      type(SpecI), pointer           :: next
      double precision, dimension(:),allocatable       :: ProbNB  ! diffu/merge/att in this array
      double precision, dimension(:),allocatable       :: ProbMerge ! Merge with 
      double precision, dimension(:),allocatable       :: ProbDecomp
      double precision               :: ProbDesorp
      double precision               :: ProbSum
  end type SpecI

  type, public :: eventnum
      integer(kind=8), dimension(:), allocatable         :: adsorb !Gas adsorption
      integer(kind=8), dimension(:), allocatable         :: diff
      integer(kind=8), dimension(:), allocatable         :: attach
      integer(kind=8), dimension(:), allocatable         :: detach
      integer(kind=8), dimension(:), allocatable         :: desorb
      integer(kind=8), dimension(:), allocatable         :: decomp
      integer(kind=8), dimension(:), allocatable         :: combine
  end type eventnum

  type(eventnum), public            :: Evt

contains

  !----------------------------------------------------------------!

  subroutine init_statistics(inp)
    implicit none
    type(InputData), intent(in)    :: inp
    integer                        :: i,j,k
    if(inp%GasNum>0) then
        allocate(Evt%adsorb(inp%GasNum))
        Evt%adsorb=0
    end if
    if(inp%SpecNum>0) then
        allocate(Evt%diff(inp%SpecNum))
        Evt%diff=0
        if(inp%FlagCore) then
            allocate(Evt%attach(inp%SpecNum))
            allocate(Evt%detach(inp%SpecNum))
            Evt%attach=0
            Evt%detach=0
        end if
    end if
    allocate(Evt%desorb(inp%DesorpNum))
    Evt%desorb=0
    if(inp%DecompNum>0) then
        allocate(Evt%decomp(inp%DecompNum))
        Evt%decomp=0
    end if
    if(inp%MergeNum>0) then
        allocate(Evt%combine(inp%MergeNum))
        Evt%combine=0
    end if
  end subroutine init_statistics

  !----------------------------------------------------------------!
  
  subroutine update_event_num(event,id)
    character(*),    intent(in)  :: event
    integer,         intent(in)  :: id
    select case(trim(adjustl(event)))
      case("adsorption")
        Evt%adsorb(id)=Evt%adsorb(id)+1
      case("diffusion")
        Evt%diff(id)=Evt%diff(id)+1
      case("attachment")
        Evt%attach(id)=Evt%attach(id)+1
      case("detachment")
        Evt%detach(id)=Evt%detach(id)+1
      case("merge")
        Evt%combine(id)=Evt%combine(id)+1
      case("decomposition")
        Evt%decomp(id)=Evt%decomp(id)+1
      case("desorption")
        Evt%desorb(id)=Evt%desorb(id)+1
    end select
  end subroutine update_event_num

  !----------------------------------------------------------------!
  
  subroutine print_event_count(u)
    integer,         intent(in)  :: u
    integer                      :: i
    if(allocated(Evt%adsorb)) then    ! adsorb
        write(u,'(A$)') "adsorption: "
        do i=1,size(Evt%adsorb)
            write(u,'(A," "$)') trim(io_adjustl(Evt%adsorb(i)))
        end do
        write(u,*)
    end if
    if(allocated(Evt%diff)) then      ! diffusion
        write(u,'(A$)') "diffusion: "
        do i=1,size(Evt%diff)
            write(u,'(A," "$)') trim(io_adjustl(Evt%diff(i)))
        end do
        write(u,*)
    end if
    if(allocated(Evt%attach)) then    ! attachment
        write(u,'(A$)') "attachment: "
        do i=1,size(Evt%attach)
            write(u,'(A," "$)') trim(io_adjustl(Evt%attach(i)))
        end do
        write(u,*)
    end if
    if(allocated(Evt%detach)) then    ! detachment
        write(u,'(A$)') "detachment: "
        do i=1,size(Evt%detach)
            write(u,'(A," "$)') trim(io_adjustl(Evt%detach(i)))
        end do
        write(u,*)
    end if
    if(allocated(Evt%combine)) then   ! merge
        write(u,'(A$)') "merge: "
        do i=1,size(Evt%combine)
            write(u,'(A," "$)') trim(io_adjustl(Evt%combine(i)))
        end do
        write(u,*)
    end if
    if(allocated(Evt%decomp)) then    ! decomposition
        write(u,'(A$)') "decomposition: "    
        do i=1,size(Evt%decomp)
            write(u,'(A," "$)') trim(io_adjustl(Evt%decomp(i)))
        end do
        write(u,*)
    end if
    if(allocated(Evt%desorb)) then    ! desorption
        write(u,'(A$)') "desorption: "    
        do i=1,size(Evt%desorb)
            write(u,'(A," "$)') trim(io_adjustl(Evt%desorb(i)))
        end do
        write(u,*)
    end if

  end subroutine print_event_count

end module statistics
