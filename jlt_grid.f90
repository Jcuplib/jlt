module jlt_grid
  use jlt_constant, only : STR_SHORT
  use jlt_grid_class, only : grid_class
  
!--------------------------------   public  ----------------------------------!
  public :: init_grid
  public :: def_grid
  public :: end_grid_def
  public :: get_num_of_my_grid
  public :: get_grid_name
  public :: get_grid_ptr

!--------------------------------   private  ---------------------------------!

  integer, parameter :: MAX_GRID = 8

  integer :: num_of_my_grid = 0
  type(grid_class), pointer :: my_grid(:)

  integer :: my_comp_id
  integer :: my_rank
  integer :: my_size
  
contains


!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine init_grid(comp_id, rank, size)
  implicit none
  integer, intent(IN) :: comp_id
  integer, intent(IN) :: rank
  integer, intent(IN) :: size

  my_comp_id = comp_id
  my_rank    = rank
  my_size    = size
  
  allocate(my_grid(MAX_GRID))
  
end subroutine init_grid
  
!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine def_grid(grid_index, comp_name, grid_name)
  implicit none
  integer, intent(IN) :: grid_index(:)
  character(len=*), intent(IN) :: comp_name
  character(len=*), intent(IN) :: grid_name

  num_of_my_grid = num_of_my_grid + 1
  my_grid(num_of_my_grid) = grid_class(my_comp_id, my_rank, my_size)
  call my_grid(num_of_my_grid)%def_grid(grid_index, trim(comp_name), trim(grid_name))

end subroutine def_grid
  
!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine end_grid_def()
  implicit none

end subroutine end_grid_def

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function get_num_of_my_grid()
  implicit none

  get_num_of_my_grid = num_of_my_grid

end function get_num_of_my_grid

!=======+=========+=========+=========+=========+=========+=========+=========+

character(STR_SHORT) function get_grid_name(grid_num)
  implicit none
  integer, intent(IN) :: grid_num

  get_grid_name = my_grid(grid_num)%get_grid_name()

end function get_grid_name

!=======+=========+=========+=========+=========+=========+=========+=========+

function get_grid_ptr(grid_name) result(res)
  use jlt_utils, only : error
  implicit none
  character(len=*), intent(IN) :: grid_name
  type(grid_class), pointer :: res
  integer :: i

  res => null()
  
  do i = 1, num_of_my_grid
     if (trim(grid_name) == trim(my_grid(i)%get_grid_name())) then
        res => my_grid(i)
        return
     end if
  end do
  
  call error("[get_grid_ptr]:no such grid name, grid_name = ", trim(grid_name))

end function get_grid_ptr

!=======+=========+=========+=========+=========+=========+=========+=========+

end module jlt_grid
