module jlt_grid_class
  use jlt_constant, only : STR_SHORT
  
!--------------------------------   public  ----------------------------------!

public :: grid_class

!--------------------------------   private  ---------------------------------!

type grid_class
   private
   integer :: comp_id
   integer :: my_rank
   integer :: my_size
   character(len=STR_SHORT) :: comp_name
   character(len=STR_SHORT) :: grid_name
   integer :: num_of_grid_index
   integer, pointer         :: grid_index(:)        ! global index of my grid points
 contains
   procedure :: def_grid
   procedure :: get_grid_name
   procedure :: get_comp_id
   procedure :: get_my_rank
   procedure :: get_my_size
   procedure :: get_grid_size
   procedure :: get_grid_index_ptr
end type grid_class

interface grid_class
   module procedure init_grid_class
end interface grid_class


contains

!=======+=========+=========+=========+=========+=========+=========+=========+

type(grid_class) function init_grid_class(comp_id, my_rank, my_size)
  implicit none
  integer, intent(IN)       :: comp_id
  integer, intent(IN)       :: my_rank
  integer, intent(IN)       :: my_size
  type(grid_class) :: my_grid

  my_grid%comp_id = comp_id
  my_grid%my_rank = my_rank
  my_grid%my_size = my_size
  my_grid%num_of_grid_index = 0

  init_grid_class = my_grid

end function init_grid_class


!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine def_grid(self, grid_index, comp_name, grid_name)
  implicit none
  class(grid_class)            :: self
  integer, intent(IN)          :: grid_index(:)
  character(len=*), intent(IN) :: comp_name
  character(len=*), intent(IN) :: grid_name

  self%comp_name = trim(comp_name)
  self%grid_name = trim(grid_name)

  self%num_of_grid_index = size(grid_index)
  allocate(self%grid_index(self%num_of_grid_index))
  self%grid_index(:) = grid_index(:)
  
end subroutine def_grid
  
!=======+=========+=========+=========+=========+=========+=========+=========+

character(len=STR_SHORT) function get_grid_name(self)
  implicit none
  class(grid_class) :: self

  get_grid_name = trim(self%grid_name)

end function get_grid_name

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function get_comp_id(self)
  implicit none
  class(grid_class) :: self

  get_comp_id = self%comp_id

end function get_comp_id

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function get_my_rank(self)
  implicit none
  class(grid_class) :: self

  get_my_rank = self%my_rank

end function get_my_rank

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function get_my_size(self)
  implicit none
  class(grid_class) :: self

  get_my_size = self%my_size

end function get_my_size

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function get_grid_size(self)
  implicit none
  class(grid_class)  :: self

  get_grid_size = self%num_of_grid_index

end function get_grid_size

!=======+=========+=========+=========+=========+=========+=========+=========+

function get_grid_index_ptr(self) result(res)
  implicit none
  class(grid_class)  :: self
  integer, pointer :: res(:)
  
  res => self%grid_index

end function get_grid_index_ptr

!=======+=========+=========+=========+=========+=========+=========+=========+

end module jlt_grid_class
