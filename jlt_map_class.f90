module jlt_map_class
  use jlt_constant, only : STR_SHORT, STR_MID, STR_LONG
  use jlt_grid_class, only : grid_class
  
!--------------------------------   public  ----------------------------------!

public :: map_class

!--------------------------------   private  ---------------------------------!


type target_map_info
     ! local exchange table
     integer                       :: num_of_exchange_rank ! number of exchange target rank
     integer, pointer              :: exchange_rank(:)     ! target rank 
     integer, pointer              :: num_of_exchange(:)   ! number of exchange index of each target rank
     integer, pointer              :: offset(:)            ! offset address 
     integer                       :: exchange_data_size   ! 
     integer, pointer              :: exchange_index(:)    ! local index array for data exchange
     integer, pointer              :: conv_table(:)        ! conersion table from grid array to exchange array (valid only send side)
end type target_map_info

type map_class
   private
   type(grid_class), pointer :: my_grid
   integer              :: comp_id                  ! my component id
   logical              :: send_recv_flag = .true.  ! send map or recv map (.true. = send map)
   integer              :: num_of_exchange_index    ! local size of exchange data array
   integer, pointer     :: exchange_index(:)        ! global index of my exchange grid points
   integer, allocatable :: conv_table(:)            ! array conversion table from exchange index to data grid index

   type(target_map_info), pointer :: tmap_info      ! exchange info of my target map

 contains
   procedure :: set_map_grid                        ! subroutine (my_grid_ptr)
   procedure :: get_comp_id                         ! integer function ()
   procedure :: get_my_rank                         ! integer function ()
   procedure :: set_send_flag                       ! subroutine ()         turn on send flag
   procedure :: set_recv_flag                       ! subroutine ()         turn on recv_flag
   procedure :: is_send_map                         ! logical function ()
   procedure :: is_recv_map                         ! logical function ()    
   procedure :: get_exchange_array_size             ! integer function ()
   procedure :: set_exchange_index                  ! subroutine (index_ptr)
   procedure :: set_target_map_info                 ! subroutine (target_rank)
   procedure :: local_2_exchange                    ! subroutine (local_data, exchange_data)
   procedure :: exchange_2_local                    ! subroutine (exchange_data, local_data)
   procedure :: get_num_of_target_pe                ! integer function ()
   procedure :: get_target_rank                     ! integer function (pe_num)
   procedure :: get_offset                          ! integer function (pe_num)
   procedure :: get_num_of_data                     ! integer function (pe_num)
   
end type map_class

interface map_class
   module procedure init_map_class
end interface map_class

contains

!=======+=========+=========+=========+=========+=========+=========+=========+

type(map_class) function init_map_class(comp_id)
  implicit none
  type(map_class) :: my_map_class
  integer, intent(IN) :: comp_id

  my_map_class%comp_id = comp_id
  my_map_class%num_of_exchange_index = 0
  my_map_class%my_grid    => null()
  my_map_class%tmap_info  => null()
  init_map_class = my_map_class
 
end function init_map_class

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine set_map_grid(self, my_grid)
  implicit none
  class(map_class) :: self
  type(grid_class), pointer :: my_grid

  self%my_grid => my_grid

end subroutine set_map_grid

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine set_send_flag(self)
  implicit none
  class(map_class) :: self

  self%send_recv_flag = .true.

end subroutine set_send_flag

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine set_recv_flag(self)
  implicit none
  class(map_class) :: self

  self%send_recv_flag = .false.

end subroutine set_recv_flag

!=======+=========+=========+=========+=========+=========+=========+=========+

logical function is_send_map(self)
  implicit none
  class(map_class) :: self

  is_send_map = self%send_recv_flag

end function is_send_map

!=======+=========+=========+=========+=========+=========+=========+=========+

logical function is_recv_map(self)
  implicit none
  class(map_class) :: self

  is_recv_map = .not.self%send_recv_flag

end function is_recv_map

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function get_comp_id(self) 
  implicit none
  class(map_class) :: self

  get_comp_id = self%comp_id

end function get_comp_id

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function get_my_rank(self) 
  implicit none
  class(map_class) :: self

  get_my_rank = self%my_grid%get_my_rank()

end function get_my_rank

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function get_exchange_array_size(self)
  implicit none
  class(map_class) :: self

  get_exchange_array_size = self%num_of_exchange_index

end function get_exchange_array_size

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine set_exchange_index(self, index_ptr)
  use jlt_utils, only : sort_int_1d, binary_search
  implicit none
  class(map_class) :: self
  integer, pointer :: index_ptr(:)
  integer, pointer :: grid_index(:)
  integer, pointer :: sorted_index(:)
  integer, pointer :: sorted_pos(:)
  integer :: grid_size
  integer :: res
  integer :: i
  
  self%num_of_exchange_index = size(index_ptr)
  
  self%exchange_index => index_ptr

  allocate(self%conv_table(self%num_of_exchange_index))

  grid_index => self%my_grid%get_grid_index_ptr()
  grid_size = size(grid_index)
  
  allocate(sorted_index(grid_size))
  allocate(sorted_pos(grid_size))

  do i = 1, grid_size
     sorted_index(i) = grid_index(i)
     sorted_pos(i) = i
  end do

  call sort_int_1d(grid_size, sorted_index, sorted_pos)

  do i = 1, self%num_of_exchange_index
     res = binary_search(sorted_index, self%exchange_index(i))
     self%conv_table(i) = sorted_pos(res)
  end do

end subroutine set_exchange_index

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine set_target_map_info(self, target_rank)
  implicit none
  class(map_class) :: self
  integer, pointer :: target_rank(:)
  integer :: i, counter

  allocate(self%tmap_info)
  
  self%tmap_info%target_rank => target_rank
  self%tmap_info%my_data_size = size(target_rank)

  self%tmap_info%num_of_pe = 1
  do i = 2, size(target_rank)
     if (target_rank(i) /= target_rank(i-1)) then
        self%tmap_info%num_of_pe = self%tmap_info%num_of_pe + 1
     end if
  end do
  
  allocate(self%tmap_info%num_of_data(self%tmap_info%num_of_pe))
  allocate(self%tmap_info%offset(self%tmap_info%num_of_pe))

  counter = 1
  self%tmap_info%num_of_data(:) =1

  do i = 2, size(target_rank)
     if (target_rank(i) == target_rank(i-1)) then
        self%tmap_info%num_of_data(counter) = self%tmap_info%num_of_data(counter) + 1
     else
        counter = counter + 1
     end if
  end do

  self%tmap_info%offset(1) = 0
  do i = 2, self%tmap_info%num_of_pe
     self%tmap_info%offset(i) = self%tmap_info%offset(i-1) + self%tmap_info%num_of_data(i-1)
  end do
  
end subroutine set_target_map_info


!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine local_2_exchange(self, local_data, exchange_data)
  implicit none
  class(map_class)       :: self
  real(kind=8), intent(IN)    :: local_data(:)
  real(kind=8), intent(INOUT) :: exchange_data(:)
  integer :: i

  do i = 1, self%num_of_exchange_index
     exchange_data(i) = local_data(self%conv_table(i))
  end do

end subroutine local_2_exchange

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine exchange_2_local(self, local_data, exchange_data)
  implicit none
  class(map_class)            :: self
  real(kind=8), intent(INOUT) :: local_data(:)
  real(kind=8), intent(IN)    :: exchange_data(:)
  integer :: i

  do i = 1, self%num_of_exchange_index
     local_data(self%conv_table(i)) = exchange_data(i)
  end do

end subroutine exchange_2_local

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function get_num_of_target_pe(self)
  implicit none
  class(map_class) :: self

  get_num_of_target_pe = self%tmap_info%num_of_pe

end function get_num_of_target_pe

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function get_target_rank(self, pe_num)
  implicit none
  class(map_class)     :: self
  integer, intent(IN) :: pe_num

  get_target_rank = self%tmap_info%target_rank(pe_num)

end function get_target_rank

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function get_offset(self, pe_num)
  implicit none
  class(map_class)     :: self
  integer, intent(IN) :: pe_num

  get_offset = self%tmap_info%offset(pe_num)

end function get_offset

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function get_num_of_data(self, pe_num)
  implicit none
  class(map_class)     :: self
  integer, intent(IN) :: pe_num

  get_num_of_data = self%tmap_info%num_of_data(pe_num)

end function get_num_of_data

!=======+=========+=========+=========+=========+=========+=========+=========+

end module jlt_map_class
