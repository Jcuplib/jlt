module jlt_exchange_class
  use jlt_constant, only   : STR_SHORT, STR_MID, STR_LONG
  use jlt_grid_class, only : grid_class
  implicit none
  private

!--------------------------------   public   ---------------------------------!

  public :: buffer_class            ! exchange buffer class
  public :: exchange_class

!--------------------------------   private  ---------------------------------!
  type buffer_class ! data buffer for MPI send/recv
     integer :: target_rank
     integer :: num_of_data
     integer :: num_of_layer
     real(kind=8), pointer :: buffer(:, :)
  end type buffer_class

  type exchange_map_info
     ! local exchange table
     integer                       :: num_of_exchange_rank = 0 ! number of exchange target rank
     integer, pointer              :: exchange_rank(:)     ! target rank 
     integer, pointer              :: num_of_exchange(:)   ! number of exchange index of each target rank
     integer, pointer              :: offset(:)            ! offset address 
     integer                       :: exchange_data_size   ! 
     integer, pointer              :: exchange_index(:)    ! local index array for data exchange
     integer, pointer              :: conv_table(:)        ! conersion table from grid array to exchange array (valid only send side)
  end type exchange_map_info

  type send_map_info
     integer                       :: intpled_data_size    ! data size of interpolation
     integer, pointer              :: intpled_index(:)     ! interpolation data index
     integer, pointer              :: conv_table(:)        ! conversion table from grid array to interpolation data array
  end type send_map_info
  
  type recv_map_info
     integer                       :: intpled_data_size    ! interpolated data size
     integer, pointer              :: intpled_index(:)     ! interpolation data index
     integer, pointer              :: conv_table(:)        ! conversion table from interpolated data array to grid array 
  end type recv_map_info
  
  type exchange_class
     private
     character(len=STR_SHORT)      :: my_name
     character(len=STR_SHORT)      :: send_comp_name
     character(len=STR_SHORT)      :: send_grid_name
     character(len=STR_SHORT)      :: recv_comp_name
     character(len=STR_SHORT)      :: recv_grid_name
     integer                       :: map_tag
     integer                       :: send_comp_id
     integer                       :: recv_comp_id
     ! local remapping table
     integer                       :: index_size           ! local_index_size
     integer, pointer              :: send_grid_index(:)   ! local index 
     integer, pointer              :: recv_grid_index(:)   ! local index
     real(kind=8), pointer         :: coef(:)              ! local coef
     integer, pointer              :: target_rank(:)       ! local target rank
     ! array conversion table for interpolation (valid only recv side)
     integer, pointer              :: send_conv_table(:)   ! conv. table from interpolation index to send grid index
     integer, pointer              :: recv_conv_table(:)   ! conv. table from interpolation index to recv grid index

     type(exchange_map_info)       :: ex_map
     type(send_map_info)           :: send_map             ! varid for send side interpolation
     type(recv_map_info)           :: my_map
     
     logical                       :: intpl_flag = .false.  ! my interpolation or not

     
   contains
     procedure :: set_mapping_table         ! subroutine (send_comp, send_grid, recv_comp, recv_grid, map_tag, send_index, recv_index, coef)
     procedure :: is_my_exchange            ! logicail function (send_comp, send_grid, recv_comp, recv_grid, map_tag)
     procedure :: get_my_name               ! character(len=STR_SHORT) function()
     procedure :: get_exchange_data_size    ! integer function ()
     procedure :: get_exchange_buffer_size  ! integer function ()
     procedure :: get_num_of_target_rank    ! integer function ()
     procedure :: get_target_rank           ! integer function (target_rank_num)
     procedure :: get_target_array_size     ! integer function (target_rank_num)
     procedure :: is_my_intpl               ! logical function ()
     procedure :: is_send_intpl             ! logical function ()
     procedure :: local_2_exchange          ! subroutine (grid_data, exchange_data)
     procedure :: exchange_2_local          ! subroutine (exchange_data, grid_data)
     procedure :: target_2_exchange_buffer  ! subroutine (target_buffer, exchange_buffer)
     procedure :: send_data_1d              ! subroutine (data, exchange_tag)
     procedure :: recv_data_1d              ! subroutine (data, exchange_tag)
     procedure :: send_data_2d              ! subroutine (data, exchange_tag)
     procedure :: recv_data_2d              ! subroutine (data, exchange_tag)
     procedure :: buffer_2_recv_data        ! subroutine (exchange_buffer, recv_data) 
     procedure :: interpolate_data          ! subroutine (send_data, recv_data, num_of_layer)
  end type exchange_class

  interface exchange_class
     module procedure init_exchange_class
  end interface exchange_class
  
contains

!=======+=========+=========+=========+=========+=========+=========+=========+

type(exchange_class) function init_exchange_class(my_name, intpl_flag)
  implicit none
  character(len=*), intent(IN) :: my_name
  logical, intent(IN)          :: intpl_flag
  type(exchange_class) :: my_exchange_class

  my_exchange_class%my_name    = trim(my_name)
  my_exchange_class%intpl_flag = intpl_flag  

  init_exchange_class = my_exchange_class

end function init_exchange_class
  
!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine set_mapping_table(self, send_comp_name, send_grid_name, recv_comp_name, recv_grid_name, &
                             map_tag, intpl_flag, send_grid, recv_grid, coef)
  use jlt_grid, only       : get_grid_ptr
  use jlt_grid_class, only : grid_class
  use jlt_comp, only       : get_comp_id_from_name
  use jlt_mpi_lib, only    : jml_BcastLocal, jml_GetComm
  use jlt_utils, only      : put_log
  use jlt_remapping, only  : make_local_mapping_table, get_target_grid_rank, reorder_index_by_target_rank, &
                             make_exchange_table, make_conversion_table
  implicit none
  class(exchange_class)        :: self
  character(len=*), intent(IN) :: send_comp_name, send_grid_name
  character(len=*), intent(IN) :: recv_comp_name, recv_grid_name
  integer, intent(IN)          :: map_tag
  logical, intent(IN)          :: intpl_flag
  integer, intent(IN)          :: send_grid(:), recv_grid(:)
  real(kind=8), intent(IN)     :: coef(:)
  
  self%send_comp_name     = trim(send_comp_name)
  self%send_grid_name     = trim(send_grid_name)
  self%recv_comp_name     = trim(recv_comp_name)
  self%recv_grid_name     = trim(recv_grid_name)
  self%map_tag            = map_tag
  self%send_comp_id       = get_comp_id_from_name(trim(send_comp_name))
  self%recv_comp_id       = get_comp_id_from_name(trim(recv_comp_name))
  self%intpl_flag         = intpl_flag
  
  if (self%is_send_intpl()) then
     call set_mapping_table_send_intpl(self, send_comp_name, send_grid_name, recv_comp_name, recv_grid_name, &
          map_tag, send_grid, recv_grid, coef)
  else
     call set_mapping_table_recv_intpl(self, send_comp_name, send_grid_name, recv_comp_name, recv_grid_name, &
          map_tag, send_grid, recv_grid, coef)
  end if
  
end subroutine set_mapping_table

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine set_mapping_table_send_intpl(self, send_comp_name, send_grid_name, recv_comp_name, recv_grid_name, &
                             map_tag, send_grid, recv_grid, coef)
  use jlt_grid, only       : get_grid_ptr
  use jlt_grid_class, only : grid_class
  use jlt_comp, only       : get_comp_id_from_name
  use jlt_mpi_lib, only    : jml_BcastLocal, jml_GetComm
  use jlt_utils, only      : put_log
  use jlt_remapping, only  : make_local_mapping_table, make_local_mapping_table_no_sort, get_target_grid_rank, reorder_index_by_target_rank, &
                             make_exchange_table, make_conversion_table
  use mpi
  implicit none
  class(exchange_class)        :: self
  character(len=*), intent(IN) :: send_comp_name, send_grid_name
  character(len=*), intent(IN) :: recv_comp_name, recv_grid_name
  integer, intent(IN)          :: map_tag
  integer, intent(IN)          :: send_grid(:), recv_grid(:)
  real(kind=8), intent(IN)     :: coef(:)
  integer                      :: dest_comp_id, source_comp_id
  integer, allocatable         :: index_array(:)
  integer                      :: int_array(1)
  character(len=STR_LONG)      :: log_str
  type(grid_class), pointer    :: my_grid
  integer :: i, j
  integer, pointer             :: my_grid_index(:)
  integer                      :: num_of_exchange, offset
  integer :: my_rank, ierr

  call mpi_comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  
  call put_log("------------------------------------------------------------------------------------------")
  call put_log("jlt_exchange_class : set_mapping_table start")
  
  if (trim(self%my_name) == trim(send_comp_name)) then
     source_comp_id = self%send_comp_id
     dest_comp_id   = self%recv_comp_id
  else
     source_comp_id = self%recv_comp_id
     dest_comp_id   = self%send_comp_id
  end if

  call put_log("jlt_exchange_class : set_mapping_table, make_local_mapping_table  start")

  
  if (trim(self%my_name) == trim(send_comp_name)) then
     my_grid => get_grid_ptr(trim(send_grid_name))
     
     call make_local_mapping_table_no_sort(send_grid, recv_grid, coef, my_grid%get_grid_index_ptr(), &
                                   self%index_size, self%send_grid_index, self%recv_grid_index, self%coef)
   
     allocate(self%target_rank(self%index_size))
     call get_target_grid_rank(recv_comp_name, recv_grid_name, self%recv_grid_index, self%target_rank)
  else
     my_grid => get_grid_ptr(trim(recv_grid_name))
     call make_local_mapping_table_no_sort(recv_grid, send_grid, coef, my_grid%get_grid_index_ptr(), &
                                   self%index_size, self%recv_grid_index, self%send_grid_index, self%coef)
     allocate(self%target_rank(self%index_size))
     call get_target_grid_rank(send_comp_name, send_grid_name, self%send_grid_index, self%target_rank)
  end if
  
  call put_log("jlt_exchange_class : set_mapping_table, make_local_mapping_table  end")


  !write(399+source_comp_id*100 + my_grid%get_my_rank(), *) "global mapping table"
  !write(399+source_comp_id*100 + my_grid%get_my_rank(), *) "send_comp = "//trim(send_comp_name)//", recv_comp = "//trim(recv_comp_name)

  do i = 1, size(send_grid)
     !write(399+source_comp_id*100 + my_grid%get_my_rank(), *) send_grid(i), recv_grid(i), coef(i)
  end do

  my_grid_index => my_grid%get_grid_index_ptr()

  !write(399+source_comp_id*100 + my_grid%get_my_rank(), *) "my grid_index"
  do i = 1, size(my_grid_index)
     !write(399+source_comp_id*100 + my_grid%get_my_rank(), *) my_grid_index(i)
  end do

  call put_log("jlt_exchange_class : set_mapping_table, reorder_index_by_target_rank start")

  call reorder_index_by_target_rank(self%send_grid_index, self%recv_grid_index, self%coef, self%target_rank)

  !write(399 + my_rank, *) "send_mapping_table_send_intpl, reorder_table ", send_comp_name, recv_comp_name
  do i = 1, size(self%target_rank)
     !write(399 + my_rank, *) self%send_grid_index(i), self%recv_grid_index(i), self%coef(i), self%target_rank(i)
  end do
  
  
  call put_log("jlt_exchange_class : set_mapping_table, reorder_index_by_target_rank end")


  !write(399+my_rank, *) "local mapping table,", self%index_size
  do i = 1, self%index_size
     !write(399+my_rank, *) self%send_grid_index(i), self%recv_grid_index(i), self%coef(i), self%target_rank(i)
  end do

  

  call put_log("jlt_exchange_class : set_mapping_table, make_exchange_table start")

  if (trim(self%my_name) == trim(send_comp_name)) then
     call make_exchange_table(self%target_rank, self%recv_grid_index, &
                              self%ex_map%num_of_exchange_rank, self%ex_map%exchange_rank, &
                              self%ex_map%num_of_exchange, self%ex_map%offset, self%ex_map%exchange_index) 
     self%ex_map%exchange_data_size = size(self%ex_map%exchange_index)
  else
     call make_exchange_table(self%target_rank, self%recv_grid_index, &
                              self%ex_map%num_of_exchange_rank, self%ex_map%exchange_rank, &
                              self%ex_map%num_of_exchange, self%ex_map%offset, self%ex_map%exchange_index) 
     self%ex_map%exchange_data_size = size(self%ex_map%exchange_index)
  end if

  !write(399+ my_rank, *) "exchange table,", send_comp_name, recv_comp_name, self%ex_map%num_of_exchange_rank
  do i = 1, self%ex_map%num_of_exchange_rank
     !write(399+my_rank, *) self%ex_map%exchange_rank(i), &
          !self%ex_map%num_of_exchange(i), self%ex_map%offset(i)
     do j = 1, self%ex_map%num_of_exchange(i)
        !write(399+my_rank, *) self%ex_map%exchange_index(self%ex_map%offset(i)+j)
     end do
  end do


  call put_log("jlt_exchange_class : set_mapping_table, make_exchange_table end")

  call put_log("jlt_exchange_class : set_mapping_table, make_conversion_table start")

  
  if (trim(self%my_name) ==  trim(send_comp_name)) then
    call put_log("jlt_exchange_class : set_mapping_table, make_conversion_table 1")
    call make_send_map_info(self%send_grid_index, self%send_map)
    call make_grid_conv_table(my_grid%get_grid_index_ptr(), self%send_map%intpled_index, self%ex_map%conv_table)
    call make_grid_conv_table(my_grid%get_grid_index_ptr(), self%send_map%intpled_index, self%send_map%conv_table)
    call put_log("jlt_exchange_class : set_mapping_table, make_conversion_table 2")
    call make_conversion_table(self%send_grid_index, self%send_map%intpled_index, self%send_conv_table)
    call make_conversion_table(self%recv_grid_index, self%ex_map%exchange_index, self%recv_conv_table)
    !call make_conversion_table(self%recv_grid_index, self%my_map%intpled_index, self%recv_conv_table)
    !call recv_recv_grid_index(self)
    !call put_log("jlt_exchange_class : set_mapping_table, make_conversion_table 3")
  else
    call put_log("jlt_exchange_class : set_mapping_table, make_conversion_table 4")
    call make_recv_map_info(self%ex_map%exchange_index, self%my_map)
    call put_log("jlt_exchange_class : set_mapping_table, make_conversion_table 5")
    call make_grid_conv_table(self%my_map%intpled_index, self%ex_map%exchange_index, self%recv_conv_table)
    call put_log("jlt_exchange_class : set_mapping_table, make_conversion_table 6")
    !call recv_recv_grid_index(self)
    call put_log("jlt_exchange_class : set_mapping_table, make_conversion_table 7")
    call make_grid_conv_table(my_grid%get_grid_index_ptr(), self%my_map%intpled_index, self%my_map%conv_table)
    call put_log("jlt_exchange_class : set_mapping_table, make_conversion_table 8")
    !call send_recv_grid_index(self, self%recv_grid_index)
  end if

  call put_log("jlt_exchange_class : set_mapping_table, make_conversion_table end")

  
  !write(399+my_rank, *) "exchange rank"
  do i = 1, self%ex_map%num_of_exchange_rank
     !write(399+my_rank, *) self%ex_map%exchange_rank(i), self%ex_map%num_of_exchange(i)
  end do

  !write(399+my_rank, *) "exchange index"
  do i = 1, size(self%ex_map%exchange_index)
     !write(399+my_rank, *) self%ex_map%exchange_index(i)
  end do

  return

  call put_log(" set_mapping_table : exchange map info start")


  call put_log(" set_mapping_table : exchange map info end")

  call put_log("jlt_exchange_class : set_mapping_table end")

  call put_log("------------------------------------------------------------------------------------------")

end subroutine set_mapping_table_send_intpl

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine set_mapping_table_recv_intpl(self, send_comp_name, send_grid_name, recv_comp_name, recv_grid_name, &
                             map_tag, send_grid, recv_grid, coef)
  use jlt_grid, only       : get_grid_ptr
  use jlt_grid_class, only : grid_class
  use jlt_comp, only       : get_comp_id_from_name
  use jlt_mpi_lib, only    : jml_BcastLocal, jml_GetComm
  use jlt_utils, only      : put_log
  use jlt_remapping, only  : make_local_mapping_table_no_sort, get_target_grid_rank, reorder_index_by_target_rank, &
                             make_exchange_table, make_conversion_table
  implicit none
  class(exchange_class)        :: self
  character(len=*), intent(IN) :: send_comp_name, send_grid_name
  character(len=*), intent(IN) :: recv_comp_name, recv_grid_name
  integer, intent(IN)          :: map_tag
  integer, intent(IN)          :: send_grid(:), recv_grid(:)
  real(kind=8), intent(IN)     :: coef(:)
  integer                      :: dest_comp_id, source_comp_id
  integer, allocatable         :: index_array(:)
  integer                      :: int_array(1)
  character(len=STR_LONG)      :: log_str
  type(grid_class), pointer    :: my_grid
  integer :: i
  integer, pointer             :: my_grid_index(:)
  
  call put_log("------------------------------------------------------------------------------------------")
  call put_log("jlt_exchange_class : set_mapping_table start")
  
  if (trim(self%my_name) == trim(send_comp_name)) then
     source_comp_id = self%send_comp_id
     dest_comp_id   = self%recv_comp_id
  else
     source_comp_id = self%recv_comp_id
     dest_comp_id   = self%send_comp_id
  end if

  call put_log("jlt_exchange_class : set_mapping_table, make_local_mapping_table  start")

  
  if (trim(self%my_name) == trim(send_comp_name)) then
     my_grid => get_grid_ptr(trim(send_grid_name))
     
     call make_local_mapping_table_no_sort(send_grid, recv_grid, coef, my_grid%get_grid_index_ptr(), &
                                   self%index_size, self%send_grid_index, self%recv_grid_index, self%coef)
   
     allocate(self%target_rank(self%index_size))
     call get_target_grid_rank(recv_comp_name, recv_grid_name, self%recv_grid_index, self%target_rank)
  else
     my_grid => get_grid_ptr(trim(recv_grid_name))
     call make_local_mapping_table_no_sort(recv_grid, send_grid, coef, my_grid%get_grid_index_ptr(), &
                                   self%index_size, self%recv_grid_index, self%send_grid_index, self%coef)
     allocate(self%target_rank(self%index_size))
     call get_target_grid_rank(send_comp_name, send_grid_name, self%send_grid_index, self%target_rank)
  end if


  
  call put_log("jlt_exchange_class : set_mapping_table, make_local_mapping_table  end")


  !write(300+source_comp_id*100 + my_grid%get_my_rank(), *) "global mapping table"
  !write(300+source_comp_id*100 + my_grid%get_my_rank(), *) "send_comp = "//trim(send_comp_name)//", recv_comp = "//trim(recv_comp_name)

  !do i = 1, size(send_grid)
  !   write(300+source_comp_id*100 + my_grid%get_my_rank(), *) send_grid(i), recv_grid(i), coef(i)
  !end do

  my_grid_index => my_grid%get_grid_index_ptr()

  !write(300+source_comp_id*100 + my_grid%get_my_rank(), *) "my grid_index"
  !do i = 1, size(my_grid_index)
  !   write(300+source_comp_id*100 + my_grid%get_my_rank(), *) my_grid_index(i)
  !end do
  
  !write(300+source_comp_id*100 + my_grid%get_my_rank(), *) "local mapping table,", self%index_size
  !do i = 1, self%index_size
  !   write(300+source_comp_id*100 + my_grid%get_my_rank(), *) self%send_grid_index(i), self%recv_grid_index(i), self%coef(i), self%target_rank(i)
  !end do
  

  call put_log("jlt_exchange_class : set_mapping_table, reorder_index_by_target_rank start")

  call reorder_index_by_target_rank(self%send_grid_index, self%recv_grid_index, self%coef, self%target_rank)

  call put_log("jlt_exchange_class : set_mapping_table, reorder_index_by_target_rank end")

  call put_log("jlt_exchange_class : set_mapping_table, make_exchange_table start")

  if (trim(self%my_name) == trim(send_comp_name)) then
     call make_exchange_table(self%target_rank, self%send_grid_index, &
                              self%ex_map%num_of_exchange_rank, self%ex_map%exchange_rank, &
                              self%ex_map%num_of_exchange, self%ex_map%offset, self%ex_map%exchange_index) 
     self%ex_map%exchange_data_size = size(self%ex_map%exchange_index)
  else
     call make_exchange_table(self%target_rank, self%send_grid_index, &
                              self%ex_map%num_of_exchange_rank, self%ex_map%exchange_rank, &
                              self%ex_map%num_of_exchange, self%ex_map%offset, self%ex_map%exchange_index) 
     self%ex_map%exchange_data_size = size(self%ex_map%exchange_index)
  end if

  call put_log("jlt_exchange_class : set_mapping_table, make_exchange_table end")

  call put_log("jlt_exchange_class : set_mapping_table, make_conversion_table start")

  
  if (trim(self%my_name) ==  trim(send_comp_name)) then
    call put_log("jlt_exchange_class : set_mapping_table, make_conversion_table 1")
    call make_grid_conv_table(my_grid%get_grid_index_ptr(), self%ex_map%exchange_index, self%ex_map%conv_table)
    call put_log("jlt_exchange_class : set_mapping_table, make_conversion_table 2")
    call send_send_grid_index(self, my_grid%get_grid_index_ptr())   
    call put_log("jlt_exchange_class : set_mapping_table, make_conversion_table 3")
  else
    call put_log("jlt_exchange_class : set_mapping_table, make_conversion_table 4")
    call make_recv_map_info(self%recv_grid_index, self%my_map)
    call put_log("jlt_exchange_class : set_mapping_table, make_conversion_table 5")
    call make_grid_conv_table(my_grid%get_grid_index_ptr(), self%my_map%intpled_index, self%my_map%conv_table)
    call put_log("jlt_exchange_class : set_mapping_table, make_conversion_table 6")
    call recv_send_grid_index(self)
    call put_log("jlt_exchange_class : set_mapping_table, make_conversion_table 7")
    call make_conversion_table(self%recv_grid_index, self%my_map%intpled_index, self%recv_conv_table)

    call sort_conversion_table(self%send_grid_index, self%send_conv_table, self%recv_grid_index, self%recv_conv_table, self%coef) ! add 20250518/mdf 20250604
    call reorder_conversion_table(self%send_grid_index, self%send_conv_table, self%recv_conv_table, self%coef) ! add 20250604
    
    
    do i = 1, size(self%recv_grid_index)
       write(600, *) self%recv_grid_index(i), self%send_grid_index(i)
    end do
    
    do i = 1, size(self%recv_conv_table)
       write(600, *) self%recv_conv_table(i), self%send_conv_table(i), self%coef(i)
    end do
    
    call put_log("jlt_exchange_class : set_mapping_table, make_conversion_table 8")
  end if

  call put_log("jlt_exchange_class : set_mapping_table, make_conversion_table end")

  !write(300+source_comp_id*100 + my_grid%get_my_rank(), *) "exchange rank"
  !do i = 1, self%ex_map%num_of_exchange_rank
  !   write(300+source_comp_id*100 + my_grid%get_my_rank(), *) self%ex_map%exchange_rank(i), self%ex_map%num_of_exchange(i)
  !end do

  !write(300+source_comp_id*100 + my_grid%get_my_rank(), *) "exchange index"
  !do i = 1, size(self%ex_map%exchange_index)
  !   write(300+source_comp_id*100 + my_grid%get_my_rank(), *) self%ex_map%exchange_index(i)
  !end do

  return
  
  call put_log(" set_mapping_table : exchange map info start")


  call put_log(" set_mapping_table : exchange map info end")

  call put_log("jlt_exchange_class : set_mapping_table end")

  call put_log("------------------------------------------------------------------------------------------")

end subroutine set_mapping_table_recv_intpl

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine make_grid_conv_table(grid_index, exchange_index, conv_table)
  use mpi
  use jlt_utils, only : sort_int_1d, binary_search
  implicit none
  integer, intent(IN)  :: grid_index(:)
  integer, intent(IN)  :: exchange_index(:)
  integer, pointer     :: conv_table(:)
  integer, allocatable :: sorted_index(:)
  integer, allocatable :: sorted_pos(:)
  integer              :: res
  integer              :: i, counter

  integer :: my_rank, ierr
  
  if (size(exchange_index) <= 0) return

  call mpi_comm_rank(MPI_COMM_WORLD, my_rank, ierr)

  
  allocate(conv_table(size(exchange_index)))

  allocate(sorted_index(size(grid_index)))
  allocate(sorted_pos(size(grid_index)))

  do i = 1, size(grid_index)
     sorted_index(i) = grid_index(i)
     sorted_pos(i)   = i
  end do
  
  call sort_int_1d(size(grid_index), sorted_index, sorted_pos)

  conv_table(:) = 0
  do i = 1, size(exchange_index)
     res = binary_search(sorted_index, exchange_index(i))
     if (res > 0) conv_table(i) = sorted_pos(res)
  end do

  !write(399 + my_rank, *) "make_grid_conf_table, "
  !do i = 1, size(grid_index)
  !   write(300 + my_rank, *) grid_index(i), exchange_index(i), conv_table(i)
  !end do
  
end subroutine make_grid_conv_table
  
!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine make_grid_conv_table_recv(grid_index, exchange_index, conv_table)
  use mpi
  use jlt_utils, only : sort_int_1d, binary_search
  implicit none
  integer, intent(IN)  :: grid_index(:)
  integer, intent(IN)  :: exchange_index(:)
  integer, intent(INOUT) :: conv_table(:)
  integer, allocatable :: sorted_index(:)
  integer, allocatable :: sorted_pos(:)
  integer              :: res
  integer              :: i, counter

  integer :: my_rank, ierr
  
  if (size(exchange_index) <= 0) return

  call mpi_comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  
  allocate(sorted_index(size(grid_index)))
  allocate(sorted_pos(size(grid_index)))

  do i = 1, size(grid_index)
     sorted_index(i) = grid_index(i)
     sorted_pos(i)   = i
  end do
  
  call sort_int_1d(size(grid_index), sorted_index, sorted_pos)

  conv_table(:) = 0
  do i = 1, size(exchange_index)
     res = binary_search(sorted_index, exchange_index(i))
     if (res > 0) conv_table(i) = sorted_pos(res)
  end do
  
  !write(399 + my_rank, *) "make_grid_conf_table, "
  do i = 1, size(grid_index)
     !write(399 + my_rank, *) grid_index(i), exchange_index(i), conv_table(i)
  end do

end subroutine make_grid_conv_table_recv
  
!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine make_send_map_info(send_grid_index, send_map)
  use jlt_remapping, only : delete_same_index
  implicit none
  integer, pointer                   :: grid_index(:)       ! data grid index
  integer, intent(IN)                :: send_grid_index(:)  ! local interpolation index
  type(send_map_info), intent(INOUT) :: send_map

  call delete_same_index(send_grid_index, send_map%intpled_index)
  send_map%intpled_data_size = size(send_map%intpled_index)
  
end subroutine make_send_map_info

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine make_recv_map_info(recv_grid_index, recv_map)
  use jlt_remapping, only : delete_same_index
  implicit none
  integer, pointer                   :: grid_index(:)       ! data grid index
  integer, intent(IN)                :: recv_grid_index(:)  ! local interpolation index
  type(recv_map_info), intent(INOUT) :: recv_map

  call delete_same_index(recv_grid_index, recv_map%intpled_index)
  recv_map%intpled_data_size = size(recv_map%intpled_index)
  
end subroutine make_recv_map_info

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine send_recv_grid_index(self, grid_index)
  use jlt_mpi_lib, only : jml_IsendModel2, jml_send_waitall, jml_GetMyrank
  implicit none
  class(exchange_class) :: self
  integer, pointer :: grid_index(:)  ! my grid index
  integer, pointer :: exchange_buffer(:)
  integer, pointer :: send_data_ptr
  integer          :: data_tag
  integer          :: i, j
  
  allocate(exchange_buffer(size(self%ex_map%exchange_index)))
  do i = 1, size(self%ex_map%exchange_index)
     exchange_buffer(i) = grid_index(self%recv_conv_table(i))
  end do

  do i = 1, self%ex_map%num_of_exchange_rank
     send_data_ptr => exchange_buffer(self%ex_map%offset(i)+1)
     !data_tag = jml_GetMyrank(self%recv_comp_id)*100000+self%ex_map%exchange_rank(i)
     data_tag = jml_GetMyrank(self%recv_comp_id)*0+self%ex_map%exchange_rank(i)
     call jml_IsendModel2(self%recv_comp_id, send_data_ptr, 1, self%ex_map%num_of_exchange(i), &
                          self%send_comp_id, self%ex_map%exchange_rank(i), data_tag)
  end do

  call jml_send_waitall()

  deallocate(exchange_buffer)
  
end subroutine send_recv_grid_index

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine recv_recv_grid_index(self)
  use jlt_mpi_lib, only : jml_IrecvModel2, jml_recv_waitall, jml_GetMyrank
  use jlt_remapping, only : make_conversion_table
  implicit none
  class(exchange_class) :: self
  integer, pointer :: exchange_buffer(:)
  integer, pointer :: recv_data_ptr
  integer          :: data_tag
  integer          :: i
  
  allocate(exchange_buffer(size(self%ex_map%exchange_index)))

  exchange_buffer(:) = 0

  do i = 1, self%ex_map%num_of_exchange_rank
     recv_data_ptr => exchange_buffer(self%ex_map%offset(i)+1)
     !data_tag = self%ex_map%exchange_rank(i)*100000+jml_GetMyrank(self%send_comp_id)
     data_tag = self%ex_map%exchange_rank(i)*0+jml_GetMyrank(self%send_comp_id)
     call jml_IrecvModel2(self%send_comp_id, recv_data_ptr, 1, self%ex_map%num_of_exchange(i), &
                          self%recv_comp_id, self%ex_map%exchange_rank(i), data_tag)
  end do

  call jml_recv_waitall()

  !write(0, *) "recv_recv_grid_index, ", self%recv_grid_index, exchange_buffer
  
  call make_conversion_table(self%recv_grid_index, exchange_buffer, self%recv_conv_table)
  
  deallocate(exchange_buffer)
  
end subroutine recv_recv_grid_index

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine send_send_grid_index(self, grid_index)
  use jlt_mpi_lib, only : jml_IsendModel2, jml_send_waitall, jml_GetMyrank
  implicit none
  class(exchange_class) :: self
  integer, pointer :: grid_index(:)  ! my grid index
  integer, pointer :: exchange_buffer(:)
  integer, pointer :: send_data_ptr
  integer          :: data_tag
  integer          :: i, j
  
  allocate(exchange_buffer(size(self%ex_map%exchange_index)))
  do i = 1, size(self%ex_map%exchange_index)
     exchange_buffer(i) = grid_index(self%ex_map%conv_table(i))
  end do

  do i = 1, self%ex_map%num_of_exchange_rank
     send_data_ptr => exchange_buffer(self%ex_map%offset(i)+1)
     !data_tag = jml_GetMyrank(self%send_comp_id)*100000+self%ex_map%exchange_rank(i)
     data_tag = jml_GetMyrank(self%send_comp_id)*0+self%ex_map%exchange_rank(i)
     call jml_IsendModel2(self%send_comp_id, send_data_ptr, 1, self%ex_map%num_of_exchange(i), &
                          self%recv_comp_id, self%ex_map%exchange_rank(i), data_tag)
  end do

  call jml_send_waitall()

  deallocate(exchange_buffer)
  
end subroutine send_send_grid_index

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine recv_send_grid_index(self)
  use jlt_mpi_lib, only : jml_IrecvModel2, jml_recv_waitall, jml_GetMyrank
  use jlt_remapping, only : make_conversion_table
  implicit none
  class(exchange_class) :: self
  integer, pointer :: exchange_buffer(:)
  integer, pointer :: recv_data_ptr
  integer          :: data_tag
  integer          :: i
  
  allocate(exchange_buffer(size(self%ex_map%exchange_index)))

  exchange_buffer(:) = 0

  do i = 1, self%ex_map%num_of_exchange_rank
     recv_data_ptr => exchange_buffer(self%ex_map%offset(i)+1)
     !data_tag = self%ex_map%exchange_rank(i)*100000+jml_GetMyrank(self%recv_comp_id)
     data_tag = self%ex_map%exchange_rank(i)*0+jml_GetMyrank(self%recv_comp_id)
     call jml_IrecvModel2(self%recv_comp_id, recv_data_ptr, 1, self%ex_map%num_of_exchange(i), &
                          self%send_comp_id, self%ex_map%exchange_rank(i), data_tag)
  end do

  call jml_recv_waitall()

  call make_conversion_table(self%send_grid_index, exchange_buffer, self%send_conv_table)
  
  deallocate(exchange_buffer)
  
end subroutine recv_send_grid_index

!=======+=========+=========+=========+=========+=========+=========+=========+
!Sort other tables using the values in the recv_conv_table as keys.
subroutine sort_conversion_table(send_index_table, send_conv_table, recv_index_table, recv_conv_table, coef)
  use jlt_utils, only : sort_int_1d, binary_search
  implicit none
  integer, pointer :: send_index_table(:)
  integer, pointer :: send_conv_table(:)
  integer, pointer :: recv_index_table(:)
  integer, pointer :: recv_conv_table(:)
  real(kind=8), pointer :: coef(:)
  integer, allocatable :: sorted_index(:), sorted_pos(:)
  real(kind=8), allocatable :: double_temp(:)
  integer :: i
  
  allocate(sorted_index(size(recv_conv_table)))
  allocate(sorted_pos(size(recv_conv_table)))

  do i = 1, size(recv_conv_table)
     sorted_index(i) = recv_conv_table(i)
     sorted_pos(i)   = i
  end do
  
  call sort_int_1d(size(recv_conv_table), sorted_index, sorted_pos)

  do i = 1, size(recv_conv_table)
     recv_conv_table(i) = sorted_index(i)
  end do
  
  sorted_index(:) = send_conv_table(:)

  do i = 1, size(recv_conv_table)
     send_conv_table(i) = sorted_index(sorted_pos(i))
  end do

  sorted_index(:) = send_index_table(:)

  do i = 1, size(recv_conv_table)
     send_index_table(i) = sorted_index(sorted_pos(i))
  end do

  sorted_index(:) = recv_index_table(:)

  do i = 1, size(recv_conv_table)
     recv_index_table(i) = sorted_index(sorted_pos(i))
  end do

  allocate(double_temp(size(recv_conv_table)))

  double_temp(:) = coef(:)
  do i = 1, size(recv_conv_table)
     coef(i) = double_temp(sorted_pos(i))
  end do
  
  deallocate(sorted_index, sorted_pos, double_temp)     

end subroutine sort_conversion_table

!=======+=========+=========+=========+=========+=========+=========+=========+
!For computations on identical grid points on the recv_conv_table (which are already sorted),
!sort other tables based on the global grid point indices on the sending side (send_index_table).
subroutine reorder_conversion_table(send_index_table, send_conv_table, recv_conv_table, coef)
  implicit none
  integer, pointer :: send_index_table(:)
  integer, pointer :: send_conv_table(:)
  integer, pointer :: recv_conv_table(:)
  real(kind=8), pointer :: coef(:)
  integer, pointer :: send_index_ptr(:)
  integer, pointer :: send_conv_ptr(:)
  real(kind=8), pointer :: coef_ptr(:)

  integer :: current_index
  integer :: index_counter
  integer :: i

  current_index = recv_conv_table(1)
  index_counter = 0
  do i = 2, size(recv_conv_table)
     if (current_index /= recv_conv_table(i)) then
        send_index_ptr => send_index_table(i-1-index_counter:i-1)
        send_conv_ptr => send_conv_table(i-1-index_counter:i-1)
        coef_ptr => coef(i-1-index_counter:i-1)
        call resort_conversion_table(send_index_ptr, send_conv_ptr, coef_ptr)
        current_index = recv_conv_table(i)
        index_counter = 0
     else
        index_counter = index_counter + 1
     end if
  end do

  i = size(recv_conv_table)
  if (index_counter > 0) then
        send_index_ptr => send_index_table(i-index_counter:i)
        send_conv_ptr => send_conv_table(i-index_counter:i)
        coef_ptr => coef(i-index_counter:i)
        call resort_conversion_table(send_index_ptr, send_conv_ptr, coef_ptr)
  end if
     
end subroutine reorder_conversion_table

subroutine resort_conversion_table(send_index, send_conv, coef)
  use jlt_utils, only : sort_int_1d, binary_search
  implicit none
  integer, pointer :: send_index(:)
  integer, pointer :: send_conv(:)
  real(kind=8), pointer :: coef(:)
  integer, allocatable :: sorted_index(:), sorted_pos(:)
  real(kind=8), allocatable :: double_temp(:)
  integer :: i

  write(0, *) "resort_conversion_table ", send_index, send_conv
  
  allocate(sorted_index(size(send_index)))
  allocate(sorted_pos(size(send_index)))

  do i = 1, size(send_index)
     sorted_index(i) = send_index(i)
     sorted_pos(i)   = i
  end do
  
  call sort_int_1d(size(send_index), sorted_index, sorted_pos)

  sorted_index(:) = send_conv(:)

  do i = 1, size(send_index)
     send_conv(i) = sorted_index(sorted_pos(i))
  end do

  allocate(double_temp(size(send_index)))

  double_temp(:) = coef(:)
  do i = 1, size(send_index)
     coef(i) = double_temp(sorted_pos(i))
  end do
  
  deallocate(sorted_index, sorted_pos, double_temp)     
  
end subroutine resort_conversion_table

!=======+=========+=========+=========+=========+=========+=========+=========+

logical function is_my_exchange(self, send_comp_name, send_grid_name, &
                                      recv_comp_name, recv_grid_name)
  implicit none
  class(exchange_class)        :: self
  character(len=*), intent(IN) :: send_comp_name, send_grid_name
  character(len=*), intent(IN) :: recv_comp_name, recv_grid_name

  is_my_exchange = .false.

  if (trim(send_comp_name) /= trim(self%send_comp_name)) return
  if (trim(send_grid_name) /= trim(self%send_grid_name)) return
  if (trim(recv_comp_name) /= trim(self%recv_comp_name)) return
  if (trim(recv_grid_name) /= trim(self%recv_grid_name)) return

  is_my_exchange = .true.
  
end function is_my_exchange
  
!=======+=========+=========+=========+=========+=========+=========+=========+

character(len=STR_SHORT) function get_my_name(self)
  implicit none
  class(exchange_class)     :: self

  get_my_name = trim(self%my_name)

end function get_my_name

  
!=======+=========+=========+=========+=========+=========+=========+=========+

function get_exchange_data_size(self) result (res)
  implicit none
  class(exchange_class) :: self
  integer :: res

  if (self%is_my_intpl()) then ! recv side
     if (self%is_send_intpl()) then
        res = self%send_map%intpled_data_size
     else
        res = self%my_map%intpled_data_size
     end if
  else                         ! send side
    res = self%ex_map%exchange_data_size
  end if

end function get_exchange_data_size

!=======+=========+=========+=========+=========+=========+=========+=========+

function get_exchange_buffer_size(self) result (res)
  implicit none
  class(exchange_class) :: self
  integer :: res

  if (self%is_my_intpl()) then ! recv side
    res = self%ex_map%exchange_data_size
  else                         ! send side
    res = self%ex_map%exchange_data_size
  end if

end function get_exchange_buffer_size

!=======+=========+=========+=========+=========+=========+=========+=========+

function get_num_of_target_rank(self) result(res)
  implicit none
  class(exchange_class) :: self
  integer :: res

  res = self%ex_map%num_of_exchange_rank

end function get_num_of_target_rank

!=======+=========+=========+=========+=========+=========+=========+=========+

function get_target_rank(self, target_rank_num) result(res)
  implicit none
  class(exchange_class) :: self
  integer, intent(IN)   :: target_rank_num
  integer :: res

  res = self%ex_map%exchange_rank(target_rank_num)

end function get_target_rank

!=======+=========+=========+=========+=========+=========+=========+=========+

function get_target_array_size(self, target_rank_num) result(res)
  implicit none
  class(exchange_class) :: self
  integer, intent(IN)   :: target_rank_num
  integer :: res

  res = self%ex_map%num_of_exchange(target_rank_num)

end function get_target_array_size

!=======+=========+=========+=========+=========+=========+=========+=========+

function is_my_intpl(self) result (res)
  implicit none
  class(exchange_class) :: self
  logical :: res

  res = self%intpl_flag

end function is_my_intpl

!=======+=========+=========+=========+=========+=========+=========+=========+

function is_send_intpl(self) result (res)
  implicit none
  class(exchange_class) :: self
  logical :: res

  if (trim(self%my_name) == trim(self%send_comp_name)) then
     res = self%intpl_flag
  else
     res = .not.self%intpl_flag
  end if
  
end function is_send_intpl

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine local_2_exchange(self, grid_data, exchange_data)
  use jlt_mpi_lib, only  : jml_GetMyrankGlobal
  implicit none
  class(exchange_class)       :: self
  real(kind=8), intent(IN)    :: grid_data(:)
  real(kind=8), intent(INOUT) :: exchange_data(:)
  integer :: i

  if (self%is_send_intpl()) then
  !write(0, *) "local to exchange ", jml_GetMyrankGlobal(), size(exchange_data), self%get_exchange_data_size(), self%send_map%conv_table, grid_data
  do i = 1, self%get_exchange_data_size()
     exchange_data(i) = grid_data(self%send_map%conv_table(i))
  end do
  else
  !write(0, *) "local to exchange ", jml_GetMyrankGlobal(), size(exchange_data), self%ex_map%exchange_data_size, self%ex_map%conv_table, grid_data
  do i = 1, self%ex_map%exchange_data_size
     exchange_data(i) = grid_data(self%ex_map%conv_table(i))
  end do

  end if

end subroutine local_2_exchange

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine exchange_2_local(self, exchange_data, grid_data)
  use mpi
  implicit none
  class(exchange_class)       :: self
  real(kind=8), intent(IN)    :: exchange_data(:)
  real(kind=8), intent(INOUT) :: grid_data(:)
  integer :: i
  integer :: my_rank, ierr

  call mpi_comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  !write(399 + my_rank, *) "exchange_2_local, ", exchange_data
  !write(399 + my_rank, *) "exchange_2_local, ", self%my_map%intpled_data_size
  
  do i = 1, self%my_map%intpled_data_size
    grid_data(self%my_map%conv_table(i)) = exchange_data(i)
    !write(399 + my_rank, *) exchange_data(i), self%my_map%conv_table(i), grid_data(i)
  end do

end subroutine exchange_2_local

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine target_2_exchange_buffer(self, target_buffer, exchange_buffer)
  implicit none
  class(exchange_class)      :: self
  type(buffer_class), intent(IN) :: target_buffer(:)
  real(kind=8), intent(INOUT)    :: exchange_buffer(:,:)
  integer :: num_of_layer
  integer :: num_of_data
  integer :: offset
  integer :: i

  num_of_layer = target_buffer(1)%num_of_layer
  
  do i = 1, self%ex_map%num_of_exchange_rank
     num_of_data = self%ex_map%num_of_exchange(i)
     offset      = self%ex_map%offset(i)
     exchange_buffer(offset+1:offset+num_of_data, 1:num_of_layer) = target_buffer(i)%buffer(1:num_of_data, 1:num_of_layer)
  end do
  
end subroutine target_2_exchange_buffer


!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine send_data_1d(self, data, exchange_buffer, intpl_tag, exchange_tag)
  use jlt_constant, only : STR_MID
  use jlt_utils, only    : put_log, get_log_level, DETAIL_LOG
  use jlt_mpi_lib, only : jml_IsendModel2, jml_send_waitall
  implicit none
  class(exchange_class)       :: self
  real(kind=8), intent(IN)    :: data(:)
  real(kind=8), pointer       :: exchange_buffer(:,:)
  integer, intent(IN)         :: intpl_tag
  integer, intent(IN)         :: exchange_tag
  real(kind=8), pointer       :: send_data(:,:)
  real(kind=8), pointer       :: intpl_data(:,:)
  character(len=STR_MID)      :: log_str
  integer                     :: target_rank
  integer                     :: num_of_data
  integer                     :: offset
  real(kind=8), pointer       :: data_ptr
  integer :: i
  

  write(log_str,'("  ",A,I5)') "[send_data_1d] send data START, exchange_tag = ", exchange_tag
  call put_log(trim(log_str))

  if (self%is_my_intpl()) then
     call self%interpolate_data(reshape(data,[size(data), 1]), &
                                exchange_buffer, 1, intpl_tag)
  else
    do i = 1, self%ex_map%exchange_data_size
      exchange_buffer(i,1) = data(i)
    end do
  end if
 
  
  if (get_log_level() == DETAIL_LOG) then
     write(log_str,'("  ",A,F15.5,A,F15.5)') "[send_data_1d] min  = ", minval(data), ", max = ", maxval(data)
     call put_log(trim(log_str))
  end if
  
  do i = 1, self%ex_map%num_of_exchange_rank
     target_rank = self%ex_map%exchange_rank(i)
     num_of_data = self%ex_map%num_of_exchange(i)
     offset      = self%ex_map%offset(i)
     data_ptr    => exchange_buffer(offset+1,1)

     if (get_log_level() == DETAIL_LOG) then
        write(log_str,'("  ",A,I8,A,I10)') "[send_data_1d] target_rank = ", target_rank, ", num_of_data = ", num_of_data
        call put_log(trim(log_str))
     end if

     call jml_ISendModel2(self%send_comp_id, data_ptr, 1, num_of_data, 1, 1, &
                          self%recv_comp_id, target_rank, exchange_tag)
  end do


  write(log_str,'("  ",A,I5)') "[send_data_1d] send data END, exchange_tag = ", exchange_tag
  call put_log(trim(log_str))
  
end subroutine send_data_1d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine send_data_2d(self, data, exchange_buffer, num_of_layer, intpl_tag, exchange_tag)
  use jlt_constant, only : STR_MID
  use jlt_utils, only    : put_log, get_log_level, DETAIL_LOG
  use jlt_mpi_lib, only : jml_IsendModel3, jml_send_waitall
  implicit none
  class(exchange_class)       :: self
  real(kind=8), intent(IN)    :: data(:,:)
  type(buffer_class), intent(INOUT) :: exchange_buffer(:)
  integer, intent(IN)         :: num_of_layer
  integer, intent(IN)         :: intpl_tag
  integer, intent(IN)         :: exchange_tag
  real(kind=8), pointer       :: send_data(:,:)
  real(kind=8), pointer       :: intpl_data(:,:)
  character(len=STR_MID)      :: log_str
  integer                     :: target_rank
  integer                     :: num_of_data
  integer                     :: offset
  integer                     :: intpl_buffer_size
  real(kind=8), pointer       :: intpl_buffer(:,:)
  type data_ptr_type
     real(kind=8), pointer       :: data_ptr
  end type data_ptr_type
  type (data_ptr_type), allocatable :: ptr_array(:)
  integer :: i, k

  write(log_str,'("  ",A,I5)') "[send_data_2d] send data START, exchange_tag = ", exchange_tag
  call put_log(trim(log_str))

  if (get_log_level() == DETAIL_LOG) then
     write(log_str,'("  ",A,F15.5,A,F15.5)') "[send_data_2d] min  = ", minval(data), ", max = ", maxval(data)
     call put_log(trim(log_str))
  end if

  allocate(ptr_array(num_of_layer))

  if (self%is_my_intpl()) then ! send interpolation
     intpl_buffer_size = 0
     do i = 1, self%ex_map%num_of_exchange_rank
        intpl_buffer_size = intpl_buffer_size + self%ex_map%num_of_exchange(i)
     end do
     allocate(intpl_buffer(intpl_buffer_size, num_of_layer))
  end if
  
  
    do k = 1, num_of_layer
       if (self%is_my_intpl()) then
           call self%interpolate_data(data, &
                                      intpl_buffer, num_of_layer, intpl_tag)
          !write(0, *) "send_data_2d, send interpolation ", intpl_buffer
       else
          !write(0, *) "send_data_2d, recv interpolation"
          !exchange_buffer(:,k) = data(:,k)
       end if
    end do
    
    do i = 1, self%ex_map%num_of_exchange_rank
         target_rank = self%ex_map%exchange_rank(i)
         num_of_data = self%ex_map%num_of_exchange(i)
         offset      = self%ex_map%offset(i)
         !!!ptr_array(k)%data_ptr    => exchange_buffer(offset + 1, k)

         if (self%is_my_intpl()) then ! send interpolation
           exchange_buffer(i)%buffer(1:num_of_data, 1:num_of_layer) = intpl_buffer(offset+1:offset+num_of_data, 1:num_of_layer)
         else
           exchange_buffer(i)%buffer(1:num_of_data, 1:num_of_layer) = data(offset+1:offset+num_of_data, 1:num_of_layer)
         end if
        
         if (get_log_level() == DETAIL_LOG) then
            write(log_str,'("  ",A,I8,A,I10, A, F15.5, F15.5)') "[send_data_2d] target_rank = ", target_rank, ", num_of_data = ", num_of_data, &
                 ", min max = ", minval(exchange_buffer(i)%buffer), maxval(exchange_buffer(i)%buffer)
            call put_log(trim(log_str))
         end if
         !call jml_ISendModel3(self%send_comp_id, exchange_buffer(offset+1:offset+num_of_data, 1:num_of_layer), 1, num_of_data, 1, num_of_layer, &
         !                     self%recv_comp_id, target_rank, exchange_tag)
         call jml_ISendModel3(self%send_comp_id, exchange_buffer(i)%buffer, 1, num_of_data, 1, num_of_layer, &
                              self%recv_comp_id, target_rank, exchange_tag)
    !!!end do
    end do

    if (self%is_my_intpl()) then
       deallocate(intpl_buffer)
    end if
    
    write(log_str,'("  ",A,I5)') "[send_data_2d] send data END, exchange_tag = ", exchange_tag
    call put_log(trim(log_str))
 
  
end subroutine send_data_2d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine send_data_2d_org(self, data, exchange_buffer, num_of_layer, intpl_tag, exchange_tag)
  use jlt_constant, only : STR_MID
  use jlt_utils, only    : put_log, get_log_level, DETAIL_LOG
  use jlt_mpi_lib, only : jml_IsendModel3, jml_send_waitall
  implicit none
  class(exchange_class)       :: self
  real(kind=8), intent(IN)    :: data(:,:)
  real(kind=8), pointer       :: exchange_buffer(:,:)
  integer, intent(IN)         :: num_of_layer
  integer, intent(IN)         :: intpl_tag
  integer, intent(IN)         :: exchange_tag
  real(kind=8), pointer       :: send_data(:,:)
  real(kind=8), pointer       :: intpl_data(:,:)
  character(len=STR_MID)      :: log_str
  integer                     :: target_rank
  integer                     :: num_of_data
  integer                     :: offset
  type data_ptr_type
     real(kind=8), pointer       :: data_ptr
  end type data_ptr_type
  type (data_ptr_type), allocatable :: ptr_array(:)
  integer :: i, k
  
  write(log_str,'("  ",A,I5)') "[send_data_2d] send data START, exchange_tag = ", exchange_tag
  call put_log(trim(log_str))

  if (get_log_level() == DETAIL_LOG) then
     write(log_str,'("  ",A,F15.5,A,F15.5)') "[send_data_2d] min  = ", minval(data), ", max = ", maxval(data)
     call put_log(trim(log_str))
  end if

  allocate(ptr_array(num_of_layer))

    do k = 1, num_of_layer
       if (self%is_my_intpl()) then
       else
          exchange_buffer(:,k) = data(:,k)
       end if

    end do
    
    do i = 1, self%ex_map%num_of_exchange_rank
         target_rank = self%ex_map%exchange_rank(i)
         num_of_data = self%ex_map%num_of_exchange(i)
         offset      = self%ex_map%offset(i)
         !!!ptr_array(k)%data_ptr    => exchange_buffer(offset + 1, k)

         if (get_log_level() == DETAIL_LOG) then
            write(log_str,'("  ",A,I8,A,I10, A, F15.5, F15.5)') "[send_data_2d] target_rank = ", target_rank, ", num_of_data = ", num_of_data, &
                 ", min max = ", minval(exchange_buffer(offset+1:offset+num_of_data, 1:num_of_layer)), &
                 maxval(exchange_buffer(offset+1:offset+num_of_data, 1:num_of_layer))
            call put_log(trim(log_str))
         end if
         call jml_ISendModel3(self%send_comp_id, exchange_buffer(offset+1:offset+num_of_data, 1:num_of_layer), 1, num_of_data, 1, num_of_layer, &
                              self%recv_comp_id, target_rank, exchange_tag)
         !!!call jml_ISendModel2(self%send_comp_id, ptr_array(k)%data_ptr, 1, num_of_data, 1, 1, &
         !!!                    self%recv_comp_id, target_rank, exchange_tag)
    !!!end do
    end do

    write(log_str,'("  ",A,I5)') "[send_data_2d] send data END, exchange_tag = ", exchange_tag
    call put_log(trim(log_str))
 
  
end subroutine send_data_2d_org

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine recv_data_1d(self, exchange_buffer, exchange_tag)
  use jlt_constant, only : STR_MID
  use jlt_utils, only    : put_log, get_log_level, DETAIL_LOG
  use jlt_mpi_lib, only : jml_IrecvModel2, jml_recv_waitall
  implicit none
  class(exchange_class)       :: self
  real(kind=8), pointer       :: exchange_buffer(:,:)
  integer, intent(IN)         :: exchange_tag
  character(len=STR_MID)      :: log_str
  integer                     :: target_rank
  integer                     :: num_of_data
  integer                     :: offset
  real(kind=8), pointer       :: data_ptr
  integer :: i
  
  write(log_str,'("  ",A,I5)') "[recv_data_1d] recv data START, exchange_tag = ", exchange_tag
  call put_log(trim(log_str))

     do i = 1, self%ex_map%num_of_exchange_rank
       target_rank = self%ex_map%exchange_rank(i)
       num_of_data = self%ex_map%num_of_exchange(i)
       offset      = self%ex_map%offset(i)
       data_ptr    => exchange_buffer(offset+1, 1)

       if (get_log_level() == DETAIL_LOG) then
          write(log_str,'("  ",A,I8,A,I10)') "[recv_data_1d] target_rank = ", target_rank, ", num_of_data = ", num_of_data
          call put_log(trim(log_str))
       end if

       call jml_IrecvModel2(self%recv_comp_id, data_ptr, 1, num_of_data, 1, 1, &
                           self%send_comp_id, target_rank, exchange_tag)
     end do

  write(log_str,'("  ",A,I5)') "[recv_data_1d] recv data END, exchange_tag = ", exchange_tag
  call put_log(trim(log_str))
  !if (get_log_level() == DETAIL_LOG) then
  !   write(log_str,'("  ",A,F15.5,A,F15.5)') "[recv_data_1d] min = ", minval(exchange_buffer(:,1)), &
  !                                                        ", max = ", maxval(exchange_buffer(:,1))
  !   call put_log(trim(log_str))
  !end if
  
end subroutine recv_data_1d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine recv_data_2d(self, exchange_target, num_of_layer, exchange_tag)
  use jlt_constant, only : STR_MID
  use jlt_utils, only    : put_log, get_log_level, DETAIL_LOG
  use jlt_mpi_lib, only : jml_IrecvModel3, jml_recv_waitall
  implicit none
  class(exchange_class)       :: self
  type(buffer_class), intent(INOUT) :: exchange_target(:)
  integer, intent(IN)         :: num_of_layer
  integer, intent(IN)         :: exchange_tag
  character(len=STR_MID)      :: log_str
  integer                     :: target_rank
  integer                     :: num_of_data
  integer                     :: offset
  type data_ptr_type
     real(kind=8), pointer    :: data_ptr
  end type data_ptr_type
  type (data_ptr_type), allocatable  :: ptr_array(:) 
  integer :: i, k
  
  write(log_str,'("  ",A,I5)') "[recv_data_2d] recv data START, exchange_tag = ", exchange_tag
  call put_log(trim(log_str))

  if (self%intpl_flag) then ! my interpolation


  else

  end if

  allocate(ptr_array(num_of_layer))
  
  do i = 1, self%ex_map%num_of_exchange_rank
       target_rank = self%ex_map%exchange_rank(i)
       num_of_data = self%ex_map%num_of_exchange(i)
       offset      = self%ex_map%offset(i)
       !!!do k = 1, num_of_layer
       !!!ptr_array(k)%data_ptr    => exchange_buffer(offset+1, k)

       if (get_log_level() == DETAIL_LOG) then
          write(log_str,'("  ",A,I8,A,I10)') "[recv_data_2d] target_rank = ", target_rank, ", num_of_data = ", num_of_data
          call put_log(trim(log_str))
       end if

       call jml_IrecvModel3(self%recv_comp_id, exchange_target(i)%buffer, 1, num_of_data, 1, num_of_layer,  &
                            self%send_comp_id, target_rank, exchange_tag)
       !exchange_buffer(offset+1:offset+num_of_data, 1:num_of_layer) = exchange_target(i)%buffer(1:num_of_data, 1:num_of_layer)
       !!!call jml_IrecvModel2(self%recv_comp_id, ptr_array(k)%data_ptr, 1, num_of_data, 1,1,  &
       !!!                     self%send_comp_id, target_rank, exchange_tag)
       !!!end do
  end do

  deallocate(ptr_array)
  
  write(log_str,'("  ",A,I5)') "[recv_data_2d] recv data END, exchange_tag = ", exchange_tag
  call put_log(trim(log_str))
  !if (get_log_level() == DETAIL_LOG) then
  !   write(log_str,'("  ",A,F15.5,A,F15.5)') "[recv_data_2d] min = ", minval(exchange_buffer(:,1:num_of_layer)), &
  !                                                        ", max = ", maxval(exchange_buffer(:,1:num_of_layer))
  !   call put_log(trim(log_str))
  !end if
  

end subroutine recv_data_2d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine recv_data_2d_org(self, exchange_buffer, num_of_layer, exchange_tag)
  use jlt_constant, only : STR_MID
  use jlt_utils, only    : put_log, get_log_level, DETAIL_LOG
  use jlt_mpi_lib, only : jml_IrecvModel3, jml_recv_waitall
  implicit none
  class(exchange_class)       :: self
  real(kind=8), pointer       :: exchange_buffer(:,:)
  integer, intent(IN)         :: num_of_layer
  integer, intent(IN)         :: exchange_tag
  character(len=STR_MID)      :: log_str
  integer                     :: target_rank
  integer                     :: num_of_data
  integer                     :: offset
  type data_ptr_type
     real(kind=8), pointer    :: data_ptr
  end type data_ptr_type
  type (data_ptr_type), allocatable  :: ptr_array(:) 
  integer :: i, k
  
  write(log_str,'("  ",A,I5)') "[recv_data_2d] recv data START, exchange_tag = ", exchange_tag
  call put_log(trim(log_str))

  if (self%intpl_flag) then ! my interpolation


  else

  end if

  allocate(ptr_array(num_of_layer))
  
  do i = 1, self%ex_map%num_of_exchange_rank
       target_rank = self%ex_map%exchange_rank(i)
       num_of_data = self%ex_map%num_of_exchange(i)
       offset      = self%ex_map%offset(i)
       !!!do k = 1, num_of_layer
       !!!ptr_array(k)%data_ptr    => exchange_buffer(offset+1, k)

       if (get_log_level() == DETAIL_LOG) then
          write(log_str,'("  ",A,I8,A,I10)') "[recv_data_2d] target_rank = ", target_rank, ", num_of_data = ", num_of_data
          call put_log(trim(log_str))
       end if

       call jml_IrecvModel3(self%recv_comp_id, exchange_buffer(offset+1:offset+num_of_data, 1:num_of_layer), 1, num_of_data, 1, num_of_layer,  &
                            self%send_comp_id, target_rank, exchange_tag)
       !!!call jml_IrecvModel2(self%recv_comp_id, ptr_array(k)%data_ptr, 1, num_of_data, 1,1,  &
       !!!                     self%send_comp_id, target_rank, exchange_tag)
       !!!end do
  end do

  deallocate(ptr_array)
  
  write(log_str,'("  ",A,I5)') "[recv_data_2d] recv data END, exchange_tag = ", exchange_tag
  call put_log(trim(log_str))
  !if (get_log_level() == DETAIL_LOG) then
  !   write(log_str,'("  ",A,F15.5,A,F15.5)') "[recv_data_2d] min = ", minval(exchange_buffer(:,1:num_of_layer)), &
  !                                                        ", max = ", maxval(exchange_buffer(:,1:num_of_layer))
  !   call put_log(trim(log_str))
  !end if
  

end subroutine recv_data_2d_org

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine buffer_2_recv_data(self, exchange_buffer, recv_data)
  use mpi
  implicit none
  class(exchange_class)       :: self
  real(kind=8), intent(IN)    :: exchange_buffer(:,:)
  real(kind=8), intent(INOUT) :: recv_data(:,:)
  integer :: recv_index
  integer :: i, k
  integer :: my_rank, ierr

  call mpi_comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  !write(399+my_rank, *) "buffer_2_recv_data"
  do k = 1, size(recv_data, 2)
     do i = 1, size(exchange_buffer, 1)
        recv_index = self%recv_conv_table(i)
        recv_data(recv_index, k) = recv_data(recv_index, k) + exchange_buffer(i, k)
        !write(399 + my_rank, *) exchange_buffer(i, k), recv_index, recv_data(recv_index, k)
     end do
  end do
  
  do k = 1, size(recv_data, 2)
     do i = 1, size(exchange_buffer, 1)
        !write(399 + my_rank, *) exchange_buffer(i, k), recv_data(i, k)
     end do
  end do

end subroutine buffer_2_recv_data

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine interpolate_data(self, send_data, recv_data, num_of_layer, intpl_tag)
  use jlt_grid, only : get_grid_ptr
  use jlt_utils, only : put_log
  use mpi
  implicit none
  class(exchange_class)       :: self
  real(kind=8), intent(IN)    :: send_data(:,:)
  real(kind=8), intent(INOUT) :: recv_data(:,:)
  integer, intent(IN)         :: num_of_layer
  integer, intent(IN)         :: intpl_tag
  integer, pointer            :: send_index(:), recv_index(:)
  integer                     :: send_grid
  integer                     :: recv_grid
  real(kind=8), allocatable   :: weight_data(:)
  integer, allocatable        :: check_data(:)
  real(kind=8) :: missing_value
  integer :: i, k, n
  type(grid_class), pointer   :: my_grid
  character(len=STR_MID)      :: log_str
  integer :: my_rank

  call MPI_comm_rank(MPI_COMM_WORLD, my_rank, i)
  
  write(log_str,'("  ",A,I5)') "[interpolate_data] interpolate data START, intpl_tag = ", intpl_tag
  call put_log(trim(log_str))

  my_grid => get_grid_ptr(trim(self%recv_grid_name))
  
  recv_data(:,:) = 0.d0

  send_index => self%send_conv_table
  recv_index => self%recv_conv_table
  
  !write(300+self%recv_comp_id*100 + my_grid%get_my_rank(), *) "interpolate_data ", num_of_layer, self%my_map%get_intpl_size(), size(recv_index), size(send_index)
  !write(300+self%recv_comp_id*100 + my_grid%get_my_rank(), *) "interpolate_data ", size(recv_data, 1), size(send_data,1), minval(send_index), maxval(send_index),&
   !                                minval(recv_index), maxval(recv_index), minval(self%coef), maxval(self%coef)

  if (intpl_tag > 9800) then
      if ((10000 >= intpl_tag).and.(intpl_tag > 9900)) then
         missing_value = -999.d0
      else if ((9900 >= intpl_tag).and.(intpl_tag > 9800)) then
         missing_value = -1.d0
      endif

      allocate(weight_data(size(recv_data(:,1))))
      allocate(check_data(size(recv_data(:,1))))
      do n = 1, num_of_layer
        weight_data(:) = 0.d0
        check_data(:) = 0
        do i = 1, size(self%coef)
          send_grid = send_index(i)
           recv_grid = recv_index(i)
           if (send_data(send_grid, n) == missing_value) then
              check_data(recv_grid) = 1
           else
              recv_data(recv_grid, n) = recv_data(recv_grid, n) + send_data(send_grid, n)*self%coef(i)
              weight_data(recv_grid) = weight_data(recv_grid) + self%coef(i)
           end if
        end do
        do i = 1, size(recv_data(:,1))
           if (weight_data(i) > 0.d0) then
              recv_data(i, n) = recv_data(i, n) / weight_data(i)
           else if (check_data(i) == 1) then
              recv_data(i, n) = missing_value
           end if
        end do
      end do

      deallocate(weight_data)
      deallocate(check_data)
  else ! scalar data
    !write(300+self%recv_comp_id*100 + my_grid%get_my_rank(), *) "interpolation"
         if (intpl_tag == 14) then
         write(300+my_grid%get_my_rank(), *) "interpolate_data test write"
         write(300+my_grid%get_my_rank(), *) size(self%coef)
         end if
    do k = 1, num_of_layer
       do i = 1, size(self%coef)
         recv_data(recv_index(i),k) = recv_data(recv_index(i),k) + send_data(send_index(i),k)*self%coef(i)
         if (intpl_tag == 14) then
         write(300+my_grid%get_my_rank(), *) i, send_index(i), recv_index(i), send_data(send_index(i),k), recv_data(recv_index(i),k), self%coef(i)
         end if
         !!write(my_rank+300, *) i, send_index(i), recv_index(i), send_data(send_index(i),k), recv_data(recv_index(i),k), self%coef(i)
         !write(300+self%recv_comp_id*100 + my_grid%get_my_rank(), *) i, send_index(i), recv_index(i), send_data(send_index(i),k), recv_data(recv_index(i),k), self%coef(i)
      end do
    end do
 end if
 
  write(log_str,'("  ",A,I5)') "[interpolate_data] interpolate data END"
  call put_log(trim(log_str))

  !write(0, *) "interpolate_data min, max = ", minval(send_data), maxval(send_data), minval(recv_data), maxval(recv_data)
end subroutine interpolate_data

!=======+=========+=========+=========+=========+=========+=========+=========+
end module jlt_exchange_class
