module jlt_map_class
  use jlt_constant, only : STR_SHORT, STR_MID, STR_LONG
  use jlt_grid_class, only : grid_class
  
!--------------------------------   public  ----------------------------------!

public :: map_class

!--------------------------------   private  ---------------------------------!


type target_map_info
   integer              :: whole_data_size          ! sum of my_data_size (significant only at root)
   integer, pointer     :: recvcounts(:)            ! array of my_data_size of earch process (significant only at root)
   integer, pointer     :: displs(:)                ! displacement for gatherv, scatterv (significant only at root)
   integer              :: my_data_size             ! data size fo my send/recv data = sum of num_of_data
   integer              :: num_of_pe                ! number of my target process
   integer, pointer     :: offset(:)                ! offset of send/recv array
   integer, pointer     :: num_of_data(:)           ! size of each send/recv data
   integer, pointer     :: target_rank(:)           ! proess number of my target
end type target_map_info

type map_class
   private
   type(grid_class), pointer :: my_grid
   integer              :: comp_id                  ! my component id
   logical              :: send_recv_flag = .true.  ! send map or recv map (.true. = send map)
   integer              :: num_of_exchange_index    ! local size of exchange data array
   integer, allocatable :: exchange_index(:)        ! global index of my exchange grid points
   integer, allocatable :: conv_table(:)            ! array conversion table 
   integer, allocatable :: recvcounts(:)            ! array of num_of_exchange_index of each process (significant only at root)
   integer, allocatable :: displs(:)                ! displacement for gatherv, scatterv (significant only at root)
   integer              :: intpl_map_size           ! size of interpolation map  (significant only at root)
   integer              :: intpl_data_size          ! size of interpolation data (significant only at root)
   integer, pointer     :: intpl_index(:)           ! table of interpolation index i to data index (significant only at root)

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
   procedure :: cal_my_exchange_index               ! subroutine (global_map_table)
   procedure :: cal_intpl_index                     ! subroutine (global_map_table)
   procedure :: send_intpl_info                     ! subroutine (my_comp_id, send_comp_id)
   procedure :: recv_intpl_info                     ! subroutine (my_comp_id, recv_comp_id)
   procedure :: local_2_exchange                    ! subroutine (local_data, exchange_data)
   procedure :: exchange_2_local                    ! subroutine (local_data, exchange_data)
   procedure :: get_intpl_size                      ! integer function ()  !return intpl_map_size
   procedure :: get_intpl_index                     ! integer(:), pointer function()
   procedure :: get_data_size                       ! integer function ()  !return intpl_data_size
   procedure :: gather_data_to_intpl                ! subroutine (local_exchange_data, global_exchange_data)
   procedure :: scatter_data_from_intpl             ! subroutine (local_exchange_data, global_exchange_data)

   procedure :: cal_local_to_local_exchange_info    ! subroutine (target_map, is_my_intpl)
   procedure :: get_my_whole_data_size              ! integer function ()                  ! return whole_data_size
   procedure :: get_my_target_data_size             ! integer function ()                  ! return my_data_size
   procedure :: get_num_of_target_pe                ! integer function ()                  ! return num_of_pe
   procedure :: get_target_rank                     ! integer function (pe_num)            ! return target_rank(pe_num)
   procedure :: get_offset                          ! integer function (pe_num)            ! return offset(pe_num)
   procedure :: get_num_of_data                     ! integer function (pe_num)            ! return num_of_data(pe_num)
   procedure :: gather_recv_data_to_intpl           ! subroutine (data, intpl_data)
   procedure :: scatter_interpolated_data_to_send   ! subroutine (send_data, intpl_data)
   
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
  my_map_class%intpl_map_size = 1  ! set default value = 1
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

subroutine cal_my_exchange_index(self, global_map_table)
  use jlt_utils, only : sort_int_1d, binary_search, put_log
  implicit none
  class(map_class)  :: self
  integer, intent(IN)    :: global_map_table(:)
  integer :: i
  integer :: res
  integer, pointer :: sorted_map(:)
  integer, pointer :: grid_index(:)
  character(len=STR_LONG) :: log_str
  
  write(log_str, *) " cal_my_exchange_index start : map_size = ", size(global_map_table),&
                    ", range = ", minval(global_map_table), maxval(global_map_table)
  call put_log(trim(log_str))
  
  allocate(sorted_map(size(global_map_table)))
  sorted_map(:) = global_map_table(:)
  
  call sort_int_1d(size(sorted_map), sorted_map)

  grid_index => self%my_grid%get_grid_index_ptr()

  self%num_of_exchange_index = 0
  do i = 1, self%my_grid%get_grid_size()
     res = binary_search(sorted_map, grid_index(i))
     if (res > 0) then
        self%num_of_exchange_index = self%num_of_exchange_index + 1
     end if
  end do

  write(log_str, *) " cal_my_exchange_index : num_of_exchange_index = ", self%num_of_exchange_index
  call put_log(trim(log_str))

  allocate(self%exchange_index(self%num_of_exchange_index))
  allocate(self%conv_table(self%num_of_exchange_index))

  self%num_of_exchange_index = 0
  do i = 1, self%my_grid%get_grid_size()
     res = binary_search(sorted_map, grid_index(i))
     if (res > 0) then
        self%num_of_exchange_index = self%num_of_exchange_index + 1
        self%exchange_index(self%num_of_exchange_index) = grid_index(i)
        self%conv_table(self%num_of_exchange_index) = i
     end if
  end do
  
  write(log_str, *) " cal_my_exchange_index : cal exchange index OK"
  call put_log(trim(log_str))

  call set_recvcounts(self)

  deallocate(sorted_map)
  
  call put_log("cal_my_exchange_index end")

end subroutine cal_my_exchange_index

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine set_recvcounts(self)
  use jlt_mpi_lib, only : jml_AllGatherLocal
  implicit none
  class(map_class) :: self
  integer :: int_array(1)
  integer :: i

  allocate(self%recvcounts(self%my_grid%get_my_size()))
  allocate(self%displs(self%my_grid%get_my_size()))
  int_array(1)= self%num_of_exchange_index     

  call jml_AllGatherLocal(self%my_grid%get_comp_id(), int_array, 1, 1, self%recvcounts)

  self%displs(1) = 0
  do i = 1, self%my_grid%get_my_size() - 1
     self%displs(i+1) = self%displs(i) + self%recvcounts(i)
  end do

end subroutine set_recvcounts

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine set_recvcounts_old(self)
  use jlt_mpi_lib, only : jml_GatherLocal
  implicit none
  class(map_class) :: self
  integer :: int_array(1)
  integer :: i
  
  if (self%my_grid%get_my_rank() == 0) then
     allocate(self%recvcounts(self%my_grid%get_my_size()))
     allocate(self%displs(self%my_grid%get_my_size()))
     int_array(1)= self%num_of_exchange_index     
  else
     allocate(self%recvcounts(1))
     allocate(self%displs(1))
     int_array(1) = self%num_of_exchange_index
  end if

  call jml_GatherLocal(self%my_grid%get_comp_id(), int_array, 1, 1, self%recvcounts)

  if (self%my_grid%get_my_rank() == 0) then
    self%displs(1) = 0
    do i = 1, self%my_grid%get_my_size() - 1
       self%displs(i+1) = self%displs(i) + self%recvcounts(i)
    end do
 end if

end subroutine set_recvcounts_old

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine cal_intpl_index(self, global_map_table)
  use jlt_mpi_lib, only : jml_GatherVLocal
  use jlt_utils, only   : sort_int_1d, binary_search, put_log, error
  implicit none
  class(map_class) :: self
  integer, intent(IN)   :: global_map_table(:)
  integer, allocatable  :: global_grid_index(:)
  integer, allocatable  :: sorted_index(:)
  character(len=STR_LONG) :: log_str
  integer :: res
  integer :: i

  call put_log("cal_intpl_index start")
  
  self%intpl_map_size  = size(global_map_table)
  self%intpl_data_size = 0
  do i = 1, self%my_grid%get_my_size()
     self%intpl_data_size = self%intpl_data_size + self%recvcounts(i)
  end do

  if (self%my_grid%get_my_rank() == 0) then

     allocate(self%intpl_index(self%intpl_map_size))
     allocate(global_grid_index(self%intpl_data_size))
     allocate(sorted_index(self%intpl_data_size))

     do i = 1, self%intpl_data_size
        sorted_index(i) = i
     end do
     
     call jml_GatherVLocal(self%my_grid%get_comp_id(), self%exchange_index, self%num_of_exchange_index, &
                           global_grid_index, self%recvcounts, self%displs)

     call sort_int_1d(self%intpl_data_size, global_grid_index, sorted_index)

     !write(0, *) "sort result ", size(global_grid_index), minval(global_grid_index), maxval(global_grid_index), &
     !     size(sorted_index), minval(sorted_index), maxval(sorted_index)
     
     do i = 1, self%intpl_map_size
        res = binary_search(global_grid_index, global_map_table(i))
        if ((res <= 0).or.(res > self%intpl_data_size)) then
           write(log_str, *) "table search error, index = ", i, ", map_index = ", global_map_table(i) 
           call error("cal_intpl_index", trim(log_str))
        end if
        self%intpl_index(i) = sorted_index(res)
     end do
     
     deallocate(global_grid_index)
     deallocate(sorted_index)
     write(log_str, *) " interpolation index calculation OK, min, max = ", minval(self%intpl_index), maxval(self%intpl_index)
     call put_log(trim(log_str))
  else
     self%intpl_map_size  = 0
     self%intpl_data_size = 0
     allocate(global_grid_index(1))
     call jml_GatherVLocal(self%my_grid%get_comp_id(), self%exchange_index, self%num_of_exchange_index, &
                           global_grid_index, self%recvcounts, self%displs)
     deallocate(global_grid_index)
  end if

  call put_log("cal_intpl_index end")
  
end subroutine cal_intpl_index

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine cal_intpl_index_old(self, global_map_table)
  use jlt_mpi_lib, only : jml_GatherVLocal
  use jlt_utils, only   : sort_int_1d, binary_search, put_log, error
  implicit none
  class(map_class) :: self
  integer, intent(IN)   :: global_map_table(:)
  integer, allocatable  :: global_grid_index(:)
  integer, allocatable  :: sorted_index(:)
  character(len=STR_LONG) :: log_str
  integer :: res
  integer :: i

  call put_log(" cal_intpl_index start")
  
  if (self%my_grid%get_my_rank() == 0) then
     self%intpl_map_size  = size(global_map_table)
     self%intpl_data_size = 0
     do i = 1, self%my_grid%get_my_size()
        self%intpl_data_size = self%intpl_data_size + self%recvcounts(i)
     end do

     allocate(self%intpl_index(self%intpl_map_size))
     allocate(global_grid_index(self%intpl_data_size))
     allocate(sorted_index(self%intpl_data_size))

     do i = 1, self%intpl_data_size
        sorted_index(i) = i
     end do
     
     call jml_GatherVLocal(self%my_grid%get_comp_id(), self%exchange_index, self%num_of_exchange_index, &
                           global_grid_index, self%recvcounts, self%displs)

     call sort_int_1d(self%intpl_data_size, global_grid_index, sorted_index)

     !write(0, *) "sort result ", size(global_grid_index), minval(global_grid_index), maxval(global_grid_index), &
     !     size(sorted_index), minval(sorted_index), maxval(sorted_index)
     
     do i = 1, self%intpl_map_size
        res = binary_search(global_grid_index, global_map_table(i))
        if ((res <= 0).or.(res > self%intpl_data_size)) then
           call error("cal_intpl_index", "table sorting error")
        end if
        self%intpl_index(i) = sorted_index(res)
     end do
     
     deallocate(global_grid_index)
     deallocate(sorted_index)
     write(log_str, *) " interpolation index calculation OK, min, max = ", minval(self%intpl_index), maxval(self%intpl_index)
     call put_log(trim(log_str))
  else
     self%intpl_map_size  = 0
     self%intpl_data_size = 0
     allocate(global_grid_index(1))
     call jml_GatherVLocal(self%my_grid%get_comp_id(), self%exchange_index, self%num_of_exchange_index, &
                           global_grid_index, self%recvcounts, self%displs)
     deallocate(global_grid_index)
  end if

  call put_log(" cal_intpl_index end")
  
end subroutine cal_intpl_index_old

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine send_intpl_info(self, my_comp_id, dest_comp_id)
  use jlt_mpi_lib, only : jml_isLocalLeader, jml_SendLeader
  use jlt_utils, only : put_log
  implicit none
  class(map_class)    :: self
  integer, intent(IN) :: my_comp_id
  integer, intent(IN) :: dest_comp_id
  integer :: int_array(1)

  call put_log("send_intpl_info start")
  
  if (jml_isLocalLeader(my_comp_id)) then
     int_array(1) = size(self%recvcounts)
     call jml_SendLeader(int_array, 1, 1, dest_comp_id - 1)
     call jml_SendLeader(self%recvcounts, 1, int_array(1), dest_comp_id - 1)
     call jml_SendLeader(self%displs, 1, int_array(1), dest_comp_id - 1)
     int_array(1) = self%intpl_map_size
     call jml_SendLeader(int_array, 1, 1, dest_comp_id - 1)
     call jml_SendLeader(self%intpl_index, 1, self%intpl_map_size, dest_comp_id - 1)
     int_array(1) = self%intpl_data_size
     call jml_SendLeader(int_array, 1, 1, dest_comp_id - 1)
  end if
  
  call put_log("send_intpl_info end")

end subroutine send_intpl_info

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine recv_intpl_info(self, my_comp_id, source_comp_id)
  use jlt_mpi_lib, only : jml_isLocalLeader, jml_RecvLeader, jml_BcastLocal
  use jlt_utils, only : put_log
  implicit none
  class(map_class)    :: self
  integer, intent(IN) :: my_comp_id
  integer, intent(IN) :: source_comp_id
  integer :: int_array(1)

  call put_log("recv_intpl_info start")
  
  if (jml_isLocalLeader(my_comp_id)) then
     call jml_RecvLeader(int_array, 1, 1, source_comp_id - 1)
     call jml_BcastLocal(my_comp_id, int_array, 1,1)
     
     allocate(self%recvcounts(int_array(1)))
     call jml_RecvLeader(self%recvcounts, 1, int_array(1), source_comp_id - 1)
     call jml_BcastLocal(my_comp_id, self%recvcounts, 1,int_array(1))

     allocate(self%displs(int_array(1)))
     call jml_RecvLeader(self%displs, 1, int_array(1), source_comp_id - 1)
     call jml_BcastLocal(my_comp_id, self%displs, 1, int_array(1))

     call jml_RecvLeader(int_array, 1, 1, source_comp_id - 1)
     self%intpl_map_size = int_array(1)
     allocate(self%intpl_index(int_array(1)))
     call jml_RecvLeader(self%intpl_index, 1, self%intpl_map_size, source_comp_id - 1)
     call jml_RecvLeader(int_array, 1, 1, source_comp_id - 1)
     self%intpl_data_size = int_array(1)
  else
     call jml_BcastLocal(my_comp_id, int_array, 1, 1)
     
     allocate(self%recvcounts(int_array(1)))
     call jml_BcastLocal(my_comp_id, self%recvcounts, 1, int_array(1))

     allocate(self%displs(int_array(1)))
     call jml_BcastLocal(my_comp_id, self%displs, 1, int_array(1))
  end if
  
  call put_log("recv_intpl_info end")

end subroutine recv_intpl_info

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function get_intpl_size(self)
  implicit none
  class(map_class) ::self

  get_intpl_size = self%intpl_map_size

end function get_intpl_size

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function get_data_size(self)
  implicit none
  class(map_class) ::self

  get_data_size = self%intpl_data_size

end function get_data_size

!=======+=========+=========+=========+=========+=========+=========+=========+

function get_intpl_index(self) result(res)
  implicit none
  class(map_class) :: self
  integer, pointer :: res(:)

  res => self%intpl_index

end function get_intpl_index

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

subroutine gather_data_to_intpl(self, data, intpl_data)
  use jlt_mpi_lib, only : jml_GatherVLocal
  implicit none
  class(map_class) :: self
  real(kind=8), intent(IN)    :: data(:)
  real(kind=8), intent(INOUT) :: intpl_data(:)

  call jml_GatherVLocal(self%my_grid%get_comp_id(), data, self%num_of_exchange_index, &
                        intpl_data, self%recvcounts, self%displs)

end subroutine gather_data_to_intpl

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine scatter_data_from_intpl(self, data, intpl_data)
  use jlt_mpi_lib, only : jml_ScatterVLocal
  implicit none
  class(map_class) :: self
  real(kind=8), intent(INOUT) :: data(:)
  real(kind=8), intent(IN)    :: intpl_data(:)

  call jml_ScatterVLocal(self%my_grid%get_comp_id(), &
                         intpl_data, self%recvcounts, self%displs, &
                         data, self%num_of_exchange_index)
  
end subroutine scatter_data_from_intpl

!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+
!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine cal_local_to_local_exchange_info(self, tmap, is_my_intpl)
  use jlt_utils, only : put_log, error
  implicit none
  class(map_class)        :: self
  class(map_class)        :: tmap
  logical, intent(IN)     :: is_my_intpl
  integer                 :: npe_self            ! number of my process
  integer                 :: npe_tmap            ! number of target process
  integer                 :: num_of_target_pe    ! number of my target process
  integer                 :: target_rank_start   ! start rank of my target
  integer                 :: target_rank_end     ! end rank of my target
  integer, allocatable    :: npe_my_target(:)
  integer                 :: my_rank
  integer                 :: counter
  integer                 :: i
  character(len=STR_LONG) :: log_str

  call put_log("------------------------------------------------------------------------------------------")

  call put_log("calculate local to local exchange info start")

  allocate(self%tmap_info)
  
  my_rank = self%get_my_rank()
  
  npe_self = size(self%recvcounts)
  npe_tmap = size(tmap%recvcounts)

  if ((is_my_intpl).and.(npe_self > npe_tmap)) then
     call error("cal_local_exchange_info", "currently my npe must <= target npe when my intpl")
  end if

  if (npe_self >= npe_tmap) then
     num_of_target_pe = 1
     allocate(npe_my_target(npe_tmap))
     npe_my_target(:) = int(npe_self/npe_tmap)
     do i = 1, npe_tmap
        if (i-1 < mod(npe_self,npe_tmap)) then
           npe_my_target(i) = npe_my_target(i) + 1
        end if
     end do
     counter = 0
     do i = 1, npe_tmap
        counter = counter + npe_my_target(i)
        if (my_rank < counter) then
           target_rank_start = i - 1
           exit
        end if
     end do
     target_rank_end = target_rank_start
  else
     allocate(npe_my_target(npe_self))
     npe_my_target(:) = int(npe_tmap/npe_self)
     do i = 1, npe_self
        if (i-1 < mod(npe_tmap,npe_self)) then
           npe_my_target(i) = npe_my_target(i) + 1
        end if
     end do

     counter = 0
     do i = 0, my_rank 
        counter = counter + npe_my_target(i+1)
     end do
     target_rank_end = counter - 1
     target_rank_start = counter - npe_my_target(my_rank + 1) 
     num_of_target_pe = target_rank_end - target_rank_start + 1
  end if
  
  self%tmap_info%num_of_pe = num_of_target_pe
  allocate(self%tmap_info%num_of_data(num_of_target_pe))
  allocate(self%tmap_info%offset(num_of_target_pe))
  allocate(self%tmap_info%target_rank(num_of_target_pe))

  self%tmap_info%my_data_size = 0
  if (is_my_intpl) then ! this means that npe_self <= npe_tmap 
    do i = 1, num_of_target_pe
       self%tmap_info%target_rank(i) = target_rank_start + i - 1
       self%tmap_info%num_of_data(i) = tmap%recvcounts(target_rank_start + i)
       self%tmap_info%my_data_size   = self%tmap_info%my_data_size + self%tmap_info%num_of_data(i)
    end do
  else                   ! not interpolation side and npe_self > npe_tmap
    do i = 1, num_of_target_pe
       self%tmap_info%target_rank(i) = target_rank_start + i - 1
       self%tmap_info%num_of_data(i) = self%recvcounts(my_rank + 1)
       self%tmap_info%my_data_size   = self%tmap_info%my_data_size + self%tmap_info%num_of_data(i)
    end do
  end if
 
  self%tmap_info%offset(1) = 0
  do i = 2, num_of_target_pe
     self%tmap_info%offset(i) = self%tmap_info%offset(i-1) + self%tmap_info%num_of_data(i-1)
  end do
  
  deallocate(npe_my_target)

  call put_log("   exchange info")
  write(log_str, *)  "     num of my target process = ", self%tmap_info%num_of_pe
  call put_log(trim(log_str))
  write(log_str, *)  "     total exchange data size = ", self%tmap_info%my_data_size
  call put_log(trim(log_str))
  do i = 1, num_of_target_pe
     write(log_str, '("      target_rank = ",I5, ", data_size = ",I8, ", offset = ", I8)') self%tmap_info%target_rank(i), &
          self%tmap_info%num_of_data(i), self%tmap_info%offset(i)
     call put_log(trim(log_str))
  end do

  call cal_whole_data_info(self)
  
  call put_log("calculate local to local exchange info end")

  call put_log("------------------------------------------------------------------------------------------")

end subroutine cal_local_to_local_exchange_info

  
!=======+=========+=========+=========+=========+=========+=========+=========+

integer function get_num_of_target_pe(self)
  implicit none
  class(map_class) :: self

  get_num_of_target_pe = self%tmap_info%num_of_pe

end function get_num_of_target_pe

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function get_my_whole_data_size(self)
  implicit none
  class(map_class) :: self

  get_my_whole_data_size = self%tmap_info%whole_data_size

end function get_my_whole_data_size

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function get_my_target_data_size(self)
  implicit none
  class(map_class) :: self

  get_my_target_data_size = self%tmap_info%my_data_size

end function get_my_target_data_size

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

subroutine cal_whole_data_info(self)
  use jlt_mpi_lib, only : jml_GatherLocal
  implicit none
  type(map_class) :: self
  integer :: int_array(1)
  integer :: i
  
  if (self%get_my_rank() == 0) then
     allocate(self%tmap_info%recvcounts(size(self%recvcounts)))
     allocate(self%tmap_info%displs(size(self%displs)))
     int_array(1) = self%tmap_info%my_data_size
     call jml_GatherLocal(self%get_comp_id(), int_array, 1, 1, self%tmap_info%recvcounts)

     self%tmap_info%whole_data_size = sum(self%tmap_info%recvcounts)

     self%tmap_info%displs(1) = 0
     do i = 2, size(self%tmap_info%recvcounts)
        self%tmap_info%displs(i) = self%tmap_info%displs(i-1) + self%tmap_info%recvcounts(i-1)
     end do
     
  else
     allocate(self%tmap_info%recvcounts(1))
     allocate(self%tmap_info%displs(1))
     int_array(1) = self%tmap_info%my_data_size
     call jml_GatherLocal(self%get_comp_id(), int_array, 1, 1, self%tmap_info%recvcounts)
  end if
  
end subroutine cal_whole_data_info

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine gather_recv_data_to_intpl(self, data, intpl_data)
  use jlt_mpi_lib, only : jml_GatherVLocal
  implicit none
  class(map_class) :: self
  real(kind=8), intent(IN)    :: data(:)
  real(kind=8), intent(INOUT) :: intpl_data(:)

  call jml_GatherVLocal(self%my_grid%get_comp_id(), data, self%tmap_info%my_data_size, &
                        intpl_data, self%tmap_info%recvcounts, self%tmap_info%displs)

end subroutine gather_recv_data_to_intpl

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine scatter_interpolated_data_to_send(self, send_data, intpl_data)
  use jlt_mpi_lib, only : jml_ScatterVLocal
  implicit none
  class(map_class) :: self
  real(kind=8), intent(INOUT) :: send_data(:)
  real(kind=8), intent(IN) :: intpl_data(:)

  call jml_ScatterVLocal(self%my_grid%get_comp_id(), &
                         intpl_data, self%tmap_info%recvcounts, self%tmap_info%displs, &
                         send_data, self%tmap_info%my_data_size)
  
end subroutine scatter_interpolated_data_to_send

!=======+=========+=========+=========+=========+=========+=========+=========+

end module jlt_map_class
