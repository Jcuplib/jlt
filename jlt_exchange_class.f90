module jlt_exchange_class
  use jlt_constant, only   : STR_SHORT, STR_MID, STR_LONG
  use jlt_grid_class, only : grid_class
  use jlt_map_class, only  : map_class
  implicit none
  private

!--------------------------------   public   ---------------------------------!

  public :: exchange_class

!--------------------------------   private  ---------------------------------!

  type exchange_class
     private
     character(len=STR_SHORT)      :: my_name
     character(len=STR_SHORT)      :: send_comp_name
     character(len=STR_SHORT)      :: send_grid_name
     character(len=STR_SHORT)      :: recv_comp_name
     character(len=STR_SHORT)      :: recv_grid_name
     integer                       :: map_tag
     integer, allocatable          :: send_grid_index(:)
     integer, allocatable          :: recv_grid_index(:)
     real(kind=8), allocatable     :: coef(:)
     
     type(map_class), pointer      :: my_map
     type(map_class), pointer      :: target_map

     logical                       :: intpl_flag = .false.  ! my interpolation or not

     integer                       :: my_exchange_data_size ! 
     
   contains
     procedure :: set_mapping_table       ! subroutine (send_comp, send_grid, recv_comp, recv_grid, map_tag, send_index, recv_index, coef)
     procedure :: is_my_exchange          ! logicail function (send_comp, send_grid, recv_comp, recv_grid, map_tag)
     procedure :: get_my_name             ! character(len=STR_SHORT) function()
     procedure :: get_my_map              ! type(map_class), pointer  function()
     procedure :: get_my_exchange_size    ! integer function ()
     procedure :: is_my_intpl             ! logical function ()
     procedure :: send_data_1d            ! subroutine (data, exchange_tag)
     procedure :: recv_data_1d            ! subroutine (data, exchange_tag)
     procedure :: send_data_2d            ! subroutine (data, exchange_tag)
     procedure :: recv_data_2d            ! subroutine (data, exchange_tag)
     procedure :: interpolate_data        ! subroutine (send_data, recv_data, num_of_layer)
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
                             map_tag, send_grid, recv_grid, coef)
  use jlt_grid, only    : get_grid_ptr
  use jlt_comp, only    : get_comp_id_from_name
  use jlt_mpi_lib, only : jml_BcastLocal, jml_GetComm
  use jlt_utils, only   : put_log
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

  call put_log("------------------------------------------------------------------------------------------")
  call put_log("jlt_exchange_class : set_mapping_table start")
  
  allocate(self%send_grid_index(size(send_grid)))
  allocate(self%recv_grid_index(size(recv_grid)))
  allocate(self%coef(size(coef)))

  self%send_comp_name     = trim(send_comp_name)
  self%send_grid_name     = trim(send_grid_name)
  self%recv_comp_name     = trim(recv_comp_name)
  self%recv_grid_name     = trim(recv_grid_name)
  self%map_tag            = map_tag
  self%send_grid_index(:) = send_grid(:)
  self%recv_grid_index(:) = recv_grid(:)
  self%coef(:)            = coef(:)

  if (trim(self%my_name) == trim(send_comp_name)) then
     source_comp_id = get_comp_id_from_name(trim(send_comp_name))
     dest_comp_id   = get_comp_id_from_name(trim(recv_comp_name))
  else
     source_comp_id = get_comp_id_from_name(trim(recv_comp_name))
     dest_comp_id   = get_comp_id_from_name(trim(send_comp_name))
  end if

  allocate(self%my_map)
  allocate(self%target_map)
  
  self%my_map     = map_class(source_comp_id)
  self%target_map = map_class(dest_comp_id)  

  
  if (trim(self%my_name) == trim(send_comp_name)) then
     call self%my_map%set_map_grid(get_grid_ptr(trim(send_grid_name)))
     call self%my_map%set_send_flag()
     call self%target_map%set_recv_flag()

     int_array(1) = size(send_grid)
     call jml_BcastLocal(source_comp_id, int_array, 1, 1, 0)

     write(log_str, *) " set_mapping_table : size of map = ", int_array(1)
     call put_log(trim(log_str))

     allocate(index_array(int_array(1)))
     if (self%my_map%get_my_rank() == 0) index_array(:) = send_grid(:)
     call jml_BcastLocal(source_comp_id, index_array, 1, int_array(1), 0)     

     write(log_str, *) " set_mapping_table : send_grid range = ", minval(index_array), maxval(index_array)
     call put_log(trim(log_str))

     call self%my_map%cal_my_exchange_index(index_array)
     call self%my_map%cal_intpl_index(index_array)
     deallocate(index_array)
  else
     call self%my_map%set_map_grid(get_grid_ptr(trim(recv_grid_name)))
     call self%my_map%set_recv_flag()
     call self%target_map%set_send_flag()

     int_array(1) = size(recv_grid)
     call jml_BcastLocal(source_comp_id, int_array, 1, 1, 0)

     write(log_str, *) " set_mapping_table : size of map = ", int_array(1)
     call put_log(trim(log_str))

     allocate(index_array(int_array(1)))
     if (self%my_map%get_my_rank() == 0) index_array(:) = recv_grid(:)
     call jml_BcastLocal(source_comp_id, index_array, 1, int_array(1), 0)

     write(log_str, *) " set_mapping_table : recv_grid range = ", minval(index_array), maxval(index_array)
     call put_log(trim(log_str))

     call self%my_map%cal_my_exchange_index(index_array)
     call self%my_map%cal_intpl_index(index_array)
     deallocate(index_array)
  end if

  call put_log(" set_mapping_table : exchange map info start")

  if (self%intpl_flag) then ! exchange map info
     call self%target_map%recv_intpl_info(self%my_map%get_comp_id(), self%target_map%get_comp_id())
     call self%my_map%send_intpl_info(self%my_map%get_comp_id(), self%target_map%get_comp_id())
  else
     call self%my_map%send_intpl_info(self%my_map%get_comp_id(), self%target_map%get_comp_id())
     call self%target_map%recv_intpl_info(self%my_map%get_comp_id(), self%target_map%get_comp_id())
  end if

  call self%my_map%cal_local_to_local_exchange_info(self%target_map, self%intpl_flag)

  !if (self%my_map%get_my_rank() == 0) then
     if (self%intpl_flag) then  ! my interpolation, send data size = target_map size
        self%my_exchange_data_size = self%my_map%get_my_target_data_size() !self%target_map%get_data_size()
     else
        self%my_exchange_data_size = self%my_map%get_my_target_data_size() !self%my_map%get_data_size()
     end if
  !else
  !   self%my_exchange_data_size = 0
  !end if
  

  call put_log(" set_mapping_table : exchange map info end")

  call put_log("jlt_exchange_class : set_mapping_table end")

  call put_log("------------------------------------------------------------------------------------------")

end subroutine set_mapping_table

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

function get_my_map(self) result(res)
  implicit none
  class(exchange_class) :: self
  type(map_class), pointer       :: res
  
  res => self%my_map

end function get_my_map

!=======+=========+=========+=========+=========+=========+=========+=========+

function get_my_exchange_size(self) result (res)
  implicit none
  class(exchange_class) :: self
  integer :: res

  res = self%my_exchange_data_size

end function get_my_exchange_size

!=======+=========+=========+=========+=========+=========+=========+=========+

function is_my_intpl(self) result (res)
  implicit none
  class(exchange_class) :: self
  logical :: res

  res = self%intpl_flag

end function is_my_intpl

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine send_data_1d(self, data, exchange_buffer, intpl_tag, exchange_tag)
  use jlt_constant, only : STR_MID
  use jlt_utils, only    : put_log
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
  
  call put_log("------------------------------------------------------------------------------------------")
  
  if ((self%intpl_flag).and.(self%my_map%get_my_rank() == 0)) then
    allocate(send_data(self%my_map%get_data_size(), 1))
    allocate(intpl_data(self%target_map%get_data_size(), 1))
  else
    allocate(send_data(1,1))
    allocate(intpl_data(1,1))   
  end if
  
  if (self%intpl_flag) then ! my interpolation

    call self%my_map%gather_data_to_intpl(data, send_data(:,1))
      
    if (self%my_map%get_my_rank() == 0) then
      call interpolate_data(self, send_data, intpl_data, 1, intpl_tag)
    end if

    call self%my_map%scatter_interpolated_data_to_send(exchange_buffer(:,1), intpl_data(:,1))
        

    do i = 1, self%my_map%get_num_of_target_pe()
       target_rank = self%my_map%get_target_rank(i)
       num_of_data = self%my_map%get_num_of_data(i)
       offset      = self%my_map%get_offset(i)
       data_ptr    => exchange_buffer(offset+1, 1)
       call jml_ISendModel2(self%my_map%get_comp_id(), data_ptr, 1, num_of_data, 1, 1, &
                             self%target_map%get_comp_id(), target_rank, exchange_tag)
    end do
    !call jml_send_waitall()
    write(log_str,'("  ",A,I5)') "[send_data_1d] interpolate and send data OK, exchange_tag = ", exchange_tag
    !call put_log(trim(log_str))
  else
    exchange_buffer(:,1) = data(:)
    do i = 1, self%my_map%get_num_of_target_pe()
       target_rank = self%my_map%get_target_rank(i)
       num_of_data = self%my_map%get_num_of_data(i)
       offset      = self%my_map%get_offset(i)
       data_ptr    => exchange_buffer(offset+1,1)
       call jml_ISendModel2(self%my_map%get_comp_id(), data_ptr, 1, num_of_data, 1, 1, &
                               self%target_map%get_comp_id(), target_rank, exchange_tag)
    end do
    
    !call jml_send_waitall()
    write(log_str,'("  ",A,I5)') "[send_data_1d] send data OK, exchange_tag = ", exchange_tag
    !call put_log(trim(log_str))
 end if
 
 deallocate(send_data)
 deallocate(intpl_data)
 
 call put_log("------------------------------------------------------------------------------------------")
  
end subroutine send_data_1d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine send_data_2d(self, data, exchange_buffer, intpl_tag, exchange_tag)
  use jlt_constant, only : STR_MID
  use jlt_utils, only    : put_log
  use jlt_mpi_lib, only : jml_IsendModel2, jml_send_waitall
  implicit none
  class(exchange_class)       :: self
  real(kind=8), intent(IN)    :: data(:,:)
  real(kind=8), pointer       :: exchange_buffer(:,:)
  integer, intent(IN)         :: intpl_tag
  integer, intent(IN)         :: exchange_tag
  real(kind=8), pointer       :: send_data(:,:)
  real(kind=8), pointer       :: intpl_data(:,:)
  integer                     :: num_of_layer
  character(len=STR_MID)      :: log_str
  integer                     :: target_rank
  integer                     :: num_of_data
  integer                     :: offset
  real(kind=8), pointer       :: data_ptr
  integer :: i, k
  
  call put_log("------------------------------------------------------------------------------------------")
  num_of_layer = size(data, 2)
  
  if ((self%intpl_flag).and.(self%my_map%get_my_rank() == 0)) then
    allocate(send_data(self%my_map%get_data_size(), num_of_layer))
    allocate(intpl_data(self%target_map%get_data_size(), num_of_layer))
  else
    allocate(send_data(1,1))
    allocate(intpl_data(1,1))   
  end if
  
  if (self%intpl_flag) then ! my interpolation

     do k = 1, num_of_layer
       call self%my_map%gather_data_to_intpl(data(:,k), send_data(:,k))
     end do
     
    if (self%my_map%get_my_rank() == 0) then
      call interpolate_data(self, send_data, intpl_data, num_of_layer, intpl_tag)
    end if

    do k = 1, num_of_layer
      call self%my_map%scatter_interpolated_data_to_send(exchange_buffer(:,k), intpl_data(:,k))
        

      do i = 1, self%my_map%get_num_of_target_pe()
         target_rank = self%my_map%get_target_rank(i)
         num_of_data = self%my_map%get_num_of_data(i)
         offset      = self%my_map%get_offset(i)
         data_ptr    => exchange_buffer(offset+1, k)
         call jml_ISendModel2(self%my_map%get_comp_id(), data_ptr, 1, num_of_data, 1, 1, &
                               self%target_map%get_comp_id(), target_rank, exchange_tag)
      end do
    end do
      !call jml_send_waitall()
    write(log_str,'("  ",A,I5)') "[send_data_1d] interpolate and send data OK, exchange_tag = ", exchange_tag
    !call put_log(trim(log_str))
  else
    do k = 1, num_of_layer
      exchange_buffer(:,k) = data(:,k)
      do i = 1, self%my_map%get_num_of_target_pe()
         target_rank = self%my_map%get_target_rank(i)
         num_of_data = self%my_map%get_num_of_data(i)
         offset      = self%my_map%get_offset(i)
         data_ptr    => exchange_buffer(offset+1,k)
         call jml_ISendModel2(self%my_map%get_comp_id(), data_ptr, 1, num_of_data, 1, 1, &
                                 self%target_map%get_comp_id(), target_rank, exchange_tag)
      end do
    end do
    !call jml_send_waitall()
    write(log_str,'("  ",A,I5)') "[send_data_1d] send data OK, exchange_tag = ", exchange_tag
    !call put_log(trim(log_str))
 end if
 
 deallocate(send_data)
 deallocate(intpl_data)
 
 call put_log("------------------------------------------------------------------------------------------")
  
end subroutine send_data_2d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine recv_data_1d(self, data, exchange_buffer, exchange_tag)
  use jlt_constant, only : STR_MID
  use jlt_utils, only    : put_log
  use jlt_mpi_lib, only : jml_IrecvModel2, jml_recv_waitall
  implicit none
  class(exchange_class)       :: self
  real(kind=8), pointer       :: data(:)
  real(kind=8), pointer       :: exchange_buffer(:,:)
  integer, intent(IN)         :: exchange_tag
  real(kind=8), pointer       :: recv_data(:,:)
  real(kind=8), pointer       :: intpl_data(:,:)
  character(len=STR_MID)      :: log_str
  integer                     :: target_rank
  integer                     :: num_of_data
  integer                     :: offset
  real(kind=8), pointer       :: data_ptr
  integer :: i
  
  call put_log("------------------------------------------------------------------------------------------")

  if (self%my_map%get_my_rank() == 0) then ! root
     allocate(recv_data(self%my_map%get_data_size(),1))
     allocate(intpl_data(self%my_map%get_my_whole_data_size(), 1))
  else
     allocate(recv_data(1,1))
     allocate(intpl_data(1,1))
  end if

  if (self%intpl_flag) then ! my interpolation
       
    do i = 1, self%my_map%get_num_of_target_pe()
       target_rank = self%my_map%get_target_rank(i)
       num_of_data = self%my_map%get_num_of_data(i)
       offset      = self%my_map%get_offset(i)
       data_ptr    => exchange_buffer(offset+1, 1)
       call jml_IrecvModel2(self%my_map%get_comp_id(), data_ptr, 1, num_of_data, 1, 1, &
                           self%target_map%get_comp_id(), target_rank, exchange_tag)
     end do

     !call jml_recv_waitall()

     !call self%my_map%gather_recv_data_to_intpl(exchange_buffer(:,1), intpl_data(:,1))

     !if (self%my_map%get_my_rank() == 0) then
     !   call interpolate_data(self, intpl_data, recv_data, 1)
     !end if
     
     !call self%my_map%scatter_data_from_intpl(data, recv_data(:,1))

  else
     do i = 1, self%my_map%get_num_of_target_pe()
       target_rank = self%my_map%get_target_rank(i)
       num_of_data = self%my_map%get_num_of_data(i)
       offset      = self%my_map%get_offset(i)
       data_ptr    => exchange_buffer(offset+1, 1)
       call jml_IrecvModel2(self%my_map%get_comp_id(), data_ptr, 1, num_of_data, 1, 1, &
                           self%target_map%get_comp_id(), target_rank, exchange_tag)
     end do
     !call jml_recv_waitall()
     !data(:) = exchange_buffer(:,1)
  end if
  
  deallocate(recv_data)
  deallocate(intpl_data)
  
  call put_log("------------------------------------------------------------------------------------------")

end subroutine recv_data_1d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine recv_data_2d(self, data, exchange_buffer, exchange_tag)
  use jlt_constant, only : STR_MID
  use jlt_utils, only    : put_log
  use jlt_mpi_lib, only : jml_IrecvModel2, jml_recv_waitall
  implicit none
  class(exchange_class)       :: self
  real(kind=8), pointer       :: data(:,:)
  real(kind=8), pointer       :: exchange_buffer(:,:)
  integer, intent(IN)         :: exchange_tag
  real(kind=8), pointer       :: recv_data(:,:)
  real(kind=8), pointer       :: intpl_data(:,:)
  integer                     :: num_of_layer
  character(len=STR_MID)      :: log_str
  integer                     :: target_rank
  integer                     :: num_of_data
  integer                     :: offset
  real(kind=8), pointer       :: data_ptr
  integer :: i, k
  
  call put_log("------------------------------------------------------------------------------------------")

  num_of_layer = size(data,2)

  if (self%my_map%get_my_rank() == 0) then ! root
     allocate(recv_data(self%my_map%get_data_size(),num_of_layer))
     allocate(intpl_data(self%my_map%get_my_whole_data_size(), num_of_layer))
  else
     allocate(recv_data(1,1))
     allocate(intpl_data(1,1))
  end if

  if (self%intpl_flag) then ! my interpolation
     do k = 1, num_of_layer   
       do i = 1, self%my_map%get_num_of_target_pe()
         target_rank = self%my_map%get_target_rank(i)
         num_of_data = self%my_map%get_num_of_data(i)
         offset      = self%my_map%get_offset(i)
         data_ptr    => exchange_buffer(offset+1, k)
         call jml_IrecvModel2(self%my_map%get_comp_id(), data_ptr, 1, num_of_data, 1, 1, &
                              self%target_map%get_comp_id(), target_rank, exchange_tag)
       end do
     end do
   
     call jml_recv_waitall()

     do k = 1, num_of_layer
       call self%my_map%gather_recv_data_to_intpl(exchange_buffer(:,k), intpl_data(:,k))
     end do
    
     if (self%my_map%get_my_rank() == 0) then
        !call interpolate_data(self, intpl_data, recv_data, num_of_layer)
     end if
     
     !write(log_str,'("  ",A,I5)') "[recv_data_1d] recv and interpolate data OK, exchange_tag = ", exchange_tag
     !call put_log(trim(log_str))
     do k = 1, num_of_layer
       call self%my_map%scatter_data_from_intpl(data(:,k), recv_data(:,k))
     end do

  else
     do k = 1, num_of_layer
       do i = 1, self%my_map%get_num_of_target_pe()
         target_rank = self%my_map%get_target_rank(i)
         num_of_data = self%my_map%get_num_of_data(i)
         offset      = self%my_map%get_offset(i)
         data_ptr    => exchange_buffer(offset+1, k)
         call jml_IrecvModel2(self%my_map%get_comp_id(), data_ptr, 1, num_of_data, 1, 1, &
                             self%target_map%get_comp_id(), target_rank, exchange_tag)
       end do
     end do
     call jml_recv_waitall()
     data(:,:) = exchange_buffer(:,:)
     
     !write(log_str,'("  ",A,I5)') "[recv_data_1d] recv data OK, exchange_tag = ", exchange_tag
     !call put_log(trim(log_str))
  end if
  
  deallocate(recv_data)
  deallocate(intpl_data)
  
  call put_log("------------------------------------------------------------------------------------------")

end subroutine recv_data_2d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine interpolate_data(self, send_data, recv_data, num_of_layer, intpl_tag)
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

  recv_data(:,:) = 0.d0
  
  if (self%my_map%is_send_map()) then
     send_index => self%my_map%get_intpl_index()
     recv_index => self%target_map%get_intpl_index()
  else
     send_index => self%target_map%get_intpl_index()
     recv_index => self%my_map%get_intpl_index()
  end if

  !write(0, *) "interpolate_data ", num_of_layer, self%my_map%get_intpl_size(), size(recv_index), size(send_index)
  !write(0, *) "interpolate_data ", size(recv_data, 1), size(send_data,1), minval(send_index), maxval(send_index),&
  !                                 minval(recv_index), maxval(recv_index), minval(self%coef), maxval(self%coef)
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
        do i = 1, self%my_map%get_intpl_size()
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
    do k = 1, num_of_layer
      do i = 1, self%my_map%get_intpl_size()
         recv_data(recv_index(i),k) = recv_data(recv_index(i),k) + send_data(send_index(i),k)*self%coef(i)
      end do
    end do
 end if
 
  !write(0, *) "interpolate_data min, max = ", minval(send_data), maxval(send_data), minval(recv_data), maxval(recv_data)
end subroutine interpolate_data

!=======+=========+=========+=========+=========+=========+=========+=========+
end module jlt_exchange_class
