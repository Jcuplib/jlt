module jlt_data_class
  use jlt_constant, only : STR_SHORT
  use jlt_exchange_class, only : exchange_class
  implicit none
  private
  
!--------------------------------   public  ----------------------------------!

public :: data_class

!--------------------------------   private  ---------------------------------!


type data_class
   private
   character(len=STR_SHORT)      :: my_name
   character(len=STR_SHORT)      :: send_data_name
   character(len=STR_SHORT)      :: recv_data_name
   type(exchange_class), pointer :: my_exchange
   logical                       :: avr_flag               ! average data or not
   integer                       :: intvl                  ! exchange interval (SEC)
   integer                       :: time_lag               ! exchange time lag (-1, 1, 0)
   integer                       :: exchange_type          ! CONCURRENT_SEND_RECV, ADVANCE_SEND_RECV, BEHIND_SEND_RECV, IMMEDIATE_SEND_RECV
   integer                       :: num_of_layer = 1       !
   real(kind=8)                  :: offset = 0.d0          ! offset value used in get_data
   real(kind=8)                  :: factor = 1.d0          ! factor value used in get_data
   integer                       :: grid_intpl_tag         ! grid interpolation tag used in the interpolation subroutine
   integer                       :: exchange_tag           ! MPI data exchange tag
   integer                       :: exchange_data_size     ! size of exchange data 
   integer                       :: exchange_buffer_size   ! size of received data before interpolation
   real(kind=8)                  :: fill_value             ! 
   real(kind=8), pointer         :: data1d(:)              ! exchange  data buffer, size = exchange_data_size
   real(kind=8), pointer         :: data1d_tmp(:)          ! averaging data buffer, size = exchange_data_size
   real(kind=8), pointer         :: data2d(:,:)            ! exchange  data buffer, size = (exchange_array_size, num_of_layer)
   real(kind=8), pointer         :: exchange_buffer(:,:)   ! received data buffer of exchange data
   integer(kind=8)               :: put_sec                ! time of put data
 contains
   procedure :: set_my_exchange     ! subroutine (send_comp, send_grid, recv_comp, recv_grid)
   procedure :: get_my_exchange     ! type(exchange_class), pointer function ()
   procedure :: get_my_name         ! character(len=STR_SHORT) function ()
   procedure :: get_send_data_name  ! character(len=STR_SHORT) function ()
   procedure :: get_recv_data_name  ! character(len=STR_SHORT) function ()
   procedure :: is_avr              ! logical function ()
   procedure :: get_intvl           ! integer function ()
   procedure :: get_time_lag        ! integer function ()
   procedure :: get_exchange_type   ! integer function ()
   procedure :: get_num_of_layer    ! integer function ()
   procedure :: get_exchange_tag    ! integer function ()
   procedure :: put_data_1d         ! subroutine (data)
   procedure :: recv_data_1d        ! subroutine ()
   procedure :: interpolate_data_1d ! subroutine ()
   procedure :: get_data_1d         ! subroutine (data)
   procedure :: put_data_2d         ! subroutine (data)
   procedure :: recv_data_2d        ! subroutine ()
   procedure :: interpolate_data_2d ! subroutine ()
   procedure :: get_data_2d         ! subroutine (data)
end type data_class

interface data_class
   module procedure init_data_class
end interface data_class


contains

!=======+=========+=========+=========+=========+=========+=========+=========+

function init_data_class(send_data_name, recv_data_name, avr_flag, &
                         intvl, time_lag, exchange_type, num_of_layer, intpl_tag, &
                         fill_value, exchange_tag, factor, offset) result(my_data_class)
  implicit none
  character(len=*), intent(IN)       :: send_data_name, recv_data_name
  logical, intent(IN)                :: avr_flag
  integer, intent(IN)                :: intvl
  integer, intent(IN)                :: time_lag
  integer, intent(IN)                :: exchange_type
  integer, intent(IN)                :: num_of_layer
  integer, intent(IN)                :: intpl_tag
  real(kind=8), intent(IN)           :: fill_value
  integer, intent(IN)                :: exchange_tag
  real(kind=8), optional, intent(IN) :: factor
  real(kind=8), optional, intent(IN) :: offset
  
  type(data_class) :: my_data_class

  my_data_class%send_data_name = trim(send_data_name)
  my_data_class%recv_data_name = trim(recv_data_name)
  my_data_class%avr_flag       = avr_flag
  my_data_class%intvl          = intvl
  my_data_class%time_lag       = time_lag
  my_data_class%exchange_type  = exchange_type
  my_data_class%num_of_layer   = num_of_layer
  my_data_class%grid_intpl_tag = intpl_tag
  my_data_class%fill_value     = fill_value
  my_data_class%exchange_tag   = exchange_tag
  my_data_class%put_sec        = -1

  if (present(factor)) then
     my_data_class%factor = factor
  end if
  
  if (present(offset)) then
     my_data_class%offset = offset
  end if
  
end function init_data_class

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine set_my_exchange(self, send_comp_name, send_grid_name, &
                                 recv_comp_name, recv_grid_name)
  use jlt_exchange, only : get_exchange_ptr
  implicit none
  class(data_class)            :: self
  character(len=*), intent(IN) :: send_comp_name, send_grid_name
  character(len=*), intent(IN) :: recv_comp_name, recv_grid_name
  integer :: global_exchange_size
  
  self%my_exchange => get_exchange_ptr(send_comp_name, send_grid_name, &
                                       recv_comp_name, recv_grid_name)

  if (trim(self%my_exchange%get_my_name()) == trim(send_comp_name)) then
     self%my_name = trim(self%send_data_name)
  else
     self%my_name = trim(self%recv_data_name)
  end if

  ! send side, exchange data size == exchange buffer size 
  self%exchange_data_size   = self%my_exchange%get_exchange_data_size()
  self%exchange_buffer_size = self%my_exchange%get_exchange_buffer_size()

  if (self%num_of_layer == 1) then
    allocate(self%data1d(self%exchange_data_size))
    self%data1d(:) = 0.d0   
  else
    allocate(self%data2d(self%exchange_data_size, self%num_of_layer))
    self%data2d(:,:) = 0.d0  
  end if
  
  if (self%is_avr()) then
     allocate(self%data1d_tmp(self%exchange_data_size))
  else
     allocate(self%data1d_tmp(1))
  end if
   
  if (self%num_of_layer == 1) then
     allocate(self%exchange_buffer(self%exchange_buffer_size,1))
  else
     allocate(self%exchange_buffer(self%exchange_buffer_size,self%num_of_layer))
  end if
  
end subroutine set_my_exchange

!=======+=========+=========+=========+=========+=========+=========+=========+

function get_my_exchange(self) result(res)
  implicit none
  class(data_class)             :: self
  type(exchange_class), pointer :: res

  res => self%my_exchange

end function get_my_exchange

!=======+=========+=========+=========+=========+=========+=========+=========+

character(len=STR_SHORT) function get_my_name(self)
  implicit none
  class(data_class) ::self

  get_my_name = self%my_name

end function get_my_name

!=======+=========+=========+=========+=========+=========+=========+=========+

character(len=STR_SHORT) function get_send_data_name(self)
  implicit none
  class(data_class) :: self

  get_send_data_name = self%send_data_name

end function get_send_data_name

!=======+=========+=========+=========+=========+=========+=========+=========+

character(len=STR_SHORT) function get_recv_data_name(self)
  implicit none
  class(data_class) :: self

  get_recv_data_name = self%recv_data_name

end function get_recv_data_name

!=======+=========+=========+=========+=========+=========+=========+=========+

logical function is_avr(self)
  implicit none
  class(data_class) :: self
  
  is_avr = self%avr_flag

end function is_avr

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function get_intvl(self)
  implicit none
  class(data_class) ::self

  get_intvl = self%intvl

end function get_intvl

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function get_time_lag(self)
  implicit none
  class(data_class) ::self

  get_time_lag = self%time_lag

end function get_time_lag

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function get_exchange_type(self)
  implicit none
  class(data_class) ::self

  get_exchange_type = self%exchange_type

end function get_exchange_type

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function get_num_of_layer(self)
  implicit none
  class(data_class) ::self

  get_num_of_layer = self%num_of_layer

end function get_num_of_layer

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function get_exchange_tag(self)
  implicit none
  class(data_class) :: self

  get_exchange_tag = self%exchange_tag

end function get_exchange_tag

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine put_data_1d(self, data, next_sec, delta_t)
  use jlt_constant, only : STR_MID, ADVANCE_SEND_RECV
  use jlt_utils, only : put_log, error
  use jlt_mpi_lib, only : jml_send_waitall
  implicit none
  class(data_class)           :: self
  real(kind=8), intent(IN)    :: data(:)
  integer(kind=8), intent(IN) :: next_sec
  integer(kind=4), intent(IN) :: delta_t
  real(kind=8) :: weight
  character(len=STR_MID) :: log_str
  integer :: i

  if (self%time_lag == 0) then ! send data now
     write(log_str,'("  ",A)') "[put_data_1d ] put and send data START , data_name = "//trim(self%my_name)
     call put_log(trim(log_str))
     call self%my_exchange%local_2_exchange(data, self%data1d)
     call self%my_exchange%send_data_1d(self%data1d, self%exchange_buffer, self%grid_intpl_tag, self%exchange_tag)
     call jml_send_waitall()
     write(log_str,'("  ",A)') "[put_data_1d ] put and send data END , data_name = "//trim(self%my_name)
     call put_log(trim(log_str))
     return
  end if
  
  if (self%put_sec == next_sec) then
     call error("put_data_1d", "data ["//trim(self%my_name)//"] double put error") 
  else
     self%put_sec = next_sec
  end if


  write(log_str,'("  ",A)') "[put_data_1d ] put data START , data_name = "//trim(self%my_name)
  call put_log(trim(log_str))

  if (self%intvl <= 0) then ! every step send
     call self%my_exchange%local_2_exchange(data, self%data1d)
     call self%my_exchange%send_data_1d(self%data1d, self%exchange_buffer, self%grid_intpl_tag, self%exchange_tag)
  else
     if (self%is_avr()) then ! average data
        call self%my_exchange%local_2_exchange(data, self%data1d_tmp)

        if (delta_t == 0) then ! first step
           weight = 1.d0
           self%data1d(:) = 0.d0 ! initialize data1d 
        else
           weight = dble(delta_t)/dble(self%intvl)
        end if
       
        do i = 1, self%exchange_data_size
           self%data1d(i) = self%data1d(i) + self%data1d_tmp(i)*weight
        end do

        write(log_str,'("  ",A,F8.5)') "[put_data_1d ] average data, data_name = "//trim(self%my_name)//", weight = ",weight
        call put_log(trim(log_str))

        if (delta_t == 0) then ! first step
           if (self%exchange_type == ADVANCE_SEND_RECV) then
              write(log_str,'("  ",A)') "[put_data_1d ] put data SKIP , data_name = "//trim(self%my_name)
              call put_log(trim(log_str))
              return ! skip first step data send
           end if
        end if
        
        if (mod(next_sec, int(self%intvl, kind=8)) == 0) then ! send step
           call self%my_exchange%send_data_1d(self%data1d, self%exchange_buffer, self%grid_intpl_tag, self%exchange_tag)
           self%data1d(:) = 0.d0 ! reset average data
        end if

     else ! snapshot data

        if (delta_t == 0) then ! first step
           if (self%exchange_type == ADVANCE_SEND_RECV) then
              write(log_str,'("  ",A)') "[put_data_1d ] put data SKIP , data_name = "//trim(self%my_name)
              call put_log(trim(log_str))
              return ! skip first step data send
           end if
        end if
     
        if (mod(next_sec, int(self%intvl, kind=8)) == 0) then ! send step        
           call self%my_exchange%local_2_exchange(data, self%data1d)
           call self%my_exchange%send_data_1d(self%data1d, self%exchange_buffer, self%grid_intpl_tag, self%exchange_tag)
        end if
     end if   

  end if
  
  write(log_str,'("  ",A)') "[put_data_1d ] put data END , data_name = "//trim(self%my_name)
  call put_log(trim(log_str))

end subroutine put_data_1d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine put_data_2d(self, data, next_sec, delta_t)
  use jlt_constant, only : STR_MID, ADVANCE_SEND_RECV
  use jlt_utils, only : put_log, error
  use jlt_mpi_lib, only : jml_send_waitall
  implicit none
  class(data_class)           :: self
  real(kind=8), intent(IN)    :: data(:,:)
  integer(kind=8), intent(IN) :: next_sec
  integer(kind=4), intent(IN) :: delta_t
  real(kind=8) :: weight
  character(len=STR_MID) :: log_str
  integer :: i, k
  
  if (self%time_lag == 0) then ! send data now
     write(log_str,'("  ",A)') "[put_data_2d ] put and send data START , data_name = "//trim(self%my_name)
     call put_log(trim(log_str))
     do k = 1, self%get_num_of_layer()
        call self%my_exchange%local_2_exchange(data(:,k), self%data2d(:,k))
     end do
     call self%my_exchange%send_data_2d(self%data2d, self%exchange_buffer, self%num_of_layer, self%grid_intpl_tag, self%exchange_tag)
     call jml_send_waitall()
     write(log_str,'("  ",A)') "[put_data_2d ] put and send data END , data_name = "//trim(self%my_name)
     call put_log(trim(log_str))
     return
  end if

  if (self%put_sec == next_sec) then
     call error("put_data_2d", "data ["//trim(self%my_name)//"] double put error") 
  else
     self%put_sec = next_sec
  end if

  write(log_str,'("  ",A)') "[put_data_2d ] put data START , data_name = "//trim(self%my_name)
  call put_log(trim(log_str))

  if (self%intvl <= 0) then ! every step send

     do k = 1, self%get_num_of_layer()
        call self%my_exchange%local_2_exchange(data(:,k), self%data2d(:,k))
     end do
     
     call self%my_exchange%send_data_2d(self%data2d, self%exchange_buffer, self%num_of_layer, self%grid_intpl_tag, self%exchange_tag)

  else
     if (self%is_avr()) then ! average data

        if (delta_t == 0) then ! first step
           weight = 1.d0
           self%data2d(:,:) = 0.d0 ! initialize data2d
        else
           weight = dble(delta_t)/dble(self%intvl)
        end if
       
        do k = 1, self%get_num_of_layer()

          do i = 1, self%exchange_data_size
             self%data2d(i,k) = self%data2d(i,k) + self%data1d_tmp(i)*weight
          end do

       end do
       
        write(log_str,'("  ",A,F8.5)') "[put_data_2d ] average data, data_name = "//trim(self%my_name)//", weight = ",weight
        call put_log(trim(log_str))

        if (delta_t == 0) then ! first step
           if (self%exchange_type == ADVANCE_SEND_RECV) then
              write(log_str,'("  ",A)') "[put_data_2d ] put data SKIP , data_name = "//trim(self%my_name)
              call put_log(trim(log_str))
              return ! skip first step data send
           end if
        end if
        
        if (mod(next_sec, int(self%intvl, kind=8)) == 0) then ! send step
           call self%my_exchange%send_data_2d(self%data2d, self%exchange_buffer, self%num_of_layer, self%grid_intpl_tag, self%exchange_tag)
           self%data2d(:,:) = 0.d0 ! reset average data
        end if

     else ! snapshot data

        if (delta_t == 0) then ! first step
           if (self%exchange_type == ADVANCE_SEND_RECV) then
              write(log_str,'("  ",A)') "[put_data_2d ] put data SKIP , data_name = "//trim(self%my_name)
              call put_log(trim(log_str))
              return ! skip first step data send
           end if
        end if
        
        if (mod(next_sec, int(self%intvl, kind=8)) == 0) then ! send step        
           do k = 1, self%get_num_of_layer()
              call self%my_exchange%local_2_exchange(data(:,k), self%data2d(:,k))
           end do
           
           call self%my_exchange%send_data_2d(self%data2d, self%exchange_buffer, self%num_of_layer, self%grid_intpl_tag, self%exchange_tag)
        end if
     end if   

  end if
  
  write(log_str,'("  ",A)') "[put_data_2d ] put data END , data_name = "//trim(self%my_name)
  call put_log(trim(log_str))

end subroutine put_data_2d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine recv_data_1d(self, current_sec)
  use jlt_constant, only : STR_MID
  use jlt_utils, only : put_log
  implicit none
  class(data_class)           ::self
  integer(kind=8), intent(IN) :: current_sec
  character(len=STR_MID) :: log_str

  if (mod(current_sec, int(self%intvl, kind=8)) == 0) then
    write(log_str,'("  ",A,I5)') "[recv_data_1d] recv data START, data_name = "//trim(self%my_name) &
                                //", exchange_tag = ", self%exchange_tag
    call put_log(trim(log_str))
    call self%my_exchange%recv_data_1d(self%exchange_buffer, self%exchange_tag)
    write(log_str,'("  ",A)') "[recv_data_1d] recv data END"
    call put_log(trim(log_str))
  end if

end subroutine recv_data_1d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine recv_data_2d(self, current_sec)
  use jlt_constant, only : STR_MID
  use jlt_utils, only : put_log
  implicit none
  class(data_class)           ::self
  integer(kind=8), intent(IN) :: current_sec
  character(len=STR_MID) :: log_str
  
  if (mod(current_sec, int(self%intvl, kind=8)) == 0) then
    write(log_str,'("  ",A,I5)') "[recv_data_2d] recv data, data_name = "//trim(self%my_name) &
                                //", exchange_tag = ", self%exchange_tag
    call put_log(trim(log_str))
    call self%my_exchange%recv_data_2d(self%exchange_buffer, self%num_of_layer, self%exchange_tag)
    write(log_str,'("  ",A)') "[recv_data_2d] recv data END"
    call put_log(trim(log_str))
  end if
 
end subroutine recv_data_2d

  
!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine interpolate_data_1d(self)
  implicit none
  class(data_class)        :: self
  real(kind=8), pointer    :: intpl_data(:,:)

  if (self%my_exchange%is_my_intpl()) then
     allocate(intpl_data(size(self%data1d), 1))
     call self%my_exchange%interpolate_data(self%exchange_buffer, intpl_data, 1, self%grid_intpl_tag, &
                                            self%factor, self%offset)
     self%data1d(:) = intpl_data(:,1)
     deallocate(intpl_data)
  end if
 
end subroutine interpolate_data_1d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine interpolate_data_2d(self)
  implicit none
  class(data_class)        :: self
  real(kind=8), pointer    :: intpl_data(:,:)
  real(kind=8), pointer    :: recv_data(:,:)
  integer                  :: num_of_layer
  integer                  :: k

  num_of_layer = self%num_of_layer
  

  if (self%my_exchange%is_my_intpl()) then
     allocate(intpl_data(size(self%data2d,1), num_of_layer))
     call self%my_exchange%interpolate_data(self%exchange_buffer, intpl_data, num_of_layer, self%grid_intpl_tag, &
                                            self%factor, self%offset)
     self%data2d(:,:) = intpl_data(:,:)
     deallocate(intpl_data)
  end if

end subroutine interpolate_data_2d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine get_data_1d(self, data, current_sec, is_get_ok)
  use jlt_constant, only : STR_MID
  use jlt_utils, only : put_log
  use jlt_mpi_lib, only : jml_recv_waitall
  implicit none
  class(data_class) :: self
  real(kind=8), intent(INOUT) :: data(:)
  integer(kind=8), intent(IN) :: current_sec
  logical, intent(INOUT)      :: is_get_ok
  character(len=STR_MID) :: log_str

  if (self%time_lag == 0) then ! recv data now
    write(log_str,'("  ",A,I5)') "[get_data_1d] recv and get data START, data_name = "//trim(self%my_name) &
                                //", exchange_tag = ", self%exchange_tag
    call put_log(trim(log_str))
    call self%my_exchange%recv_data_1d(self%exchange_buffer, self%exchange_tag)
    call jml_recv_waitall()

    call self%interpolate_data_1d()
    
    data(:) = self%fill_value  ! set fill value
    call self%my_exchange%exchange_2_local(self%data1d, data)
    is_get_ok = .true.
    write(log_str,'("  ",A)') "[get_data_1d] recv and get data END"
    call put_log(trim(log_str))
  end if
  
  if (mod(current_sec, int(self%intvl, kind=8)) == 0) then

     data(:) = self%fill_value  ! set fill value

     call self%my_exchange%exchange_2_local(self%data1d, data)

     is_get_ok = .true.

     write(log_str,'("  ",A,F15.5,F15.5)') "[get_data_1d ] get data OK, data_name = "//trim(self%my_name), minval(data), maxval(data)
     call put_log(trim(log_str))
  else
     write(log_str,'("  ",A)') "[get_data_1d ] get data SKIP"
     call put_log(trim(log_str))
     is_get_ok = .false.
   end if
   
end subroutine get_data_1d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine get_data_2d(self, data, current_sec, is_get_ok)
  use jlt_constant, only : STR_MID
  use jlt_utils, only : put_log
  use jlt_mpi_lib, only : jml_recv_waitall
  implicit none
  class(data_class) :: self
  real(kind=8), intent(INOUT) :: data(:,:)
  integer(kind=8), intent(IN) :: current_sec
  logical, intent(INOUT)      :: is_get_ok
   character(len=STR_MID) :: log_str
  integer :: k
  
  if (self%time_lag == 0) then ! recv data now
    write(log_str,'("  ",A,I5)') "[get_data_2d] recv and get data START, data_name = "//trim(self%my_name) &
                                //", exchange_tag = ", self%exchange_tag
    call put_log(trim(log_str))
    call self%my_exchange%recv_data_2d(self%exchange_buffer, self%num_of_layer, self%exchange_tag)
    call jml_recv_waitall()

    call self%interpolate_data_2d()
    
    data(:,1:self%get_num_of_layer()) = self%fill_value  ! set fill value

    do k = 1, self%get_num_of_layer()
       call self%my_exchange%exchange_2_local(self%data2d(:,k), data(:,k))
    end do

    is_get_ok = .true.

    write(log_str,'("  ",A)') "[get_data_2d] recv and get data END"
    call put_log(trim(log_str))
  end if

  if (mod(current_sec, int(self%intvl, kind=8)) == 0) then

     data(:,1:self%get_num_of_layer()) = self%fill_value  ! set fill value

     do k = 1, self%get_num_of_layer()
        call self%my_exchange%exchange_2_local(self%data2d(:,k), data(:,k))
     end do
     
     is_get_ok = .true.

     write(log_str,'("  ",A,F15.5,F15.5)') "[get_data_2d ] get data OK, data_name = "//trim(self%my_name), &
                                           minval(data(:,1:self%get_num_of_layer())), maxval(data(:,1:self%get_num_of_layer()))
     call put_log(trim(log_str))
  else
     write(log_str,'("  ",A)') "[get_data_2d ] get data SKIP"
     call put_log(trim(log_str))
     is_get_ok = .false.
   end if
   
end subroutine get_data_2d

!=======+=========+=========+=========+=========+=========+=========+=========+

end module jlt_data_class
