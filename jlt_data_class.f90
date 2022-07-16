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
   integer                       :: num_of_layer = 1       !
   integer                       :: grid_intpl_tag         ! grid interpolation tag used in the interpolation subroutine
   integer                       :: exchange_tag           ! MPI data exchange tag
   integer                       :: exchange_array_size    !
   real(kind=8)                  :: fill_value             ! 
   real(kind=8), pointer         :: data1d(:)              ! exchange  data buffer, size = exchange_array_size
   real(kind=8), pointer         :: data1d_tmp(:)          ! averaging data buffer, size = exchange_array_size
   real(kind=8), pointer         :: data2d(:,:)            ! exchange  data buffer, size = (exchange_array_size, num_of_layer)
   real(kind=8), pointer         :: exchange_buffer(:,:)   ! global data buffer of exchange data
   integer(kind=8)               :: put_sec                ! time of put data
 contains
   procedure :: set_my_exchange     ! subroutine (send_comp, send_grid, recv_comp, recv_grid)
   procedure :: get_my_exchange     ! type(exchange_class), pointer function ()
   procedure :: get_my_name         ! character(len=STR_SHORT) function ()
   procedure :: get_send_data_name  ! character(len=STR_SHORT) function ()
   procedure :: get_recv_data_name  ! character(len=STR_SHORT) function ()
   procedure :: is_avr              ! logical function ()
   procedure :: get_intvl           ! integer function ()
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
                         intvl, num_of_layer, intpl_tag, &
                         fill_value, exchange_tag) result(my_data_class)
  implicit none
  character(len=*), intent(IN) :: send_data_name, recv_data_name
  logical, intent(IN)          :: avr_flag
  integer, intent(IN)          :: intvl
  integer, intent(IN)          :: num_of_layer
  integer, intent(IN)          :: intpl_tag
  real(kind=8), intent(IN)     :: fill_value
  integer, intent(IN)          :: exchange_tag
  type(data_class) :: my_data_class

  my_data_class%send_data_name = trim(send_data_name)
  my_data_class%recv_data_name = trim(recv_data_name)
  my_data_class%avr_flag       = avr_flag
  my_data_class%intvl          = intvl
  my_data_class%num_of_layer   = num_of_layer
  my_data_class%grid_intpl_tag = intpl_tag
  my_data_class%fill_value     = fill_value
  my_data_class%exchange_tag   = exchange_tag
  my_data_class%put_sec        = -1
  
end function init_data_class

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine set_my_exchange(self, send_comp_name, send_grid_name, &
                                 recv_comp_name, recv_grid_name)
  use jlt_exchange, only : get_exchange_ptr
  use jlt_map_class, only : map_class
  implicit none
  class(data_class)            :: self
  character(len=*), intent(IN) :: send_comp_name, send_grid_name
  character(len=*), intent(IN) :: recv_comp_name, recv_grid_name
  type(map_class), pointer     :: my_map
  integer :: global_exchange_size
  
  self%my_exchange => get_exchange_ptr(send_comp_name, send_grid_name, &
                                       recv_comp_name, recv_grid_name)

  if (trim(self%my_exchange%get_my_name()) == trim(send_comp_name)) then
     self%my_name = trim(self%send_data_name)
  else
     self%my_name = trim(self%recv_data_name)
  end if

  my_map => self%my_exchange%get_my_map()

  self%exchange_array_size = my_map%get_exchange_array_size()

  if (self%num_of_layer == 1) then
    allocate(self%data1d(self%exchange_array_size))
    self%data1d(:) = 0.d0   
  else
    allocate(self%data2d(self%exchange_array_size, self%num_of_layer))
    self%data2d(:,:) = 0.d0  
  end if
  
  if (self%is_avr()) then
     allocate(self%data1d_tmp(self%exchange_array_size))
  else
     allocate(self%data1d_tmp(1))
  end if

  !if (my_map%get_my_rank() == 0) then ! root only
     
     if (self%num_of_layer == 1) then
        allocate(self%exchange_buffer(self%my_exchange%get_my_exchange_size(),1))
     else
        allocate(self%exchange_buffer(self%my_exchange%get_my_exchange_size(),self%num_of_layer))
     end if
     
  !else
  !   allocate(self%exchange_buffer(1,1))
  !end if
  
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
  use jlt_constant, only : STR_MID
  use jlt_utils, only : put_log, error
  use jlt_map_class, only : map_class
  implicit none
  class(data_class)           :: self
  real(kind=8), intent(IN)    :: data(:)
  integer(kind=8), intent(IN) :: next_sec
  integer(kind=4), intent(IN) :: delta_t
  type(map_class), pointer :: my_map
  real(kind=8) :: weight
  character(len=STR_MID) :: log_str
  integer :: i

  if (self%put_sec == next_sec) then
     call error("put_data_1d", "data ["//trim(self%my_name)//"] double put error") 
  else
     self%put_sec = next_sec
  end if

  if (self%intvl <= 0) then ! every step send
     my_map => self%my_exchange%get_my_map()
     call my_map%local_2_exchange(data, self%data1d)
     call self%my_exchange%send_data_1d(self%data1d, self%exchange_buffer, self%grid_intpl_tag, self%exchange_tag)
     write(log_str,'("  ",A)') "[put_data_1d ] send    the data, data_name = "//trim(self%my_name)
     call put_log(trim(log_str))
  else
     if (self%is_avr()) then ! average data
        my_map => self%my_exchange%get_my_map()
        call my_map%local_2_exchange(data, self%data1d_tmp)

        if (delta_t == 0) then ! first step
           weight = 1.d0
           self%data1d(:) = 0.d0 ! initialize data1d 
        else
           weight = dble(delta_t)/dble(self%intvl)
        end if
       
        do i = 1, self%exchange_array_size
           self%data1d(i) = self%data1d(i) + self%data1d_tmp(i)*weight
        end do

        write(log_str,'("  ",A,F8.5)') "[put_data_1d ] average the data, data_name = "//trim(self%my_name)//", weight = ",weight
        call put_log(trim(log_str))

        if (mod(next_sec, self%intvl) == 0) then ! send step
           call self%my_exchange%send_data_1d(self%data1d, self%exchange_buffer, self%grid_intpl_tag, self%exchange_tag)
           self%data1d(:) = 0.d0 ! reset average data
           write(log_str,'("  ",A)') "[put_data_1d ] send    the data, data_name = "//trim(self%my_name)
           call put_log(trim(log_str))
        end if

     else ! snapshot data
        if (mod(next_sec, self%intvl) == 0) then ! send step        
           my_map => self%my_exchange%get_my_map()
           call my_map%local_2_exchange(data, self%data1d)
           call self%my_exchange%send_data_1d(self%data1d, self%exchange_buffer, self%grid_intpl_tag, self%exchange_tag)
           write(log_str,'("  ",A)') "[put_data_1d ] send    the data, data_name = "//trim(self%my_name)
           call put_log(trim(log_str))
        end if
     end if   

  end if
  
end subroutine put_data_1d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine put_data_2d(self, data, next_sec, delta_t)
  use jlt_constant, only : STR_MID
  use jlt_utils, only : put_log, error
  use jlt_map_class, only : map_class
  implicit none
  class(data_class)           :: self
  real(kind=8), intent(IN)    :: data(:,:)
  integer(kind=8), intent(IN) :: next_sec
  integer(kind=4), intent(IN) :: delta_t
  type(map_class), pointer :: my_map
  real(kind=8) :: weight
  character(len=STR_MID) :: log_str
  integer :: i, k
  
  if (self%put_sec == next_sec) then
     call error("put_data_2d", "data ["//trim(self%my_name)//"] double put error") 
  else
     self%put_sec = next_sec
  end if

  if (self%intvl <= 0) then ! every step send
     my_map => self%my_exchange%get_my_map()

     do k = 1, self%get_num_of_layer()
        call my_map%local_2_exchange(data(:,k), self%data2d(:,k))
     end do
     
     call self%my_exchange%send_data_2d(self%data2d, self%exchange_buffer, self%grid_intpl_tag, self%exchange_tag)

     write(log_str,'("  ",A)') "[put_data_2d ] send    the data, data_name = "//trim(self%my_name)
     call put_log(trim(log_str))
  else
     if (self%is_avr()) then ! average data
        my_map => self%my_exchange%get_my_map()

        if (delta_t == 0) then ! first step
           weight = 1.d0
           self%data2d(:,:) = 0.d0 ! initialize data2d
        else
           weight = dble(delta_t)/dble(self%intvl)
        end if
       

        do k = 1, self%get_num_of_layer()

          call my_map%local_2_exchange(data(:,k), self%data1d_tmp)
          
          do i = 1, self%exchange_array_size
             self%data2d(i,k) = self%data2d(i,k) + self%data1d_tmp(i)*weight
          end do

       end do
       
        write(log_str,'("  ",A,F8.5)') "[put_data_2d ] average the data, data_name = "//trim(self%my_name)//", weight = ",weight
        call put_log(trim(log_str))

        if (mod(next_sec, self%intvl) == 0) then ! send step
           call self%my_exchange%send_data_2d(self%data2d, self%exchange_buffer, self%grid_intpl_tag, self%exchange_tag)
           self%data2d(:,:) = 0.d0 ! reset average data
           write(log_str,'("  ",A)') "[put_data_1d ] send    the data, data_name = "//trim(self%my_name)
           call put_log(trim(log_str))
        end if

     else ! snapshot data
        if (mod(next_sec, self%intvl) == 0) then ! send step        
           my_map => self%my_exchange%get_my_map()

           do k = 1, self%get_num_of_layer()
              call my_map%local_2_exchange(data(:,k), self%data2d(:,k))
           end do
           
           call self%my_exchange%send_data_2d(self%data2d, self%exchange_buffer, self%grid_intpl_tag, self%exchange_tag)
           write(log_str,'("  ",A)') "[put_data_2d ] send    the data, data_name = "//trim(self%my_name)
           call put_log(trim(log_str))
        end if
     end if   

  end if
  
end subroutine put_data_2d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine recv_data_1d(self, current_sec)
  use jlt_constant, only : STR_MID
  use jlt_utils, only : put_log
  implicit none
  class(data_class)           ::self
  integer(kind=8), intent(IN) :: current_sec
  character(len=STR_MID) :: log_str

  if (mod(current_sec, self%intvl) == 0) then
    call self%my_exchange%recv_data_1d(self%data1d, self%exchange_buffer, self%exchange_tag)
    write(log_str,'("  ",A,I5)') "[recv_data_1d] recv    the data, data_name = "//trim(self%my_name) &
                                //", exchange_tag = ", self%exchange_tag
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
  
  if (mod(current_sec, self%intvl) == 0) then
    call self%my_exchange%recv_data_2d(self%data2d, self%exchange_buffer, self%exchange_tag)
    write(log_str,'("  ",A,I5)') "[recv_data_2d] recv    the data, data_name = "//trim(self%my_name) &
                                //", exchange_tag = ", self%exchange_tag
    call put_log(trim(log_str))
  end if
 
end subroutine recv_data_2d

  
!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine interpolate_data_1d(self)
  use jlt_map_class, only : map_class
  implicit none
  class(data_class)        :: self
  type(map_class), pointer :: my_map
  real(kind=8), pointer    :: intpl_data(:,:)
  real(kind=8), pointer    :: recv_data(:,:)
  
  my_map => self%my_exchange%get_my_map()

  if (self%my_exchange%is_my_intpl()) then
    if (my_map%get_my_rank() == 0) then ! root
       allocate(recv_data(my_map%get_data_size(),1))
       allocate(intpl_data(my_map%get_my_whole_data_size(), 1))
    else
       allocate(recv_data(1,1))
       allocate(intpl_data(1,1))
    end if

    call my_map%gather_recv_data_to_intpl(self%exchange_buffer(:,1), intpl_data(:,1))

    if (my_map%get_my_rank() == 0) then
        call self%my_exchange%interpolate_data(intpl_data, recv_data, 1, self%grid_intpl_tag)
    end if
     
    call my_map%scatter_data_from_intpl(self%data1d, recv_data(:,1))

    deallocate(intpl_data)
    deallocate(recv_data)
  else
    self%data1d(:) = self%exchange_buffer(:,1)
  end if
 
end subroutine interpolate_data_1d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine interpolate_data_2d(self)
  use jlt_map_class, only : map_class
  implicit none
  class(data_class)        :: self
  type(map_class), pointer :: my_map
  real(kind=8), pointer    :: intpl_data(:,:)
  real(kind=8), pointer    :: recv_data(:,:)
  integer                  :: num_of_layer
  integer                  :: k

  num_of_layer = size(self%data2d, 2)
  
  my_map => self%my_exchange%get_my_map()

  if (self%my_exchange%is_my_intpl()) then
    if (my_map%get_my_rank() == 0) then ! root
       allocate(recv_data(my_map%get_data_size(),num_of_layer))
       allocate(intpl_data(my_map%get_my_whole_data_size(), num_of_layer))
    else
       allocate(recv_data(1,1))
       allocate(intpl_data(1,1))
    end if

    do k = 1, num_of_layer
      call my_map%gather_recv_data_to_intpl(self%exchange_buffer(:,k), intpl_data(:,k))
   end do
   
    if (my_map%get_my_rank() == 0) then
      call self%my_exchange%interpolate_data(intpl_data, recv_data, num_of_layer, self%grid_intpl_tag)
    end if
     
    do k = 1, num_of_layer
      call my_map%scatter_data_from_intpl(self%data2d(:,k), recv_data(:,k))
    end do
   
    deallocate(intpl_data)
    deallocate(recv_data)
  else
    self%data2d(:,:) = self%exchange_buffer(:,:)
  end if
 
end subroutine interpolate_data_2d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine get_data_1d(self, data, current_sec, is_get_ok)
  use jlt_constant, only : STR_MID
  use jlt_utils, only : put_log
  use jlt_map_class, only : map_class
  implicit none
  class(data_class) :: self
  real(kind=8), intent(INOUT) :: data(:)
  integer(kind=8), intent(IN) :: current_sec
  logical, intent(INOUT)      :: is_get_ok
  type(map_class), pointer :: my_map
  character(len=STR_MID) :: log_str

  my_map => self%my_exchange%get_my_map()

  if (mod(current_sec, self%intvl) == 0) then

     data(:) = self%fill_value  ! set fill value

     write(log_str,'("  ",A,F)') "[get_data_1d ] get data, data_name = "//trim(self%my_name), data(1)
     call put_log(trim(log_str))

     call my_map%exchange_2_local(data, self%data1d)

     is_get_ok = .true.

     write(log_str,'("  ",A,F)') "[get_data_1d ] get data, data_name = "//trim(self%my_name), data(1)
     call put_log(trim(log_str))
  else
     is_get_ok = .false.
   end if
   
end subroutine get_data_1d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine get_data_2d(self, data, current_sec, is_get_ok)
  use jlt_constant, only : STR_MID
  use jlt_utils, only : put_log
  use jlt_map_class, only : map_class
  implicit none
  class(data_class) :: self
  real(kind=8), intent(INOUT) :: data(:,:)
  integer(kind=8), intent(IN) :: current_sec
  logical, intent(INOUT)      :: is_get_ok
  type(map_class), pointer :: my_map
  character(len=STR_MID) :: log_str
  integer :: k
  
  my_map => self%my_exchange%get_my_map()

  if (mod(current_sec, self%intvl) == 0) then

     data(:,:) = self%fill_value  ! set fill value

     do k = 1, self%get_num_of_layer()
        call my_map%exchange_2_local(data(:,k), self%data2d(:,k))
     end do
     
     is_get_ok = .true.

     write(log_str,'("  ",A)') "[get_data_2d ] get data, data_name = "//trim(self%my_name)
     call put_log(trim(log_str))
  else
     is_get_ok = .false.
   end if
   
end subroutine get_data_2d

!=======+=========+=========+=========+=========+=========+=========+=========+

end module jlt_data_class
