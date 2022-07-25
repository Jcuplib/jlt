!
!Copyright (c) 2011, arakawa@rist.jp
!All rights reserved.
!
module jlt_data
  use jlt_constant, only : STR_SHORT, STR_LONG
  use jlt_data_class, only : data_class
  implicit none
  private

!--------------------------------   public   ---------------------------------!

public :: init_data                         ! subroutine (my_comp_name) 
public :: set_send_data                     ! subroutine (send_comp, send_grid, send_data_name, recv_comp, recv_grid, recv_data_name, 
                                            !             is_avr, intvl, num_of_layer, &
                                            !             grid_intpl_tag, fill_value, exchange_tag)
public :: set_recv_data                     ! subroutine (send_comp, send_grid, send_data_name, recv_comp, recv_grid, recv_data_name, 
                                            !             is_avr, intvl, num_of_layer, &
                                            !             grid_intpl_tag, fill_value, exchange_tag)
!public :: set_my_data                       ! subroutine ()
public :: get_num_of_send_data              ! integer function ()
public :: get_num_of_recv_data              ! integer function ()
public :: get_send_data_name                ! character(len=STR_SHORT) function (data_num)
public :: get_recv_data_name                ! character(len=STR_SHORT) function (data_num)
public :: get_send_data_intvl               ! integer function (data_num)
public :: get_recv_data_intvl               ! integer function (data_num)
public :: recv_my_data                      ! subroutine (data_name)
public :: interpolate_recv_data             ! subroutine (data_name)
public :: put_data_1d                       ! subroutine (data_name, data)
public :: get_data_1d                       ! subroutine (data_name, data)
public :: put_data_2d                       ! subroutine (data_name, data)
public :: get_data_2d                       ! subroutine (data_name, data)

!--------------------------------   private  ---------------------------------!

integer, parameter :: MAX_DATA = 40

integer :: num_of_send_data = 0
type(data_class), pointer :: send_data(:)
integer :: num_of_recv_data = 0
type(data_class), pointer :: recv_data(:)

character(len=STR_SHORT) :: my_comp

contains

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine init_data(my_comp_name)
  implicit none
  character(len=*), intent(IN) :: my_comp_name
  
  allocate(send_data(MAX_DATA))
  allocate(recv_data(MAX_DATA))

  my_comp = trim(my_comp_name)

end subroutine init_data

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine set_send_data(send_comp, send_grid, send_data_name, recv_comp, recv_grid, recv_data_name, is_avr, intvl, num_of_layer, &
                         grid_intpl_tag, fill_value, exchange_tag)
  use jlt_utils, only : put_log
  implicit none
  character(len=*), intent(IN) :: send_comp, send_grid, send_data_name
  character(len=*), intent(IN) :: recv_comp, recv_grid, recv_data_name
  logical, intent(IN)          :: is_avr
  integer, intent(IN)          :: intvl
  integer, intent(IN)          :: num_of_layer
  integer, intent(IN)          :: grid_intpl_tag
  real(kind=8), intent(IN)     :: fill_value
  integer, intent(IN)          :: exchange_tag
  character(len=STR_LONG)      :: log_str

  call put_log("jlt_set_send_data", 2)
  call put_log("  send_comp = "//trim(send_comp), 2)
  call put_log("  send_grid = "//trim(send_grid), 2)
  call put_log("  send_data = "//trim(send_data_name), 2)
  call put_log("  recv_comp = "//trim(recv_comp), 2)
  call put_log("  recv_grid = "//trim(recv_grid), 2)
  call put_log("  recv_data = "//trim(recv_data_name), 2)
  
  num_of_send_data = num_of_send_data + 1

  send_data(num_of_send_data) = data_class(send_data_name, recv_data_name, is_avr, intvl, num_of_layer, &
                                           grid_intpl_tag, fill_value, exchange_tag)
  call send_data(num_of_send_data)%set_my_exchange(send_comp, send_grid, &
                                                   recv_comp, recv_grid)
  
end subroutine set_send_data

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine set_recv_data(send_comp, send_grid, send_data_name, recv_comp, recv_grid, recv_data_name, is_avr, intvl, num_of_layer, &
                         grid_intpl_tag, fill_value, exchange_tag)
  use jlt_utils, only : put_log
  implicit none
  character(len=*), intent(IN) :: send_comp, send_grid, send_data_name
  character(len=*), intent(IN) :: recv_comp, recv_grid, recv_data_name
  logical, intent(IN)          :: is_avr
  integer, intent(IN)          :: intvl
  integer, intent(IN)          :: num_of_layer
  integer, intent(IN)          :: grid_intpl_tag
  real(kind=8), intent(IN)     :: fill_value
  integer, intent(IN)          :: exchange_tag
  
  call put_log("jlt_set_recv_data", 2)
  call put_log("  send_comp = "//trim(send_comp), 2)
  call put_log("  send_grid = "//trim(send_grid), 2)
  call put_log("  send_data = "//trim(send_data_name), 2)
  call put_log("  recv_comp = "//trim(recv_comp), 2)
  call put_log("  recv_grid = "//trim(recv_grid), 2)
  call put_log("  recv_data = "//trim(recv_data_name), 2)

  num_of_recv_data = num_of_recv_data + 1

  recv_data(num_of_recv_data) = data_class(send_data_name, recv_data_name, is_avr, intvl, num_of_layer, &
                                           grid_intpl_tag, fill_value, exchange_tag)
  call recv_data(num_of_recv_data)%set_my_exchange(send_comp, send_grid, &
                                                   recv_comp, recv_grid)
  
end subroutine set_recv_data

!=======+=========+=========+=========+=========+=========+=========+=========+

!subroutine set_my_data()
!  use jlt_namelist, only : exchange_data_type, exchange_conf_type, exchange_data_type, &
!                           get_num_of_exchange_conf, get_exchange_conf_ptr_num, &
!                           get_num_of_exchange_data_conf, get_exchange_data_conf_ptr
!  use jlt_comp, only : is_model_running
!  implicit none
!  type(exchange_conf_type), pointer :: exchange_conf_ptr
!  type(exchange_data_type), pointer :: exchange_data_ptr
!  integer :: i, j

!  do i = 1, get_num_of_exchange_conf()
!    exchange_conf_ptr => get_exchange_conf_ptr_num(i)
!    if ((trim(my_comp) == trim(exchange_conf_ptr%send_comp))) then ! set send data config
!       if (is_model_running(trim(exchange_conf_ptr%recv_comp))) then ! if target comp is alive
!         do j =1, get_num_of_exchange_data_conf(i)
!            num_of_send_data = num_of_send_data + 1
!            exchange_data_ptr => get_exchange_data_conf_ptr(i, j)
!            send_data(num_of_send_data) = data_class(exchange_data_ptr%send_data, exchange_data_ptr%recv_data, &
!                                                     exchange_data_ptr%is_avr, &
!                                                     exchange_data_ptr%intvl, exchange_data_ptr%num_of_layer, &
!                                                     exchange_data_ptr%grid_intpl_tag, &
!                                                     exchange_data_ptr%fill_value, &
!                                                     exchange_data_ptr%exchange_tag)
!            call send_data(num_of_send_data)%set_my_exchange(exchange_conf_ptr%send_comp, exchange_conf_ptr%send_grid, &
!                                                             exchange_conf_ptr%recv_comp, exchange_conf_ptr%recv_grid)
          
!         end do
!       end if
!    end if
!  end do
  
!  do i = 1, get_num_of_exchange_conf()
!    exchange_conf_ptr => get_exchange_conf_ptr_num(i)
!    if ((trim(my_comp) == trim(exchange_conf_ptr%recv_comp))) then ! set recv data config
!       if (is_model_running(trim(exchange_conf_ptr%send_comp))) then
!         do j =1, get_num_of_exchange_data_conf(i)
!            num_of_recv_data = num_of_recv_data + 1
!            exchange_data_ptr => get_exchange_data_conf_ptr(i, j)
!            recv_data(num_of_recv_data) = data_class(exchange_data_ptr%send_data, exchange_data_ptr%recv_data, &
!                                                     exchange_data_ptr%is_avr, &
!                                                     exchange_data_ptr%intvl, exchange_data_ptr%num_of_layer, &
!                                                     exchange_data_ptr%grid_intpl_tag, &
!                                                     exchange_data_ptr%fill_value, &
!                                                     exchange_data_ptr%exchange_tag)
!            call recv_data(num_of_recv_data)%set_my_exchange(exchange_conf_ptr%send_comp, exchange_conf_ptr%send_grid, &
!                                                             exchange_conf_ptr%recv_comp, exchange_conf_ptr%recv_grid)
          
!         end do
!       end if
!    end if
!  end do
  
!  call write_data_config()
  
!end subroutine set_my_data

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine write_data_config()
  use jlt_utils, only : put_log
  implicit none
  character(len=STR_LONG) :: log_str
  integer :: i
  
  write(log_str, *) "send data config : num_of_send_data = ", num_of_send_data
  call put_log(trim(log_str))

  do i = 1, num_of_send_data
     write(log_str, *) "   my_name, recv_data_name = ", trim(send_data(i)%get_my_name()) &
                       //", "//trim(send_data(i)%get_recv_data_name())
     call put_log(trim(log_str))
     write(log_str, *)  "   avr_flag = ", send_data(i)%is_avr(), ", intvl = ", send_data(i)%get_intvl() 
     call put_log(trim(log_str))
     write(log_str, *)  "   num_of_layer = ", send_data(i)%get_num_of_layer(), ", tag = ", send_data(i)%get_exchange_tag()
     call put_log(trim(log_str))
  end do

  write(log_str, *) "send data config : num_of_send_data = ", num_of_send_data
  do i = 1, num_of_recv_data
     write(log_str, *) "   my_name, send_data_name = ", trim(recv_data(i)%get_my_name()) &
                       //", "//trim(recv_data(i)%get_send_data_name())
     call put_log(trim(log_str))
     write(log_str, *)  "   avr_flag = ", send_data(i)%is_avr(), ", intvl = ", send_data(i)%get_intvl() 
     call put_log(trim(log_str))
     write(log_str, *)  "   num_of_layer = ", send_data(i)%get_num_of_layer(), ", tag = ", send_data(i)%get_exchange_tag()
     call put_log(trim(log_str))
  end do
  
end subroutine write_data_config

  

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function get_num_of_send_data()
  implicit none

  get_num_of_send_data = num_of_send_data

end function get_num_of_send_data

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function get_num_of_recv_data()
  implicit none

  get_num_of_recv_data = num_of_recv_data

end function get_num_of_recv_data

!=======+=========+=========+=========+=========+=========+=========+=========+

character(len=STR_SHORT) function get_send_data_name(data_num)
  implicit none
  integer, intent(IN) :: data_num

  get_send_data_name = send_data(data_num)%get_my_name()

end function get_send_data_name

!=======+=========+=========+=========+=========+=========+=========+=========+

character(len=STR_SHORT) function get_recv_data_name(data_num)
  implicit none
  integer, intent(IN) :: data_num

  get_recv_data_name = recv_data(data_num)%get_my_name()

end function get_recv_data_name

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function get_send_data_intvl(data_num)
  implicit none
  integer, intent(IN) :: data_num

  get_send_data_intvl = send_data(data_num)%get_intvl()

end function get_send_data_intvl

!=======+=========+=========+=========+=========+=========+=========+=========+

integer function get_recv_data_intvl(data_num)
  implicit none
  integer, intent(IN) :: data_num

  get_recv_data_intvl = recv_data(data_num)%get_intvl()

end function get_recv_data_intvl

!=======+=========+=========+=========+=========+=========+=========+=========+

function get_send_data_ptr(data_name) result(res)
  use jlt_utils, only : error
  
  implicit none
  character(len=*), intent(IN) :: data_name
  type(data_class), pointer :: res
  integer :: i

  do i = 1, num_of_send_data
     if (trim(data_name) == trim(send_data(i)%get_my_name())) then
        res => send_data(i)
        return
     end if
  end do

  call error("[jlt_data:get_send_data_ptr]", "no such data, data name = "//trim(data_name))
  
end function get_send_data_ptr

!=======+=========+=========+=========+=========+=========+=========+=========+

function get_recv_data_ptr(data_name) result(res)
  use jlt_utils, only : error
  
  implicit none
  character(len=*), intent(IN) :: data_name
  type(data_class), pointer :: res
  integer :: i

  do i = 1, num_of_recv_data
     if (trim(data_name) == trim(recv_data(i)%get_my_name())) then
        res => recv_data(i)
        return
     end if
  end do

  call error("[jlt_data:get_recv_data_ptr]", "no such data, data name = "//trim(data_name))
  
end function get_recv_data_ptr

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine put_data_1d(data_name, data, next_sec, delta_t)
  use jlt_utils, only : error
  implicit none
  character(len=*), intent(IN) :: data_name
  real(kind=8), intent(IN)     :: data(:)
  integer(kind=8), intent(IN)  :: next_sec  ! next step seconds
  integer(kind=4), intent(IN)  :: delta_t   ! current step delta_t
  type(data_class), pointer :: data_ptr

  data_ptr => get_send_data_ptr(data_name)

  if (data_ptr%get_num_of_layer() > 1) then
     call error("[jlt_data:put_data_1d]", "data dimension mismatch")
  end if
  
  call data_ptr%put_data_1d(data, next_sec, delta_t)
  
end subroutine put_data_1d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine put_data_2d(data_name, data, next_sec, delta_t)
  use jlt_utils, only : error
  implicit none
  character(len=*), intent(IN) :: data_name
  real(kind=8), intent(IN)     :: data(:,:)
  integer(kind=8), intent(IN)  :: next_sec  ! next step seconds
  integer(kind=4), intent(IN)  :: delta_t   ! current step delta_t
  type(data_class), pointer :: data_ptr

  data_ptr => get_send_data_ptr(data_name)

  if (data_ptr%get_num_of_layer() <= 1) then
     call error("[jlt_data:put_data_2d]", "data dimension mismatch")
  end if
  
  call data_ptr%put_data_2d(data, next_sec, delta_t)
  
end subroutine put_data_2d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine recv_my_data(data_name, current_sec)
  implicit none
  character(len=*), intent(IN) :: data_name
  integer(kind=8),  intent(IN) :: current_sec  ! current step seconds
  type(data_class), pointer :: data_ptr

  data_ptr => get_recv_data_ptr(data_name)

  if (data_ptr%get_num_of_layer() == 1) then
    call data_ptr%recv_data_1d(current_sec)
  else
    call data_ptr%recv_data_2d(current_sec)
  end if
  
end subroutine recv_my_data

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine interpolate_recv_data(data_name)
  implicit none
  character(len=*), intent(IN) :: data_name
  type(data_class), pointer :: data_ptr

  data_ptr => get_recv_data_ptr(data_name)

  if (data_ptr%get_num_of_layer() == 1) then
    call data_ptr%interpolate_data_1d()
  else
    call data_ptr%interpolate_data_2d()
  end if
  
end subroutine interpolate_recv_data

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine get_data_1d(data_name, data, current_sec, is_get_ok)
  use jlt_utils, only : error
  implicit none
  character(len=*), intent(IN) :: data_name
  real(kind=8), intent(INOUT)  :: data(:)
  integer(kind=8), intent(IN)  :: current_sec
  logical, intent(INOUT)       :: is_get_ok
  type(data_class), pointer :: data_ptr

  data_ptr => get_recv_data_ptr(data_name)

  if (data_ptr%get_num_of_layer() > 1) then
     call error("[jlt_data:get_data_1d]", "data dimension mismatch")
  end if
  
  call data_ptr%get_data_1d(data, current_sec, is_get_ok)
  
end subroutine get_data_1d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine get_data_2d(data_name, data, current_sec, is_get_ok)
  use jlt_utils, only : error
  implicit none
  character(len=*), intent(IN) :: data_name
  real(kind=8), intent(INOUT)  :: data(:,:)
  integer(kind=8), intent(IN)  :: current_sec
  logical, intent(INOUT)       :: is_get_ok
  type(data_class), pointer :: data_ptr

  data_ptr => get_recv_data_ptr(data_name)

  if (data_ptr%get_num_of_layer() <= 1) then
     call error("[jlt_data:get_data_2d]", "data dimension mismatch")
  end if
  
  call data_ptr%get_data_2d(data, current_sec, is_get_ok)
  
end subroutine get_data_2d

!=======+=========+=========+=========+=========+=========+=========+=========+

end module jlt_data
