module jlt_exchange
  use jlt_exchange_class, only : exchange_class
  implicit none
  private

!--------------------------------   public   ---------------------------------!

  public :: init_exchange                     ! subroutine ()
  public :: set_mapping_table                 ! subroutine (my_name, send_comp_name, send_grid_name, recv_comp_name, recv_grid_name,
                                              !             map_tag, intpl_flag, send_grid, recv_grid, coef)
  public :: get_num_of_exchange               ! integer function ()
  public :: is_exchange_assigned              ! logical function (send_comp_name, send_grid_name, recv_comp_name, recv_grid_name)
  public :: get_exchange_ptr                  ! type(exchange_type), pointer, function (send_comp_name, send_grid_name,
                                              !                                         recv_comp_name, recv_grid_name, intpl_tag)
  public :: write_exchange                    ! subroutine (fid)
  public :: read_exchange                     ! subroutine (fid)

 !--------------------------------   private  ---------------------------------!

  integer, parameter :: MAX_EXCHANGE = 8
  integer :: num_of_exchange
  
  type(exchange_class), allocatable, target :: exchange(:)

  interface get_exchange_ptr
     module procedure get_exchange_ptr_name, get_exchange_ptr_num
  end interface
  
contains

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine init_exchange()
  implicit none

  num_of_exchange = 0
  allocate(exchange(MAX_EXCHANGE))
  
end subroutine init_exchange

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine set_mapping_table(my_name, send_comp_name, send_grid_name, recv_comp_name, recv_grid_name, &
                             map_tag, intpl_flag, intpl_mode, send_grid, recv_grid, coef)
  use jlt_grid, only : get_grid_ptr
  use jlt_utils, only : error
  
  implicit none
  character(len=*), intent(IN) :: my_name
  character(len=*), intent(IN) :: send_comp_name, send_grid_name
  character(len=*), intent(IN) :: recv_comp_name, recv_grid_name
  integer, intent(IN)          :: map_tag
  logical, intent(IN)          :: intpl_flag
  integer, intent(IN)          :: intpl_mode
  integer, intent(IN)          :: send_grid(:), recv_grid(:)
  real(kind=8), intent(IN)     :: coef(:)

  num_of_exchange = num_of_exchange + 1

  if (num_of_exchange > MAX_EXCHANGE) then
     call error("set_mapping_table", "num_of_exchange exceeded MAX_EXCHANGE. Please increase MAX_EXCHANGE")
  end if
  
  exchange(num_of_exchange) = exchange_class(trim(my_name), intpl_flag, intpl_mode)

  call exchange(num_of_exchange)%set_mapping_table(trim(send_comp_name), trim(send_grid_name), &
                                                   trim(recv_comp_name), trim(recv_grid_name), &
                                                   map_tag, intpl_flag, send_grid, recv_grid, coef)

end subroutine set_mapping_table

!=======+=========+=========+=========+=========+=========+=========+=========+

function get_num_of_exchange() result(res)
  implicit none
  integer :: res

  res = num_of_exchange

end function get_num_of_exchange

!=======+=========+=========+=========+=========+=========+=========+=========+

function is_exchange_assigned(send_comp_name, send_grid_name, &
                              recv_comp_name, recv_grid_name, map_tag) result (res)
  implicit none
  character(len=*), intent(IN) :: send_comp_name, send_grid_name
  character(len=*), intent(IN) :: recv_comp_name, recv_grid_name
  integer, intent(IN)          :: map_tag
  logical :: res
  integer :: i

  do i = 1, num_of_exchange
     if (exchange(i)%is_my_exchange(send_comp_name, send_grid_name, &
                                    recv_comp_name, recv_grid_name, map_tag)) then
        res = .true.
        return
     end if
  end do

  res = .false.
  
end function is_exchange_assigned

!=======+=========+=========+=========+=========+=========+=========+=========+

function get_exchange_ptr_name(send_comp_name, send_grid_name, &
                          recv_comp_name, recv_grid_name, map_tag) result (res)
  use jlt_utils, only : error
  implicit none
  character(len=*), intent(IN) :: send_comp_name, send_grid_name
  character(len=*), intent(IN) :: recv_comp_name, recv_grid_name
  integer, intent(IN)          :: map_tag
  type(exchange_class), pointer :: res
  integer :: i

  do i = 1, num_of_exchange
     if (exchange(i)%is_my_exchange(send_comp_name, send_grid_name, &
                                    recv_comp_name, recv_grid_name, map_tag)) then
        res => exchange(i)
        return
     end if
  end do

  call error("[jlt_exchange:get_exchange_ptr]","no such name:"//trim(send_comp_name)//":"// &
             trim(send_grid_name)//":"//trim(recv_comp_name)//":"//trim(recv_grid_name))
  
  res => null()
  
end function get_exchange_ptr_name

!=======+=========+=========+=========+=========+=========+=========+=========+

function get_exchange_ptr_num(exchange_num) result (res)
  use jlt_utils, only : error
  implicit none
  integer, intent(IN) :: exchange_num
  type(exchange_class), pointer :: res

  if (exchange_num > num_of_exchange) then
     call error("[jlt_exchange:get_exchange_ptr]","exchange num exceeded the defiend exchange pattern")
  end if
  
  res => exchange(exchange_num)
  
end function get_exchange_ptr_num

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine write_exchange(fid)
  implicit none
  integer, intent(IN) :: fid
  integer :: i

  write(fid) num_of_exchange
  
  do i = 1, num_of_exchange
     call exchange(i)%write_exchange_class(fid)
  end do

end subroutine write_exchange

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine read_exchange(fid)
  implicit none
  integer, intent(IN) :: fid
  integer :: i

  read(fid) num_of_exchange

  do i = 1, num_of_exchange
     call exchange(i)%read_exchange_class(fid)
  end do

end subroutine read_exchange

!=======+=========+=========+=========+=========+=========+=========+=========+

end module jlt_exchange
