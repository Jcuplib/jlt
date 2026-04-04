module jlt_interface
  use jlt_constant, only : NUM_OF_EXCHANGE_DATA, STR_SHORT, STR_MID
  use jlt_mpi_lib, only : jlt_set_world => jml_set_global_comm, &
                          jlt_get_world => jml_get_global_comm
  use jlt_comp, only : jlt_get_num_of_component   => get_num_of_total_component, &
                       jlt_get_component_name     => get_component_name, &
                       jlt_get_comp_num_from_name => get_comp_id_from_name, &
                       jlt_is_my_component        => is_my_component, &
                       jlt_is_model_running       => is_model_running
  use jlt_data, only : jlt_set_data               => set_data
  implicit none
  private
  
!--------------------------------   public  ----------------------------------!
  public :: jlt_varp_type
  public :: jlt_varg_type
  
  public :: jlt_set_world                   ! subroutine (global_comm)
  public :: jlt_get_world                   ! integer function ()
  public :: jlt_set_new_comp                ! subroutine (comp_name)
  public :: jlt_initialize                  ! subroutine (comp_name, time_unit, log_level, log_strerr)
  public :: jlt_coupling_end                ! subroutine (time_array, is_call_mpi_finalize)
  public :: jlt_get_mpi_parameter           ! subroutine (comp_name, comm, group, size, rank)
  public :: jlt_def_grid                    ! subroutine (grid_index, comp_name, grid_name, num_of_vlayer)  
  public :: jlt_end_grid_def                ! subroutine () ! dummy
  public :: jlt_def_varp                    ! subroutine ! dummy subroutine
  public :: jlt_def_varg                    ! subroutine ! dummy subroutine
  public :: jlt_end_var_def                 ! subroutnne () ! dummy subroutine
  public :: jlt_set_data                    ! subroutine (send_comp, send_grid, send_data_name, recv_comp, recv_grid, recv_data_name, map_tag, 
                                            !             is_avr, intvl, time_lag, num_of_layer, &
                                            !             grid_intpl_tag, fill_value, exchange_tag)
  !public :: jlt_set_recv_data               ! subroutine (send_comp, send_grid, send_data_name, recv_comp, recv_grid, recv_data_name, 
  !                                          !             is_avr, intvl, num_of_layer, &
  !                                          !             grid_intpl_tag, fill_value, exchange_tag)
  public :: jlt_set_fill_value              ! subroutine (fill_value)  ! dummy 
  public :: jlt_get_fill_value              ! real(kind=8) function () ! dummy
  public :: jlt_set_mapping_table           ! subroutine (my_name, send_com, send_grid, recv_comp, recv_grid, map_tag,
                                            !             is_recv_intpl, intpl_mode, send_index, recv_index, coef)
  public :: jlt_init_time                   ! subroutine (time_array)
  public :: jlt_set_time                    ! subroutine (current_time, delta_t)
  public :: jlt_put_data                    ! subroutine (data_name, data)
  public :: jlt_get_data                    ! subroutine (data_name, data, is_get_ok)

  public :: jlt_error                       ! subroutine (sub_name, err_str)
  public :: jlt_get_num_of_component        ! integer function ()
  public :: jlt_get_myrank_global           ! integer function ()
  public :: jlt_get_component_name          ! character(len=STR_LEN) function (comp_id)
  public :: jlt_get_comp_num_from_name      ! integer function (comp_name)
  public :: jlt_is_my_component             ! logical function (comp_id) or (comp_name)
  public :: jlt_is_model_running            ! logical function (comp_name)
  public :: jlt_bcast_local                 ! subroutine (my_comp, data)
  public :: jlt_send_array                  ! subroutine (my_comp, recv_comp, send_array)
  public :: jlt_recv_array                  ! subroutine (my_comp, send_comp, recv_array)
  public :: jlt_log                         ! subroutine (sub_name, log_str, log_level)  
  public :: jlt_inc_calendar                ! subroutine (date, delta_t)
  public :: jlt_inc_time                    ! subroutine (comp_name, time_array)
  
  public :: jlt_write_mapping_table         ! subroutine (fid) ! dummy
  public :: jlt_read_mapping_table          ! subroutine (fid) ! dummy

  public :: jlt_set_user_interpolation      ! subroutine (send_comp, send_grid, recv_comp, recv_grid, map_tag, user_func)
  
!--------------------------------   private  ---------------------------------!

  character(len=STR_SHORT) :: my_comp_name = ""
  integer :: my_comp_id                    = -1
  integer :: num_of_comp                   = -1
  integer :: max_num_of_exchange_data      = -1
  
  integer :: my_comm                       = -1
  integer :: my_group                      = -1
  integer :: my_size                       = -1
  integer :: my_rank                       = -1

  interface jlt_put_data
     module procedure jlt_put_data_1d, jlt_put_data_25d
     module procedure jlt_put_data_1d_Ptr, jlt_put_data_25d_ptr
  end interface jlt_put_data

  interface jlt_get_data
     module procedure jlt_get_data_1d, jlt_get_data_25d
     module procedure jlt_get_data_1d_ptr, jlt_get_data_25d_ptr
  end interface jlt_get_data

  interface jlt_bcast_local
     module procedure bcast_array_local_int
  end interface jlt_bcast_local
  
  interface jlt_send_array
    module procedure send_array_to_recv_model_str
    module procedure send_array_to_recv_model_int
    module procedure send_array_to_recv_model_real
    module procedure send_array_to_recv_model_dbl
  end interface

  interface jlt_recv_array
    module procedure recv_array_from_send_model_str
    module procedure jlt_put_data_1d, jlt_put_data_25d
    module procedure recv_array_from_send_model_int
    module procedure recv_array_from_send_model_real
    module procedure recv_array_from_send_model_dbl
  end interface


  integer(kind=8) :: current_sec = 0
  integer(kind=8) :: next_sec  = 0
  integer(kind=4) :: current_delta_t = 0

  type jlt_varp_type
    character(STR_SHORT) :: comp_name
    character(STR_SHORT) :: grid_name
    character(STR_SHORT) :: data_name    
  end type jlt_varp_type

  type jlt_varg_type
    character(STR_SHORT) :: comp_name
    character(STR_SHORT) :: grid_name
    character(STR_SHORT) :: data_name    
  end type jlt_varg_type
  
contains

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine jlt_set_new_comp(component_name)
  use jlt_comp, only : set_my_component
  implicit none
  character(len=*), intent(IN) :: component_name

  my_comp_name = trim(component_name)
  call set_my_component(component_name)

end subroutine jlt_set_new_comp

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine jlt_initialize(model_name, default_time_unit, log_level, log_stderr)
  use jlt_comp, only : init_model_process, get_num_of_total_component, &
                       is_my_component, get_component_name, get_comp_id_from_name
  use jlt_utils, only : set_log_level, init_log, put_log, IntToStr
  use jlt_time, only : time_type, init_all_time, init_each_time, set_time_data, TU_SEC, TU_MIL, TU_MCR, &
                        set_time_unit, get_time_unit
  use jlt_mpi_lib, only : jml_abort, jml_AllreduceMin, jml_AllreduceMax, jml_AllreduceMaxLocal
  use jlt_grid, only : init_grid
  use jlt_exchange, only : init_exchange
  use jlt_data, only : init_data
  use jlt_remapping, only : init_remapping
  implicit none
  character(len=*), intent(IN) :: model_name ! main component name of my task 
  character(len=3), optional, intent(IN) :: default_time_unit ! 2014/07/03 [ADD]
  integer, optional, intent(IN) :: log_level ! 0, 1, 2
  logical, optional, intent(IN) :: log_stderr
  integer :: num_of_comp 
  integer :: mdl
  integer :: opt_log_level
  logical :: opt_log_stderr
  integer :: my_time_unit, tu_min, tu_max
  integer :: my_comp, max_comp
  integer :: i

  if (my_comp_name == "") then ! jlt_set_new_comp not called
     my_comp_name = trim(model_name)
     call jlt_set_new_comp(model_name)
  end if

  
  call init_model_process() ! 2014/08/27 [MOD]

  ! set time unit
  if (present(default_time_unit)) then
    select case(default_time_unit)
    case("SEC")
      call set_time_unit(TU_SEC)
    case("MIL")
      call set_time_unit(TU_MIL)
    case("MCR")
      call set_time_unit(TU_MCR)
    case default
      write(0, *) "jlt_initialize, default_time_unit setting error!!!, "//trim(default_time_unit)
      stop 999
    end select
  else
    call set_time_unit(TU_SEC)
  end if

  ! check time unit of all components
  my_time_unit = get_time_unit()
  
  call jml_AllreduceMin(my_time_unit, tu_min)
  call jml_AllreduceMax(my_time_unit, tu_max)
  
  if (tu_min /= tu_max) then
    write(0,*) "jlt_intialize, the time unit must be same in whole component"
    call jml_abort()
    stop 9999 
  end if

  call jlt_get_mpi_parameter(trim(model_name), my_comm, my_group, my_size, my_rank)
  
  num_of_comp = get_num_of_total_component()

  max_num_of_exchange_data = NUM_OF_EXCHANGE_DATA ! set initial value 2013/04/02 

  call init_all_time(num_of_comp)

  do mdl = 1, num_of_comp
    call init_each_time(mdl, 1) ! the number of domain is set to 1
  end do


  if (present(log_level)) then
    opt_log_level = log_level
  else
    opt_log_level = 0 ! default no output log
  end if

  if (present(log_stderr)) then
    opt_log_stderr = log_stderr
  else
    opt_log_stderr = .false. ! default no output stderr
  end if

  call set_log_level(opt_log_level, opt_log_stderr)

  call init_log(trim(model_name))

  my_comp_id   = get_comp_id_from_name(trim(model_name))

  call init_grid(my_comp_id, my_rank, my_size)
  call init_exchange()

  !call set_time_data(current_time, 0, 0, 0, 0, 0, int(0, kind=8))

  !current_time%delta_t = 0

  call put_log("coupler initialization completed ", 1)

  do mdl = 1, get_num_of_total_component()
    if (is_my_component(mdl)) then
      call put_log("assigned component name : "//trim(get_component_name(mdl))//", comp_id = "//trim(IntToStr(mdl)))
    end if
  end do

  call init_data(trim(model_name))
  
  call init_remapping(my_comm)
  
end subroutine jlt_initialize

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine jlt_coupling_end(time_array, isCallFinalize)
  use jlt_mpi_lib, only : jml_finalize!, jml_Send1D_m2c, jml_destruct_window
  use jlt_utils, only : finalize_log, put_log
  use jlt_grid, only : get_num_of_my_grid, get_grid_name
  use jlt_remapping, only : delete_grid_rank_file
  implicit none
  integer, intent(IN)           :: time_array(:)
  logical, optional, intent(IN) :: isCallFinalize
  integer :: i

  
  call put_log("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ", 1)
  call put_log("!!!!!!!!!!!!!   COUPLER FINALIZE  START !!!!!!!!!!!! ", 1)
  call put_log("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ", 1)

  do i = 1, get_num_of_my_grid()
     call delete_grid_rank_file(trim(my_comp_name), trim(get_grid_name(i)))
  end do
  
  call jml_finalize(isCallFinalize)

  call put_log("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ", 1)
  call put_log("!!!!!!!!!!!!!!!!  COUPLING COMPLETED !!!!!!!!!!!!!!! ", 1)
  call put_log("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ", 1)

  call finalize_log()

end subroutine jlt_coupling_end

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine jlt_get_mpi_parameter(comp_name, my_comm, my_group, my_size, my_rank)
  use jlt_mpi_lib, only : jml_GetComm, jml_GetMyGroup, jml_GetCommSizeLocal, jml_GetMyRank, jml_GetCommNULL, &
                          jml_GetCommGlobal, jml_GetMyRankGlobal, jml_GetCommSizeGlobal
  use jlt_comp, only : get_comp_id_from_name, is_my_component
  implicit none
  character(len=*), intent(IN) :: comp_name
  integer, intent(OUT) :: my_comm, my_group, my_size, my_rank
  integer :: comp_id


  if (trim(comp_name) == "GLOBAL") then
     my_comm = jml_GetCommGlobal()
     my_group = 0
     my_size = jml_GetCommSizeGlobal()
     my_rank = jml_GetMyrankGlobal()
     return
  end if
  
  comp_id = get_comp_id_from_name(comp_name)

  if (is_my_component(comp_id)) then
    my_comm = jml_GetComm(comp_id) 
  else
    my_comm = jml_GetCommNULL()
  end if

  my_group = jml_GetMyGroup(comp_id) 
  my_size  = jml_GetCommSizeLocal(comp_id) 
  my_rank  = jml_GetMyRank(comp_id)

end subroutine jlt_get_mpi_parameter


!=======+=========+=========+=========+=========+=========+=========+=========+

integer function jlt_get_myrank_global()
  use jlt_mpi_lib, only : jml_GetMyrankGlobal
  implicit none

  jlt_get_myrank_global = jml_GetMyrankGlobal()

end function jlt_get_myrank_global

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine jlt_def_grid(grid_index, model_name, grid_name, num_of_vgrid)
  use jlt_comp, only : get_num_of_my_component, is_my_component
  use jlt_grid, only : def_grid
  use jlt_utils, only : error, put_log, IntToStr
  use jlt_remapping, only : make_grid_rank_file
  implicit none
  integer, intent(IN) :: grid_index(:)
  character(len=*), intent(IN) :: model_name ! model (component) name
  character(len=*), intent(IN) :: grid_name
  integer, optional, intent(IN) :: num_of_vgrid
  integer :: i

  if (.not.is_my_component(model_name)) then
    call error("jlt_def_grid", "Component name : "//trim(model_name)//" is not defined")
  end if
  
  if (minval(grid_index) <= 0) then
    call error("jlt_def_grid", "grid index must be >= 1")
  end if

  if (present(num_of_vgrid)) then
    if (num_of_vgrid > max_num_of_exchange_data) then
      max_num_of_exchange_data = num_of_vgrid
    end if
  end if

  call def_grid(grid_index, model_name, grid_name)

  call put_log("jlt_def_grid : component name : "//trim(model_name)//", grid name : "//trim(grid_name)//", grid size : " &
                                //trim(IntToStr(size(grid_index))) &
               //", min : "//trim(IntToStr(minval(grid_index)))//", max : "//trim(IntToStr(maxval(grid_index))))

  call put_log("jlt_def_grid : make_grid_rank_file start")

  call make_grid_rank_file(model_name, grid_name, grid_index)
  
  call put_log("jlt_def_grid : make_grid_rank_file finish")

end subroutine jlt_def_grid

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine jlt_end_grid_def()
  use mpi
  implicit none
  integer :: ierror

  call mpi_barrier(MPI_COMM_WORLD, ierror)
  !write(0, *) "MPI_barrier OK"
  
end subroutine jlt_end_grid_def

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine jlt_def_varp(varp_ptr, comp_name, data_name, grid_name, num_of_layer)
  implicit none
  type(jlt_varp_type), pointer :: varp_ptr
  character(len=*), intent(IN) :: comp_name
  character(len=*), intent(IN) :: data_name
  character(len=*), intent(IN) :: grid_name
  integer, intent(IN)          :: num_of_layer

  allocate(varp_ptr)

  varp_ptr%comp_name = trim(comp_name)
  varp_ptr%grid_name = trim(grid_name)
  varp_ptr%data_name = trim(data_name)

end subroutine jlt_def_varp

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine jlt_def_varg(varg_ptr, comp_name, data_name, grid_name, num_of_data, &
     send_model_name, send_data_name, recv_mode, interval, time_lag, &
     mapping_tag, exchange_tag, time_intpl_tag)
  implicit none
  type(jlt_varg_type), pointer :: varg_ptr
  character(len=*), intent(IN) :: comp_name
  character(len=*), intent(IN) :: data_name
  character(len=*), intent(IN) :: grid_name
  integer, optional, intent(IN) :: num_of_data
  character(len=*), intent(IN) :: send_model_name
  character(len=*), intent(IN) :: send_data_name
  character(len=3), optional, intent(IN) :: recv_mode
  integer, optional, intent(IN) :: interval
  integer, optional, intent(IN) :: time_lag
  integer, optional, intent(IN) :: mapping_tag
  integer, optional, intent(IN) :: exchange_tag
  integer, optional, intent(IN) :: time_intpl_tag ! time interpolation tag 2018/07/25

  allocate(varg_ptr)

  varg_ptr%comp_name = trim(comp_name)
  varg_ptr%grid_name = trim(grid_name)
  varg_ptr%data_name = trim(data_name)

end subroutine jlt_def_varg

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine jlt_end_var_def()
  implicit none

end subroutine jlt_end_var_def

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine jlt_set_fill_value(fill_value)
  implicit none
  real(kind=8), intent(IN) :: fill_value

end subroutine jlt_set_fill_value

!=======+=========+=========+=========+=========+=========+=========+=========+

function jlt_get_fill_value() result(fill_value)
  implicit none
  real(kind=8) :: fill_value

  fill_value = -9999.d0

end function jlt_get_fill_value

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine jlt_set_mapping_table(my_model_name, &
                                 send_model_name, send_grid_name, recv_model_name,  recv_grid_name, &
                                 map_tag, is_recv_intpl, intpl_mode, send_grid, recv_grid, coef)
  use jlt_constant, only : NUM_OF_EXCHANGE_GRID, MAX_GRID, NO_GRID, STR_LONG
  use jlt_constant, only : INTPL_SERIAL_FAST, INTPL_SERIAL_SAFE, INTPL_PARALLEL, INTPL_USER
  use jlt_mpi_lib, only : jml_isLocalLeader, jml_BcastLocal, jml_SendLeader, jml_RecvLeader, jml_GetMyrank, &
                           jml_GetLeaderRank
  use jlt_utils, only : put_log, IntToStr, error
  use jlt_comp, only : get_comp_id_from_name,is_my_component
  use jlt_exchange, only : set_mapping_table
  use jlt_data, only : get_num_of_send_data, is_my_send_data, set_my_send_exchange, &
                       get_num_of_recv_data, is_my_recv_data, set_my_recv_exchange
  implicit none
  character(len=*), intent(IN)  :: my_model_name
  character(len=*), intent(IN)  :: send_model_name, send_grid_name
  character(len=*), intent(IN)  :: recv_model_name, recv_grid_name
  integer, intent(IN)           :: map_tag
  logical, intent(IN)           :: is_recv_intpl
  character(len=*), intent(IN)  :: intpl_mode ! "FAST", "SAFE", "PARALLEL", "USER"
  integer, intent(IN)           :: send_grid(:), recv_grid(:)
  real(kind=8), intent(IN)      :: coef(:)
  logical :: is_my_intpl
  integer :: intpl_mode_int
  integer :: my_model_id, send_model_id, recv_model_id
  integer :: i
  
  call put_log("------------------------------------------------------------------------------------------")
  call put_log("set mapping table start : "//trim(send_model_name)//":"//trim(recv_model_name) &
              //", grid = "//trim(send_grid_name)//":"//trim(recv_grid_name),1)


  my_model_id   = get_comp_id_from_name(trim(my_model_name))
  send_model_id = get_comp_id_from_name(trim(send_model_name))
  recv_model_id = get_comp_id_from_name(trim(recv_model_name))

  if (my_model_id == recv_model_id) then
     is_my_intpl = is_recv_intpl
  else
     is_my_intpl = .not.is_recv_intpl
  end if

  select case (trim(intpl_mode))
  case ("FAST", "fast", "Fast")
     intpl_mode_int = INTPL_SERIAL_FAST
  case("SAFE", "safe", "Safe")
     intpl_mode_int = INTPL_SERIAL_SAFE
  case("PARALLEL", "parallel", "Parallel")
     intpl_mode_int = INTPL_PARALLEL
  case("USER", "user", "User")
     intpl_mode_int = INTPL_USER
  case default
     call error("jlt_set_mapping_table", "intpl_mode msut be FAST or SAFE or PARALLEL or USER")
  end select
  
  call set_mapping_table(trim(my_model_name), trim(send_model_name), trim(send_grid_name), &
                                              trim(recv_model_name), trim(recv_grid_name), &
                                              map_tag, is_my_intpl, intpl_mode_int, send_grid, recv_grid, coef)
                                              
  call put_log("set mapping table end : "//trim(send_model_name)//":"//trim(recv_model_name),1)

  
  do i = 1, get_num_of_send_data()
     if (is_my_send_data(i, send_model_name, send_grid_name, recv_model_name, recv_grid_name)) then
        call set_my_send_exchange(i, send_model_name, send_grid_name, recv_model_name, recv_grid_name, map_tag)
     end if
  end do
  
  do i = 1, get_num_of_recv_data()
     if (is_my_recv_data(i, send_model_name, send_grid_name, recv_model_name, recv_grid_name)) then
        call set_my_recv_exchange(i, send_model_name, send_grid_name, recv_model_name, recv_grid_name, map_tag)
     end if
  end do
  
  
  call put_log("------------------------------------------------------------------------------------------")

  return

end subroutine jlt_set_mapping_table

!=======+=========+=========+=========+=========+=========+=========+=========+

!subroutine jlt_set_mapping_table(my_model_name, &
!                                 send_model_name, send_grid_name, recv_model_name,  recv_grid_name, &
!                                 send_grid, recv_grid, coef)
!  use jlt_constant, only : NUM_OF_EXCHANGE_GRID, MAX_GRID, NO_GRID, STR_LONG
!  use jlt_mpi_lib, only : jml_isLocalLeader, jml_BcastLocal, jml_SendLeader, jml_RecvLeader, jml_GetMyrank, &
!                           jml_GetLeaderRank
!  use jlt_utils, only : put_log, IntToStr, error
!  use jlt_comp, only : get_comp_id_from_name,is_my_component
!  use jlt_exchange, only : set_mapping_table
!  use jlt_namelist, only : exchange_conf_type, get_exchange_conf_ptr_name
!  implicit none
!  character(len=*), intent(IN)  :: my_model_name
!  character(len=*), intent(IN)  :: send_model_name, send_grid_name
!  character(len=*), intent(IN)  :: recv_model_name, recv_grid_name
!  integer, intent(IN)           :: send_grid(:), recv_grid(:)
!  real(kind=8), intent(IN)      :: coef(:)
!  type(exchange_conf_type), pointer :: my_exchange
  
!  integer :: my_model_id, send_model_id, recv_model_id

!  call put_log("------------------------------------------------------------------------------------------")
!  call put_log("------------------------------------------------------------------------------------------")
!  call put_log("set mapping table start : "//trim(send_model_name)//":"//trim(recv_model_name) &
!              //", grid = "//trim(send_grid_name)//":"//trim(recv_grid_name),1)

 
!  my_model_id   = get_comp_id_from_name(trim(my_model_name))
!  send_model_id = get_comp_id_from_name(trim(send_model_name))
!  recv_model_id = get_comp_id_from_name(trim(recv_model_name))

!  my_exchange => get_exchange_conf_ptr_name(send_model_name, send_grid_name, &
!                                            recv_model_name, recv_grid_name)

!  call set_mapping_table(trim(my_model_name), trim(send_model_name), trim(send_grid_name), &
!                                              trim(recv_model_name), trim(recv_grid_name), &
!                                              my_exchange%map_tag, is_my_intpl, send_grid, recv_grid, coef)
                                              
!  call put_log("set mapping table end : "//trim(send_model_name)//":"//trim(recv_model_name),1)

!  call put_log("------------------------------------------------------------------------------------------")
!  call put_log("------------------------------------------------------------------------------------------")

!  return

!end subroutine jlt_set_mapping_table

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine jlt_init_time(time_array)
  use jlt_utils, only : put_log
  !use jlt_data, only : set_my_data
  implicit none
  integer, intent(IN) :: time_array(6)
  
  call put_log("------------------------------------------------------------------------------------------")

  !call set_my_data()

  current_sec = 0
  next_sec  = 0

  call put_log("------------------------------------------------------------------------------------------")

end subroutine jlt_init_time

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine jlt_set_time(comp_name, current_time, delta_t)
  use jlt_constant, only : STR_MID
  use jlt_utils, only : put_log
  use jlt_data, only : get_num_of_recv_data, get_recv_data_name, recv_my_data, &
                       interpolate_recv_data
  use jlt_mpi_lib, only : jml_send_waitall, jml_recv_waitall
  implicit none
  character(len=*), intent(IN) :: comp_name
  integer, intent(IN) :: current_time(6)
  integer, intent(IN) :: delta_t
  character(len=STR_MID) :: log_str
  integer :: i
  
  current_sec = next_sec
 
  next_sec = current_sec + delta_t

  current_delta_t = delta_t

  call put_log("------------------------------------------------------------------------------------------")

  write(log_str,'("  ",A, I8, A, I8, A, I5)') "[jlt_set_time] set time START, current_sec =  ", current_sec, &
                                              ", next_sec = ", next_sec, ", delta_t = ", delta_t
  call put_log(trim(log_str))
  
  

  do i = 1, get_num_of_recv_data()
     call recv_my_data(get_recv_data_name(i), current_sec)
  end do

  call jml_send_waitall()
  call jml_recv_waitall()

  do i = 1, get_num_of_recv_data()
     call interpolate_recv_data(get_recv_data_name(i), current_sec)
  end do

  write(log_str,'("  ",A)') "[jlt_set_time] set time END"
  call put_log(trim(log_str))

  call put_log("------------------------------------------------------------------------------------------")

end subroutine jlt_set_time

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine jlt_put_data_1d(data_name, data, data_vector)
  use jlt_data, only : put_data_1d
  use jlt_utils, only : put_log
  implicit none
  character(len=*), intent(IN)       :: data_name
  real(kind=8), intent(IN)           :: data(:)
  real(kind=8), optional, intent(IN) :: data_vector(:)
  character(len=STR_MID) :: log_str

  call put_log("------------------------------------------------------------------------------------------")

  write(log_str,'("  ",A)') "[jlt_put_data_1d ] put data START , data_name = "//trim(data_name)
  call put_log(trim(log_str))

  call put_data_1d(data_name, data, next_sec, current_delta_t)

  write(log_str,'("  ",A)') "[jlt_put_data_1d ] put data END , data_name = "//trim(data_name)
  call put_log(trim(log_str))

  call put_log("------------------------------------------------------------------------------------------")

end subroutine jlt_put_data_1d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine jlt_put_data_1d_ptr(varp_ptr, data, data_vector)
  use jlt_data, only : put_data_1d
  implicit none
  type(jlt_varp_type), pointer       :: varp_ptr
  real(kind=8), intent(IN)           :: data(:)
  real(kind=8), optional, intent(IN) :: data_vector(:)
  
  call put_data_1d(trim(varp_ptr%data_name), data, next_sec, current_delta_t)

end subroutine jlt_put_data_1d_ptr

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine jlt_get_data_1d(data_name, data, data_vector, is_recv_ok)
  use jlt_data, only : get_data_1d
  use jlt_utils, only : put_log
  implicit none
  character(len=*),intent(IN)           :: data_name
  real(kind=8), intent(INOUT)           :: data(:)
  real(kind=8), optional, intent(INOUT) :: data_vector(:)
  logical, intent(INOUT)                :: is_recv_ok
  character(len=STR_MID) :: log_str

  call put_log("------------------------------------------------------------------------------------------")

  write(log_str,'("  ",A)') "[jlt_get_data_1d ] get data START , data_name = "//trim(data_name)
  call put_log(trim(log_str))

  call get_data_1d(data_name, data, current_sec, is_recv_ok)

  write(log_str,'("  ",A)') "[jlt_get_data_1d ] get data END , data_name = "//trim(data_name)
  call put_log(trim(log_str))

  call put_log("------------------------------------------------------------------------------------------")

end subroutine jlt_get_data_1d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine jlt_get_data_1d_ptr(varg_ptr, data, data_vector, is_recv_ok)
  use jlt_data, only : get_data_1d
  implicit none
  type(jlt_varg_type), pointer          :: varg_ptr
  real(kind=8), intent(INOUT)           :: data(:)
  real(kind=8), optional, intent(INOUT) :: data_vector(:)
  logical, intent(INOUT)                :: is_recv_ok

  call get_data_1d(trim(varg_ptr%data_name), data, current_sec, is_recv_ok)

end subroutine jlt_get_data_1d_ptr


!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine jlt_put_data_25d(data_name, data, data_vector)
  use jlt_data, only : put_data_2d
  use jlt_utils, only : put_log
  implicit none
  character(len=*), intent(IN)       :: data_name
  real(kind=8), intent(IN)           :: data(:,:)
  real(kind=8), optional, intent(IN) :: data_vector(:)
  character(len=STR_MID) :: log_str
  
  call put_log("------------------------------------------------------------------------------------------")

  write(log_str,'("  ",A)') "[jlt_put_data_25d ] put data START , data_name = "//trim(data_name)
  call put_log(trim(log_str))

  !write(0, *) "jlt_put_data_25d ", trim(data_name), ", ", minval(data), maxval(data)
  
  call put_data_2d(data_name, data, next_sec, current_delta_t)

  write(log_str,'("  ",A)') "[jlt_put_data_25d ] put data END , data_name = "//trim(data_name)
  call put_log(trim(log_str))

  call put_log("------------------------------------------------------------------------------------------")

end subroutine jlt_put_data_25d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine jlt_put_data_25d_ptr(varp_ptr, data, data_vector)
  use jlt_data, only : put_data_2d
  implicit none
  type(jlt_varp_type), pointer       :: varp_ptr
  real(kind=8), intent(IN)           :: data(:,:)
  real(kind=8), optional, intent(IN) :: data_vector(:)
  
  call put_data_2d(trim(varp_ptr%data_name), data, next_sec, current_delta_t)

end subroutine jlt_put_data_25d_ptr

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine jlt_get_data_25d(data_name, data, data_vector, is_recv_ok)
  use jlt_data, only : get_data_2d
  use jlt_utils, only : put_log
  implicit none
  character(len=*),intent(IN)           :: data_name
  real(kind=8), intent(INOUT)           :: data(:,:)
  real(kind=8), optional, intent(INOUT) :: data_vector(:)
  logical, intent(INOUT)                :: is_recv_ok
  character(len=STR_MID) :: log_str
  
  call put_log("------------------------------------------------------------------------------------------")

  write(log_str,'("  ",A)') "[jlt_get_data_25d ] get data START , data_name = "//trim(data_name)
  call put_log(trim(log_str))

  call get_data_2d(data_name, data, current_sec, is_recv_ok)

  write(log_str,'("  ",A)') "[jlt_get_data_25d ] get data END , data_name = "//trim(data_name)
  call put_log(trim(log_str))

  call put_log("------------------------------------------------------------------------------------------")

end subroutine jlt_get_data_25d

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine jlt_get_data_25d_ptr(varg_ptr, data, data_vector, is_recv_ok)
  use jlt_data, only : get_data_2d
  implicit none
  type(jlt_varg_type), pointer          :: varg_ptr
  real(kind=8), intent(INOUT)           :: data(:,:)
  real(kind=8), optional, intent(INOUT) :: data_vector(:)
  logical, intent(INOUT)                :: is_recv_ok
  
  call get_data_2d(trim(varg_ptr%data_name), data, current_sec, is_recv_ok)

end subroutine jlt_get_data_25d_ptr


!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine jlt_error(sub_name, error_str)
  use jlt_utils, only : Error
  character(len=*), intent(IN) :: sub_name
  character(len=*), intent(IN) :: error_str

  call Error(sub_name, error_str)

end subroutine jlt_error


!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine bcast_array_local_int(my_comp_name,  data)
  use jlt_mpi_lib, only : jml_BcastLocal
  use jlt_comp, only : get_comp_id_from_name
  implicit none
  integer, intent(INOUT) :: data(:)
  character(len=*), intent(IN) :: my_comp_name
  integer :: my_model

  my_model = get_comp_id_from_name(my_comp_name)

  call jml_BcastLocal(my_model, data, 1, size(data))

end subroutine bcast_array_local_int

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine send_array_to_recv_model_str(my_comp_name, recv_comp_name, array)
  use jlt_mpi_lib, only : jml_isLocalLeader, jml_SendLeader
  use jlt_comp, only : get_comp_id_from_name
  implicit none
  character(len=*), intent(IN) :: my_comp_name, recv_comp_name
  character(len=*), intent(IN) :: array
  integer :: my_model
  integer :: recv_model
  integer :: array_size(1)

  my_model = get_comp_id_from_name(my_comp_name)
  recv_model = get_comp_id_from_name(recv_comp_name)

  if (jml_isLocalLeader(my_model)) then
    call jml_SendLeader(array, recv_model-1) 
  end if

end subroutine send_array_to_recv_model_str

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine recv_array_from_send_model_str(my_comp_name, send_comp_name, array, bcast_flag)
  use jlt_mpi_lib, only : jml_isLocalLeader, jml_RecvLeader, jml_BcastLocal
  use jlt_comp, only : get_comp_id_from_name
  implicit none
  character(len=*), intent(IN) :: my_comp_name, send_comp_name
  character(len=*), intent(INOUT) :: array
  logical, optional, intent(IN) :: bcast_flag
  integer :: my_model
  integer  :: send_model
  logical :: is_bcast = .true.

  is_bcast = .true.

  if (present(bcast_flag)) is_bcast = bcast_flag

  my_model = get_comp_id_from_name(my_comp_name)
  send_model = get_comp_id_from_name(send_comp_name)

  if (jml_isLocalLeader(my_model)) then
    call jml_RecvLeader(array, send_model-1)
  end if

  if (is_bcast) call jml_BcastLocal(my_model, array)

end subroutine recv_array_from_send_model_str

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine send_array_to_recv_model_int(my_comp_name, recv_comp_name, array)
  use jlt_mpi_lib, only : jml_isLocalLeader, jml_SendLeader
  use jlt_comp, only : get_comp_id_from_name
  implicit none
  character(len=*), intent(IN) :: my_comp_name, recv_comp_name
  integer, intent(IN) :: array(:)
  integer :: my_model
  integer :: recv_model
  integer :: array_size(1)

  my_model = get_comp_id_from_name(my_comp_name)
  recv_model = get_comp_id_from_name(recv_comp_name)

  if (jml_isLocalLeader(my_model)) then
    call jml_SendLeader(array, 1, size(array), recv_model-1) 
  end if

end subroutine send_array_to_recv_model_int

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine recv_array_from_send_model_int(my_comp_name, send_comp_name, array, bcast_flag)
  use jlt_mpi_lib, only : jml_isLocalLeader, jml_RecvLeader, jml_BcastLocal
  use jlt_comp, only : get_comp_id_from_name
  implicit none
  character(len=*), intent(IN) :: my_comp_name, send_comp_name
  integer, intent(INOUT) :: array(:)
  logical, optional, intent(IN) :: bcast_flag
  integer :: my_model
  integer  :: send_model
  logical :: is_bcast = .true.

  is_bcast = .true.
  
  if (present(bcast_flag)) is_bcast = bcast_flag

  my_model = get_comp_id_from_name(my_comp_name)
  send_model = get_comp_id_from_name(send_comp_name)

  if (jml_isLocalLeader(my_model)) then
    call jml_RecvLeader(array, 1, size(array), send_model-1)
  end if

  if (is_bcast) call jml_BcastLocal(my_model, array, 1, size(array))

end subroutine recv_array_from_send_model_int

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine send_array_to_recv_model_real(my_comp_name, recv_comp_name, array)
  use jlt_mpi_lib, only : jml_isLocalLeader, jml_SendLeader
  use jlt_comp, only : get_comp_id_from_name
  implicit none
  character(len=*), intent(IN) :: my_comp_name, recv_comp_name
  real(kind=4), intent(IN) :: array(:)
  integer :: my_model
  integer :: recv_model
  integer :: array_size(1)

  my_model = get_comp_id_from_name(my_comp_name)
  recv_model = get_comp_id_from_name(recv_comp_name)

  if (jml_isLocalLeader(my_model)) then
    call jml_SendLeader(array, 1, size(array), recv_model-1) 
  end if

end subroutine send_array_to_recv_model_real

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine recv_array_from_send_model_real(my_comp_name, send_comp_name, array, bcast_flag)
  use jlt_mpi_lib, only : jml_isLocalLeader, jml_RecvLeader, jml_BcastLocal
  use jlt_comp, only : get_comp_id_from_name
  implicit none
  character(len=*), intent(IN) :: my_comp_name, send_comp_name
  real(kind=4), intent(INOUT) :: array(:)
  logical, optional, intent(IN) :: bcast_flag
  integer :: my_model
  integer :: send_model
  logical :: is_bcast = .true.

  is_bcast = .true.

  if (present(bcast_flag)) is_bcast = bcast_flag

  my_model = get_comp_id_from_name(my_comp_name)
  send_model = get_comp_id_from_name(send_comp_name)

  if (jml_isLocalLeader(my_model)) then
    call jml_RecvLeader(array, 1, size(array), send_model-1)
  end if

  if (is_bcast) call jml_BcastLocal(my_model, array, 1, size(array))

end subroutine recv_array_from_send_model_real

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine send_array_to_recv_model_dbl(my_comp_name, recv_comp_name, array)
  use jlt_mpi_lib, only : jml_isLocalLeader, jml_SendLeader
  use jlt_comp, only : get_comp_id_from_name
  implicit none
  character(len=*), intent(IN) :: my_comp_name, recv_comp_name
  real(kind=8), intent(IN) :: array(:)
  integer :: my_model
  integer :: recv_model
  integer :: array_size(1)

  my_model = get_comp_id_from_name(my_comp_name)
  recv_model = get_comp_id_from_name(recv_comp_name)

  if (jml_isLocalLeader(my_model)) then
    call jml_SendLeader(array, 1, size(array), recv_model-1) 
  end if

end subroutine send_array_to_recv_model_dbl

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine recv_array_from_send_model_dbl(my_comp_name, send_comp_name, array, bcast_flag)
  use jlt_mpi_lib, only : jml_isLocalLeader, jml_RecvLeader, jml_BcastLocal
  use jlt_comp, only : get_comp_id_from_name
  implicit none
  character(len=*), intent(IN) :: my_comp_name, send_comp_name
  real(kind=8), intent(INOUT) :: array(:)
  logical, optional, intent(IN) :: bcast_flag
  integer :: my_model
  integer  :: send_model
  logical :: is_bcast = .true.

  is_bcast = .true.

  if (present(bcast_flag)) is_bcast = bcast_flag

  my_model = get_comp_id_from_name(my_comp_name)
  send_model = get_comp_id_from_name(send_comp_name)

  if (jml_isLocalLeader(my_model)) then
    call jml_RecvLeader(array, 1, size(array), send_model-1)
  end if

  if (is_bcast) call jml_BcastLocal(my_model, array, 1, size(array))

end subroutine recv_array_from_send_model_dbl


!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine jlt_log(sub_name, error_str, log_level)
  use jlt_utils, only : put_log
  character(len=*), intent(IN) :: sub_name
  character(len=*), intent(IN) :: error_str
  integer, optional, intent(IN) :: log_level
  integer :: ll

  ll = 2
  if (present(log_level)) ll = log_level

  call put_log("Sub["//trim(sub_name)//"] : "//trim(error_str), ll)

end subroutine jlt_log

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine jlt_inc_calendar(itime, del_t)
  use jlt_time, only : get_current_time, time_type, get_delta_t, &
                        inc_calendar, get_time_unit, TU_SEC, TU_MIL, TU_MCR
  use jlt_comp, only : get_comp_id_from_name
  use jlt_utils, only : error
  implicit none
  integer, intent(INOUT) :: itime(:)
  integer, intent(IN) :: del_t
  type(time_type) :: time
  integer(kind=8) :: time_sec

  !write(0,*) "jcup_inc_time 1 "//trim(component_name)//", ",itime
  time%yyyy = itime(1)
  time%mo   = itime(2)
  time%dd   = itime(3)
  time%hh   = itime(4)
  time%mm   = itime(5)
  time%ss   = itime(6)
  
  select case(get_time_unit())
  case(TU_SEC)
  case(TU_MIL)
    if (size(itime) < 7) call error("jlt_inc_time", "array size of itime must be >= 7")
    time%milli_sec = itime(7)
  case(TU_MCR)
    if (size(itime) < 8) call error("jlt_inc_time", "array size of itime must be >= 8")
    time%milli_sec = itime(7)
    time%micro_sec = itime(8)
  case default
    call error("jlt_inc_time", "time unit parameter error")
  end select

  call inc_calendar(time, del_t)

  itime(1) = time%yyyy
  itime(2) = time%mo
  itime(3) = time%dd
  itime(4) = time%hh
  itime(5) = time%mm
  itime(6) = time%ss

  select case(get_time_unit())
  case(TU_SEC)
  case(TU_MIL)
    if (size(itime) < 7) call error("jlt_inc_time", "array size of itime must be >= 7")
    itime(7) = time%milli_sec
  case(TU_MCR)
    if (size(itime) < 8) call error("jlt_inc_time", "array size of itime must be >= 8")
    itime(7) = time%milli_sec
    itime(8) = time%micro_sec
  case default
    call error("jlt_inc_time", "time unit parameter error")
  end select

  !write(0,*) "jlt_inc_time 2 "//trim(component_name)//", ",itime

end subroutine jlt_inc_calendar

!=======+=========+=========+=========+=========+=========+=========+=========+
subroutine jlt_inc_time(component_name, itime)
  use jlt_time, only : get_current_time, time_type, get_delta_t, &
                        inc_time, get_time_unit, TU_SEC, TU_MIL, TU_MCR
  use jlt_comp, only : get_comp_id_from_name
  use jlt_utils, only : error
  implicit none
  character(len=*), intent(IN) :: component_name
  integer, intent(INOUT) :: itime(:)
  type(time_type) :: time
  integer(kind=8) :: time_sec
  integer :: del_t
  integer :: comp_id

  comp_id = get_comp_id_from_name(component_name)

  call get_current_time(comp_id, 1, time)
  call get_delta_t(comp_id, 1, del_t)

  !write(0,*) "jlt_inc_time 1 "//trim(component_name)//", ",itime

  !call inc_time(time, del_t)

  itime(1) = time%yyyy
  itime(2) = time%mo
  itime(3) = time%dd
  itime(4) = time%hh
  itime(5) = time%mm
  itime(6) = time%ss

  select case(get_time_unit())
  case(TU_SEC)
  case(TU_MIL)
    if (size(itime) < 7) call error("jlt_inc_time", "array size of itime must be >= 7")
    itime(7) = time%milli_sec
  case(TU_MCR)
    if (size(itime) < 8) call error("jlt_inc_time", "array size of itime must be >= 8")
    itime(7) = time%milli_sec
    itime(8) = time%micro_sec
  case default
    call error("jlt_inc_time", "time unit parameter error")
  end select


end subroutine jlt_inc_time

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine jlt_write_mapping_table(fid)
  use jlt_exchange, only : write_exchange
  implicit none
  integer, intent(IN) :: fid

  call write_exchange(fid)
  
end subroutine jlt_write_mapping_table

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine jlt_read_mapping_table(fid)
  use jlt_constant, only : STR_SHORT
  use jlt_exchange, only : read_exchange, get_num_of_exchange, get_exchange_ptr
  use jlt_exchange_class, only : exchange_class
  use jlt_data, only : get_num_of_send_data, is_my_send_data, set_my_send_exchange, &
                       get_num_of_recv_data, is_my_recv_data, set_my_recv_exchange
  implicit none
  integer, intent(IN) :: fid
  character(len=STR_SHORT) :: send_model_name
  character(len=STR_SHORT) :: send_grid_name
  character(len=STR_SHORT) :: recv_model_name
  character(len=STR_SHORT) :: recv_grid_name
  integer                  :: map_tag
  class(exchange_class), pointer :: exchange_ptr
  integer :: i, j
  
  call read_exchange(fid)

  do i = 1, get_num_of_exchange()
     exchange_ptr => get_exchange_ptr(i)
     send_model_name = trim(exchange_ptr%get_send_comp_name())
     send_grid_name  = trim(exchange_ptr%get_send_grid_name())
     recv_model_name = trim(exchange_ptr%get_recv_comp_name())
     recv_grid_name  = trim(exchange_ptr%get_recv_grid_name())
     map_tag         = exchange_ptr%get_map_tag()
     
     do j = 1, get_num_of_send_data()
        if (is_my_send_data(j, send_model_name, send_grid_name, recv_model_name, recv_grid_name)) then
           call set_my_send_exchange(j, send_model_name, send_grid_name, recv_model_name, recv_grid_name, map_tag)
        end if
     end do
  
     do j = 1, get_num_of_recv_data()
        if (is_my_recv_data(j, send_model_name, send_grid_name, recv_model_name, recv_grid_name)) then
           call set_my_recv_exchange(j, send_model_name, send_grid_name, recv_model_name, recv_grid_name, map_tag)
        end if
     end do
  end do
  
end subroutine jlt_read_mapping_table

!=======+=========+=========+=========+=========+=========+=========+=========+

subroutine jlt_set_user_interpolation(send_comp, send_grid, recv_comp, recv_grid, map_tag, user_func)
  use jlt_exchange_class, only : exchange_class, interpolation_user_ifc
  use jlt_exchange, only : get_exchange_ptr
  implicit none
  character(len=*), intent(IN) :: send_comp
  character(len=*), intent(IN) :: send_grid
  character(len=*), intent(IN) :: recv_comp
  character(len=*), intent(IN) :: recv_grid
  integer, intent(IN)          :: map_tag
  procedure(interpolation_user_ifc), pointer, intent(IN) :: user_func
  type(exchange_class), pointer :: exchange_ptr

  exchange_ptr => get_exchange_ptr(send_comp, send_grid, recv_comp, recv_grid, map_tag)

  call exchange_ptr%set_user_interpolation(user_func)
  
end subroutine jlt_set_user_interpolation

  
!=======+=========+=========+=========+=========+=========+=========+=========+

end module jlt_interface
