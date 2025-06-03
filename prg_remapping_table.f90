program remapping_table
  use mpi
  use jlt_remapping
  implicit none
  integer, parameter :: GRID_SIZE = 10
  integer, parameter :: GLOBAL_MAP_SIZE = 20
  integer :: grid_index(GRID_SIZE)
  integer :: grid_rank(GRID_SIZE)
  integer :: global_index(GLOBAL_MAP_SIZE)
  integer :: global_target(GLOBAL_MAP_SIZE)
  real(kind=8) :: global_coef(GLOBAL_MAP_SIZE)
  integer, pointer :: local_index(:)
  integer, pointer :: local_target(:)
  real(kind=8), pointer :: local_coef(:)
  integer          :: local_size
  integer          :: local_rank
  integer :: ierror
  integer :: i

  call mpi_init(ierror)
  call mpi_comm_size(MPI_COMM_WORLD, local_size, ierror)
  call mpi_comm_rank(MPI_COMM_WORLD, local_rank, ierror)
  call init_remapping(MPI_COMM_WORLD)

  do i = 1, GRID_SIZE
     grid_index(i) = 100 + i - (GRID_SIZE * local_rank)
  end do

  call make_grid_rank_file("AGCM", "agcm_grid", grid_index)
  call get_target_grid_rank("AGCM", "agcm_grid", grid_index, grid_rank)
  
  call mpi_finalize(ierror)
  stop
  
  do i = 1, GLOBAL_MAP_SIZE
     global_index(i) = 2*i
     global_target(i) = i
     global_coef(i)  = i
  end do

  call make_local_mapping_table(global_index, global_target, global_coef, &
                                grid_index, &
                                local_size, local_index, local_target, local_coef)

  write(*,*) local_size
  write(*,*) local_index
  
end program remapping_table
