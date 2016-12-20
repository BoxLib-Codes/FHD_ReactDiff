module checkpoint_module

  use parallel
  use bl_types
  use multifab_module
  use ml_layout_module
  use bl_IO_module
  use fab_module
  use bl_rng_module
  use bl_random_module
  use fabio_module, only: fabio_mkdir, fabio_ml_multifab_write_d
  use probin_common_module, only: dim_in
  use probin_reactdiff_module, only: nspecies, use_bl_rng

  implicit none

  private

  public :: checkpoint_write, checkpoint_read

contains

  subroutine checkpoint_write(mla,n_new,time,dt,istep_to_write)
    
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: n_new(:)
    real(kind=dp_t), intent(in   ) :: time,dt
    integer        , intent(in   ) :: istep_to_write

    type(multifab), pointer :: chkdata(:)

    integer :: n,nlevs,dm

    character(len=11) :: sd_name
    character(len=40) :: rand_name

    type(bl_prof_timer), save :: bpt

    call build(bpt, "checkpoint_write")

    nlevs = mla%nlevel
    dm = mla%dim

    allocate(chkdata(nlevs))
    do n = 1,nlevs
       call multifab_build(chkdata(n), mla%la(n), nspecies, 0)
       call multifab_copy_c(chkdata(n), 1, n_new(n), 1, nspecies)
    end do
    write(unit=sd_name,fmt='("chk",i8.8)') istep_to_write

    call checkpoint_write_doit(nlevs, sd_name, chkdata, mla%mba%rr, time, dt)
    
    ! random state
    if (use_bl_rng) then

       ! engines
       rand_name = sd_name//'/rng_eng_diff'
       call bl_rng_save_engine(rng_eng_diffusion, rand_name)
       rand_name = sd_name//'/rng_eng_react'
       call bl_rng_save_engine(rng_eng_reaction, rand_name)

    end if

    do n = 1,nlevs
       call multifab_destroy(chkdata(n))
    end do
    deallocate(chkdata)

    call destroy(bpt)

  contains

    subroutine checkpoint_write_doit(nlevs_in, dirname, mfs, rrs, time_in, dt_in)
      
      integer         , intent(in) :: nlevs_in
      type(multifab)  , intent(in) :: mfs(:)
      integer         , intent(in) :: rrs(:,:)
      character(len=*), intent(in) :: dirname
      real(kind=dp_t) , intent(in) :: time_in, dt_in

      integer :: n
      character(len=128) :: header, sd_name
      integer :: un

      integer         :: nlevs
      real(kind=dp_t) :: time, dt

      type(bl_prof_timer), save :: bpt

      namelist /chkpoint/ time
      namelist /chkpoint/ dt
      namelist /chkpoint/ nlevs

      call build(bpt,"checkpoint_write_doit")

      if ( parallel_IOProcessor() ) call fabio_mkdir(dirname)

      call parallel_barrier() ! All CPUs have to wait till the directory is built.

      write(unit=sd_name, fmt='(a,"/State")') trim(dirname)
      call fabio_ml_multifab_write_d(mfs, rrs(:,1), sd_name)

      if (parallel_IOProcessor()) then
         print *,'Writing cc state to checkpoint file ',trim(sd_name)
         print *,' '
      end if

      time  = time_in
      dt = dt_in
      nlevs = nlevs_in

      if (parallel_IOProcessor()) then
         header = "Header"
         un = unit_new()
         open(unit=un, &
              file = trim(dirname) // "/" // trim(header), &
              form = "formatted", access = "sequential", &
              status = "replace", action = "write")
         write(unit=un, nml = chkpoint)
         do n = 1,nlevs-1
            write(unit=un,fmt=*) rrs(n,1)
         end do
         close(un)
      end if

      call destroy(bpt)

    end subroutine checkpoint_write_doit

  end subroutine checkpoint_write

  subroutine checkpoint_read(mfs, dirname, rrs_out, time_out, dt_out, nlevs_out)

    use bl_IO_module
    use fab_module
    use fabio_module, only: fabio_ml_multifab_read_d
    use parallel

    type(multifab)  ,                pointer :: mfs(:)
    character(len=*), intent(in   )          :: dirname
    integer         , intent(  out)          :: nlevs_out
    real(kind=dp_t) , intent(  out)          :: time_out, dt_out
    integer         , intent(  out)          :: rrs_out(:)

    character(len=128) :: header, sd_name

    integer :: n, un, nlevs
    integer, pointer :: rrs(:)

    real(kind=dp_t) :: time, dt

    type(bl_prof_timer), save :: bpt

    namelist /chkpoint/ nlevs
    namelist /chkpoint/ time
    namelist /chkpoint/ dt

    call build(bpt,"checkpoint_read")

!   First read the header information
    header = "Header"
    un = unit_new()
    open(unit=un, &
         file = trim(dirname) // "/" // trim(header), &
         status = "old", &
         action = "read")
    read(unit=un, nml = chkpoint)
    allocate(rrs(nlevs-1))
    do n = 1,nlevs-1
       read(unit=un,fmt=*) rrs(n)
    end do
    close(un)

     time_out = time
       dt_out = dt
    nlevs_out = nlevs

!   Read the state data into a multilevel multifab.
    write(unit=sd_name, fmt='(a,"/State")') trim(dirname)
    call fabio_ml_multifab_read_d(mfs, sd_name)

    rrs_out(1:nlevs-1) = rrs(1:nlevs-1)

    deallocate(rrs)

    call destroy(bpt)

  end subroutine checkpoint_read

end module checkpoint_module
