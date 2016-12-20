module write_plotfile_n_module

  use ml_layout_module
  use multifab_module
  use probin_common_module , only : prob_lo, prob_hi
  use probin_reactdiff_module, only: nspecies
  use fabio_module

  implicit none

  private

  public :: write_plotfile_n

contains

  subroutine write_plotfile_n(mla,n_in,dx,time,istep_to_write)

    type(ml_layout)  , intent(in   ) :: mla
    type(multifab)   , intent(in   ) :: n_in(:)
    integer          , intent(in   ) :: istep_to_write
    real(kind=dp_t)  , intent(in   ) :: dx(:,:), time

    integer :: i,n,n_plot_comps

    character(len=11) :: sd_name

    character(len=20), allocatable :: plot_names(:)

    ! These are only used if you want to coarsen your plotdata before writing
    ! Start crse
    type(multifab), allocatable :: plotdata(:)          ! cell-centered

    integer :: nlevs,dm

    type(bl_prof_timer), save :: bpt

    call build(bpt,"write_plotfile_n")

    nlevs = mla%nlevel
    dm = mla%dim

    n_plot_comps = nspecies

    allocate(plot_names(n_plot_comps))

    do n=1,nspecies
       write(plot_names(n),'(a,i0)') "n", n
    enddo

    allocate(plotdata(nlevs))

    do n=1,nlevs
       call multifab_build(plotdata(n), mla%la(n), n_plot_comps, 0)
    end do

    do n = 1,nlevs
       call multifab_copy_c(plotdata(n),1,n_in(n),1,nspecies)
    end do

    write(unit=sd_name,fmt='("plt",i8.8)') istep_to_write
    if ( parallel_IOProcessor() ) then
      write(*,'(2A)') "Saving PLOT FILEs to directory ", trim(sd_name)
      write(*,*)
    end if

    ! cell-centered plotfile
    call fabio_ml_multifab_write_d(plotdata, mla%mba%rr(:,1), sd_name, plot_names, &
                                   mla%mba%pd(1), prob_lo, prob_hi, time, dx(1,:))

    do n=1,nlevs
      call multifab_destroy(plotdata(n))
    end do
    deallocate(plotdata)
    deallocate(plot_names)

    call destroy(bpt)

  end subroutine write_plotfile_n

end module write_plotfile_n_module
