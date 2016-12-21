module histogram_module

  use multifab_module
  use ml_layout_module
  use probin_common_module, only: n_cells
  use probin_reactdiff_module, only: nspecies

  implicit none

  private
  public :: init_histogram, destroy_histogram, write_histogram

  type(multifab), save ::  n_serial
  type(layout), save :: la_serial

contains

  subroutine init_histogram(mla)

    type(ml_layout), intent(in   ) :: mla

    type(box)  :: bx_serial
    type(boxarray)  :: bxa_serial
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: dm

    dm = mla%dim

    ! first set lo and hi to the problem domain
    lo = 0
    hi(1:dm) = n_cells(1:dm) - 1

    ! build n_serial
    ! set bx_serial to reduced dimension problem domain
    call box_build_2(bx_serial,lo,hi)
    ! build a boxarray containing only one box
    call boxarray_build_bx(bxa_serial,bx_serial)
    call layout_build_ba(la_serial, bxa_serial, bx_serial, mla%pmask, &
                         explicit_mapping=(/parallel_IOProcessorNode()/) )
    call destroy(bxa_serial)
    call multifab_build(n_serial, la_serial, nspecies, 0)

  end subroutine init_histogram

  subroutine destroy_histogram()

    call multifab_destroy(n_serial)
    call destroy(la_serial)

  end subroutine destroy_histogram

  subroutine write_histogram(n_in,time,step)

    type(multifab) , intent(in) :: n_in(:)
    real(kind=dp_t), intent(in) :: time
    integer        , intent(in) :: step

    real(kind=dp_t), pointer :: np(:,:,:,:)
    integer :: i,j,k

    character(len=12) :: hist_name

    write(unit=hist_name,fmt='("hist",i8.8)') step

    OPEN(100,FILE=hist_name)

    call copy(n_serial,1,n_in(1),1,nspecies)

    if (parallel_IOProcessor()) then

       np => dataptr(n_serial,1)

       do k=lbound(np,3), ubound(np,3)
          do j=lbound(np,2), ubound(np,2)
             do i=lbound(np,1), ubound(np,1)     
                write(100,*) real(time), real(np(i,j,k,:)) 
             end do
          end do
       end do
       
    end if

  end subroutine write_histogram


end module histogram_module
