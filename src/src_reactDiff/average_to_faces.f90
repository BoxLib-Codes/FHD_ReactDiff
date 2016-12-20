module average_to_faces_module

  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  use probin_reactdiff_module, only: avg_type

  implicit none

  private

  public :: average_to_faces

contains

  subroutine average_to_faces(mla,n_cc,n_fc,incomp,outcomp,numcomp,dv)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: n_cc(:)
    type(multifab) , intent(inout) :: n_fc(:,:)
    integer        , intent(in   ) :: incomp,outcomp,numcomp
    real(kind=dp_t), intent(in   ) :: dv

    integer :: i,dm,n,nlevs,lo(mla%dim),hi(mla%dim)
    integer :: ng_c,ng_f

    real(kind=dp_t), pointer :: cp(:,:,:,:)
    real(kind=dp_t), pointer :: fx(:,:,:,:)
    real(kind=dp_t), pointer :: fy(:,:,:,:)
    real(kind=dp_t), pointer :: fz(:,:,:,:)

    type(mfiter) :: mfi
    type(box) :: xnodalbox, ynodalbox, znodalbox
    integer :: xlo(mla%dim), xhi(mla%dim)
    integer :: ylo(mla%dim), yhi(mla%dim)
    integer :: zlo(mla%dim), zhi(mla%dim)

    type(bl_prof_timer),save :: bpt

    call build(bpt,"average_to_faces")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_c = n_cc(1)%ng
    ng_f = n_fc(1,1)%ng
    
    !$omp parallel private(mfi,n,i,xnodalbox,ynodalbox,znodalbox,xlo,ylo,zlo) &
    !$omp private(xhi,yhi,zhi,cp,fx,fy,fz,lo,hi)
    do n=1,nlevs
       call mfiter_build(mfi, n_cc(n), tiling=.true.)
       do while (more_tile(mfi))
          i = get_fab_index(mfi)

          xnodalbox = get_nodaltilebox(mfi,1)
          xlo = lwb(xnodalbox)
          xhi = upb(xnodalbox)
          ynodalbox = get_nodaltilebox(mfi,2)
          ylo = lwb(ynodalbox)
          yhi = upb(ynodalbox)
          znodalbox = get_nodaltilebox(mfi,3)
          zlo = lwb(znodalbox)
          zhi = upb(znodalbox)

          cp => dataptr(n_cc(n),i)
          fx => dataptr(n_fc(n,1),i)
          fy => dataptr(n_fc(n,2),i)
          lo = lwb(get_box(n_cc(n),i))
          hi = upb(get_box(n_cc(n),i))
          select case (dm)
          case (2)
             call average_to_faces_2d(cp(:,:,1,:),ng_c, &
                                      fx(:,:,1,:),fy(:,:,1,:),ng_f, lo,hi, &
                                      xlo,xhi,ylo,yhi, incomp,outcomp,numcomp,dv)
          case (3)
             fz => dataptr(n_fc(n,3),i)
             call average_to_faces_3d(cp(:,:,:,:),ng_c, &
                                      fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:),ng_f, lo,hi, &
                                      xlo,xhi,ylo,yhi,zlo,zhi, incomp,outcomp,numcomp,dv)
          end select
       end do
    end do
    !$omp end parallel

    ! sync the n_fces at the boundaries
    do n=1,nlevs
       do i=1,dm
          call multifab_fill_boundary(n_fc(n,i))  
       end do
    end do

    call destroy(bpt)

  end subroutine average_to_faces

  subroutine average_to_faces_2d(n_cc,ng_c,n_fcx,n_fcy,ng_f,glo,ghi, &
                                 xlo,xhi,ylo,yhi,incomp,outcomp,numcomp,dv)

    integer        , intent(in   ) :: glo(:),ghi(:),xlo(:),xhi(:),ylo(:),yhi(:)
    integer        , intent(in   ) :: ng_c,ng_f,incomp,outcomp,numcomp
    real(kind=dp_t), intent(in   ) ::  n_cc(glo(1)-ng_c:,glo(2)-ng_c:,:)
    real(kind=dp_t), intent(inout) :: n_fcx(glo(1)-ng_f:,glo(2)-ng_f:,:)
    real(kind=dp_t), intent(inout) :: n_fcy(glo(1)-ng_f:,glo(2)-ng_f:,:)
    real(kind=dp_t), intent(in   ) :: dv

    integer :: i,j,comp

    do comp=0,numcomp-1

       ! x-faces
       do j=xlo(2),xhi(2)
          do i=xlo(1),xhi(1)
             n_fcx(i,j,outcomp+comp) = average_values(n_cc(i-1,j,incomp+comp), n_cc(i,j,incomp+comp),dv)
          end do
       end do

       ! y-faces
       do j=ylo(2),yhi(2)
          do i=ylo(1),yhi(1)
             n_fcy(i,j,outcomp+comp) = average_values(n_cc(i,j-1,incomp+comp), n_cc(i,j,incomp+comp),dv)
          end do
       end do

    end do

  end subroutine average_to_faces_2d

  subroutine average_to_faces_3d(n_cc,ng_c,n_fcx,n_fcy,n_fcz,ng_f,glo,ghi, &
                                 xlo,xhi,ylo,yhi,zlo,zhi,incomp,outcomp,numcomp,dv)

    integer        , intent(in   ) :: glo(:),ghi(:),xlo(:),xhi(:),ylo(:),yhi(:),zlo(:),zhi(:)
    integer        , intent(in   ) :: ng_c,ng_f,incomp,outcomp,numcomp
    real(kind=dp_t), intent(in   ) ::  n_cc(glo(1)-ng_c:,glo(2)-ng_c:,glo(3)-ng_c:,:)
    real(kind=dp_t), intent(inout) :: n_fcx(glo(1)-ng_f:,glo(2)-ng_f:,glo(3)-ng_f:,:)
    real(kind=dp_t), intent(inout) :: n_fcy(glo(1)-ng_f:,glo(2)-ng_f:,glo(3)-ng_f:,:)
    real(kind=dp_t), intent(inout) :: n_fcz(glo(1)-ng_f:,glo(2)-ng_f:,glo(3)-ng_f:,:)
    real(kind=dp_t), intent(in   ) :: dv

    integer :: i,j,k,comp

    do comp=0,numcomp-1

       ! x-faces
       do k=xlo(3),xhi(3)
          do j=xlo(2),xhi(2)
             do i=xlo(1),xhi(1)
                n_fcx(i,j,k,outcomp+comp) = average_values(n_cc(i-1,j,k,incomp+comp), n_cc(i,j,k,incomp+comp),dv)
             end do
          end do
       end do

       ! y-faces
       do k=ylo(3),yhi(3)
          do j=ylo(2),yhi(2)
             do i=ylo(1),yhi(1)
                n_fcy(i,j,k,outcomp+comp) = average_values(n_cc(i,j-1,k,incomp+comp), n_cc(i,j,k,incomp+comp),dv)
             end do
          end do
       end do

       ! z-faces
       do k=zlo(3),zhi(3)
          do j=zlo(2),zhi(2)
             do i=zlo(1),zhi(1)
                n_fcz(i,j,k,outcomp+comp) = average_values(n_cc(i,j,k-1,incomp+comp), n_cc(i,j,k,incomp+comp),dv)
             end do
          end do
       end do

    end do

  end subroutine average_to_faces_3d
  
  function average_values(value1,value2,dv) result(av)
      real(kind=dp_t) :: av
      real(kind=dp_t), intent(in) :: value1, value2
      real(kind=dp_t), intent(in) :: dv
      real(kind=dp_t) :: tmp,tmp1,tmp2
  
      select case(avg_type)
      case(1) ! Arithmetic with a C0-smoothed Heaviside
         if ( (value1 .le. 0.d0) .or. (value2 .le. 0.d0) ) then
            av=0.d0
         else
            tmp1=min(dv*value1,1.d0)
            tmp2=min(dv*value2,1.d0)
            av=(value1+value2)/2.d0*tmp1*tmp2
         end if
      case(2) ! Geometric
         av=sqrt(max(value1,0.d0)*max(value2,0.d0))
      case(3) ! Harmonic
         ! What we want here is the harmonic mean of max(value1,0) and max(value2,0)
         ! Where we define the result to be zero if either one is zero
         ! But numerically we want to avoid here division by zero
         if ( (value1 .le. 10.0d0*tiny(1.d0)) .or. (value2 .le. 10.0d0*tiny(1.d0)) ) then
            av=0.d0
         else
            av=2.d0 / (1.d0/value1 + 1.d0/value2)
         end if
      case(10) ! Arithmetic with (discontinuous) Heaviside
         if ( (value1 .le. 0.d0) .or. (value2 .le. 0.d0) ) then
            av=0.d0
         else
            av=(value1+value2)/2.d0
         end if
      case(11) ! Arithmetic with C1-smoothed Heaviside
         if ( (value1 .le. 0.d0) .or. (value2 .le. 0.d0) ) then
            av=0.d0
         else
            tmp1=dv*value1
            if (tmp1<1.d0) then
               tmp1=(3.d0-2.d0*tmp1)*tmp1**2
            else
               tmp1=1.d0
            end if
            tmp2=dv*value2
            if (tmp2<1.d0) then
               tmp2=(3.d0-2.d0*tmp2)*tmp2**2
            else
               tmp2=1.d0
            end if
            av=(value1+value2)/2.d0*tmp1*tmp2
         endif
      case(12) ! Arithmetic with C2-smoothed Heaviside
         if ( (value1 .le. 0.d0) .or. (value2 .le. 0.d0) ) then
            av=0.d0
         else
            tmp1=dv*value1
            if (tmp1<1.d0) then
               tmp1=(10.d0-15.d0*tmp1+6.d0*tmp1**2)*tmp1**3
            else
               tmp1=1.d0
            end if
            tmp2=dv*value2
            if (tmp2<1.d0) then
               tmp2=(10.d0-15.d0*tmp2+6.d0*tmp2**2)*tmp2**3
            else
               tmp2=1.d0
            end if
            av=(value1+value2)/2.d0*tmp1*tmp2
         endif
      case default
         call bl_error("average_to_faces_2d: invalid avg_type")   
      end select   
  
  end function average_values

end module average_to_faces_module
