module stochastic_n_fluxdiv_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  use multifab_fill_random_module
  use multifab_physbc_module
  use average_to_faces_module
  use div_and_grad_module
  use bl_rng_module
  use probin_common_module, only: variance_coef_mass, initial_variance, density_weights
  use probin_reactdiff_module, only: nspecies, cross_section, use_bl_rng

  implicit none

  private

  public :: stochastic_n_fluxdiv, fill_mass_stochastic, &
       init_mass_stochastic, destroy_mass_stochastic, &
       add_init_n_fluctuations

  ! stochastic fluxes for mass densities are face-centered
  type(multifab), allocatable, save :: stoch_W_fc(:,:,:)

  integer, save :: n_rngs ! how many random number stages
  
contains
  
  ! This computes div (sqrt(2*variance*D*n_cc / (dt*dV)) * W ) where W is a collection of i.i.d. standard normal variates
  ! It can easily be generalized to compute
  !  div (sqrt(2*variance*D*n_cc / (dt*dV)) * (alpha*W_1+beta*W_2) ) if n_rng=2
  !  but this is not used or implemented at present
  subroutine stochastic_n_fluxdiv(mla,n_cc,diff_coef_face,stoch_fluxdiv,dx,dt, &
                                  the_bc_tower,increment_in)

    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(in   )  :: n_cc(:)
    type(multifab) , intent(in   )  :: diff_coef_face(:,:)
    type(multifab) , intent(inout)  :: stoch_fluxdiv(:)
    real(kind=dp_t), intent(in   )  :: dx(:,:)
    real(kind=dp_t), intent(in   )  :: dt
    type(bc_tower) , intent(in   )  :: the_bc_tower
    logical  , intent(in), optional :: increment_in ! Increment or overwrite stoch_fluxdiv argument?

    integer :: i,dm,n,nlevs
    real(kind=dp_t)  :: variance, dv

    type(multifab) :: flux(mla%nlevel,mla%dim)

    logical :: increment_div

    type(bl_prof_timer), save :: bpt

    call build(bpt,"stochastic_n_fluxdiv")

    dm = mla%dim
    nlevs = mla%nlevel

    dv = product(dx(1,1:dm))*cross_section

    increment_div = .false.
    if (present(increment_in)) increment_div = increment_in

    ! single cell case set stochastic mass fluxdiv to zero 
    ! (or its increment if increment_in=T) and return
    if ((multifab_volume(n_cc(1))/nspecies)<=1) then
       if (.not. increment_div) then
          do n=1,nlevs
             call multifab_setval(stoch_fluxdiv(n),0.d0,all=.true.)
          end do
       end if
       return
    end if

    do n=1,nlevs
       do i=1,dm
          call multifab_build_edge(flux(n,i),mla%la(n),nspecies,0,i)
       end do
    end do

    ! average n to faces, store in "flux" so as to avoid an extra multifab
    ! alternatively one could multiply diff_coeff_face with this number
    call average_to_faces(mla,n_cc,flux,1,1,nspecies,dv)

    ! assumble fluxes on faces, sqrt(2*D_k*n_k / (dt*dV)) * random_normal
    variance = sqrt(2.d0*variance_coef_mass/(dv*dt))        
    call assemble_stoch_n_fluxes(mla,n_cc,diff_coef_face,flux)
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_s(flux(n,i), variance, 0)
       end do
    end do      

    ! take flux divergence
    call compute_div(mla,flux,stoch_fluxdiv,dx,1,1,nspecies, &
                     increment_in=increment_div)

    do n=1,nlevs
       do i=1,dm
          call multifab_destroy(flux(n,i))
       end do
    end do

    call destroy(bpt)

  end subroutine stochastic_n_fluxdiv

  subroutine assemble_stoch_n_fluxes(mla,n_cc,diff_coef_face,flux)

    ! note: n averaged to faces is stored in "flux" on entry to this subroutine

    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(in   )  :: n_cc(:)
    type(multifab) , intent(in   )  :: diff_coef_face(:,:)
    type(multifab) , intent(inout)  :: flux(:,:) ! On input this contains the face-averaged number densities

    integer :: i,dm,n,nlevs,lo(mla%dim),hi(mla%dim)
    integer :: ng_n,ng_d,ng_f,ng_s

    real(kind=dp_t), pointer :: np(:,:,:,:)
    real(kind=dp_t), pointer :: dx(:,:,:,:)
    real(kind=dp_t), pointer :: dy(:,:,:,:)
    real(kind=dp_t), pointer :: dz(:,:,:,:)
    real(kind=dp_t), pointer :: fx(:,:,:,:)
    real(kind=dp_t), pointer :: fy(:,:,:,:)
    real(kind=dp_t), pointer :: fz(:,:,:,:)
    real(kind=dp_t), pointer :: sx(:,:,:,:)
    real(kind=dp_t), pointer :: sy(:,:,:,:)
    real(kind=dp_t), pointer :: sz(:,:,:,:)

    type(mfiter) :: mfi
    type(box) :: xnodalbox, ynodalbox, znodalbox
    integer :: xlo(mla%dim), xhi(mla%dim)
    integer :: ylo(mla%dim), yhi(mla%dim)
    integer :: zlo(mla%dim), zhi(mla%dim)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"assemble_stoch_n_fluxes")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_n = n_cc(1)%ng
    ng_d = diff_coef_face(1,1)%ng
    ng_f = flux(1,1)%ng
    ng_s = stoch_W_fc(1,1,1)%ng
    
    !$omp parallel private(mfi,n,i,xnodalbox,ynodalbox,znodalbox,xlo,ylo,zlo) &
    !$omp private(xhi,yhi,zhi,np,dx,dy,dz,fx,fy,fz,sx,sy,sz,lo,hi)
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

          np => dataptr(n_cc(n),i)
          dx => dataptr(diff_coef_face(n,1),i)
          dy => dataptr(diff_coef_face(n,2),i)
          fx => dataptr(flux(n,1),i)
          fy => dataptr(flux(n,2),i)
          sx => dataptr(stoch_W_fc(n,1,1),i)
          sy => dataptr(stoch_W_fc(n,2,1),i)
          lo = lwb(get_box(n_cc(n),i))
          hi = upb(get_box(n_cc(n),i))
          select case (dm)
          case (2)
             call assemble_stoch_n_fluxes_2d(np(:,:,1,:),ng_n, &
                                             dx(:,:,1,:),dy(:,:,1,:),ng_d, &
                                             fx(:,:,1,:),fy(:,:,1,:),ng_f, &
                                             sx(:,:,1,:),sy(:,:,1,:),ng_s, lo,hi, &
                                             xlo,xhi,ylo,yhi)
          case (3)
             dz => dataptr(diff_coef_face(n,3),i)
             fz => dataptr(flux(n,3),i)
             sz => dataptr(stoch_W_fc(n,3,1),i)
             call assemble_stoch_n_fluxes_3d(np(:,:,:,:),ng_n, &
                                             dx(:,:,:,:),dy(:,:,:,:),dz(:,:,:,:),ng_d, &
                                             fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:),ng_f, &
                                             sx(:,:,:,:),sy(:,:,:,:),sz(:,:,:,:),ng_s, lo,hi, &
                                             xlo,xhi,ylo,yhi,zlo,zhi)
          end select
       end do
    end do
    !$omp end parallel

    ! sync the fluxes at the boundaries
    do n=1,nlevs
       do i=1,dm
          call multifab_internal_sync(flux(n,i))
          call multifab_fill_boundary(flux(n,i))  
       end do
    end do

    call destroy(bpt)

  end subroutine assemble_stoch_n_fluxes

  subroutine assemble_stoch_n_fluxes_2d(n_cc,ng_n,coefx,coefy,ng_d,fluxx,fluxy,ng_f, &
                                        stochx,stochy,ng_s,glo,ghi, &
                                        xlo,xhi,ylo,yhi)

    integer        , intent(in   ) :: glo(:),ghi(:),ng_n,ng_d,ng_f,ng_s
    integer        , intent(in   ) :: xlo(:),xhi(:),ylo(:),yhi(:)
    real(kind=dp_t), intent(in   ) ::   n_cc(glo(1)-ng_n:,glo(2)-ng_n:,:)
    real(kind=dp_t), intent(in   ) ::  coefx(glo(1)-ng_d:,glo(2)-ng_d:,:)
    real(kind=dp_t), intent(in   ) ::  coefy(glo(1)-ng_d:,glo(2)-ng_d:,:)
    real(kind=dp_t), intent(inout) ::  fluxx(glo(1)-ng_f:,glo(2)-ng_f:,:)
    real(kind=dp_t), intent(inout) ::  fluxy(glo(1)-ng_f:,glo(2)-ng_f:,:)
    real(kind=dp_t), intent(in   ) :: stochx(glo(1)-ng_s:,glo(2)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: stochy(glo(1)-ng_s:,glo(2)-ng_s:,:)

    integer :: i,j

    ! x-fluxes
    do j=xlo(2),xhi(2)
       do i=xlo(1),xhi(1)
          fluxx(i,j,1:nspecies) = &
               sqrt(coefx(i,j,1:nspecies)*fluxx(i,j,1:nspecies))*stochx(i,j,1:nspecies)
       end do
    end do

    ! y-fluxes
    do j=ylo(2),yhi(2)
       do i=ylo(1),yhi(1)
          fluxy(i,j,1:nspecies) = &
               sqrt(coefy(i,j,1:nspecies)*fluxy(i,j,1:nspecies))*stochy(i,j,1:nspecies)
       end do
    end do

  end subroutine assemble_stoch_n_fluxes_2d

  subroutine assemble_stoch_n_fluxes_3d(n_cc,ng_n,coefx,coefy,coefz,ng_d,fluxx,fluxy,fluxz,ng_f, &
                                        stochx,stochy,stochz,ng_s,glo,ghi, &
                                        xlo,xhi,ylo,yhi,zlo,zhi)

    integer        , intent(in   ) :: glo(:),ghi(:),ng_n,ng_d,ng_f,ng_s
    integer        , intent(in   ) :: xlo(:),xhi(:),ylo(:),yhi(:),zlo(:),zhi(:)
    real(kind=dp_t), intent(in   ) ::   n_cc(glo(1)-ng_n:,glo(2)-ng_n:,glo(3)-ng_n:,:)
    real(kind=dp_t), intent(in   ) ::  coefx(glo(1)-ng_d:,glo(2)-ng_d:,glo(3)-ng_d:,:)
    real(kind=dp_t), intent(in   ) ::  coefy(glo(1)-ng_d:,glo(2)-ng_d:,glo(3)-ng_d:,:)
    real(kind=dp_t), intent(in   ) ::  coefz(glo(1)-ng_d:,glo(2)-ng_d:,glo(3)-ng_d:,:)
    real(kind=dp_t), intent(inout) ::  fluxx(glo(1)-ng_f:,glo(2)-ng_f:,glo(3)-ng_f:,:)
    real(kind=dp_t), intent(inout) ::  fluxy(glo(1)-ng_f:,glo(2)-ng_f:,glo(3)-ng_f:,:)
    real(kind=dp_t), intent(inout) ::  fluxz(glo(1)-ng_f:,glo(2)-ng_f:,glo(3)-ng_f:,:)
    real(kind=dp_t), intent(in   ) :: stochx(glo(1)-ng_s:,glo(2)-ng_s:,glo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: stochy(glo(1)-ng_s:,glo(2)-ng_s:,glo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: stochz(glo(1)-ng_s:,glo(2)-ng_s:,glo(3)-ng_s:,:)

    integer :: i,j,k

    ! x-fluxes
    do k=xlo(3),xhi(3)
       do j=xlo(2),xhi(2)
          do i=xlo(1),xhi(1)
             fluxx(i,j,k,1:nspecies) = &
                  sqrt(coefx(i,j,k,1:nspecies)*fluxx(i,j,k,1:nspecies))*stochx(i,j,k,1:nspecies)
          end do
       end do
    end do

    ! y-fluxes
    do k=ylo(3),yhi(3)
       do j=ylo(2),yhi(2)
          do i=ylo(1),yhi(1)
             fluxy(i,j,k,1:nspecies) = &
                  sqrt(coefy(i,j,k,1:nspecies)*fluxy(i,j,k,1:nspecies))*stochy(i,j,k,1:nspecies)
          end do
       end do
    end do

    ! z-fluxes
    do k=zlo(3),zhi(3)
       do j=zlo(2),zhi(2)
          do i=zlo(1),zhi(1)
             fluxz(i,j,k,1:nspecies) = &
                  sqrt(coefz(i,j,k,1:nspecies)*fluxz(i,j,k,1:nspecies))*stochz(i,j,k,1:nspecies)
          end do
       end do
    end do

  end subroutine assemble_stoch_n_fluxes_3d

  ! call this once at the beginning of simulation to allocate multifabs
  ! that will hold random numbers
  subroutine init_mass_stochastic(mla,n_rngs_in)

    type(ml_layout), intent(in   ) :: mla
    integer        , intent(in   ) :: n_rngs_in

    ! local
    integer :: n,nlevs,i,dm,rng
    
    type(bl_prof_timer), save :: bpt

    call build(bpt,"init_mass_stochastic")

    n_rngs = n_rngs_in

    nlevs = mla%nlevel
    dm = mla%dim

    allocate(stoch_W_fc(mla%nlevel, mla%dim, n_rngs))

    do n=1,nlevs
       do rng=1,n_rngs
          do i=1,dm
             ! we need one face-centered flux for each concentration
             call multifab_build_edge(stoch_W_fc(n,i,rng),mla%la(n),nspecies,0,i)
          end do
       end do ! end loop over n_rngs
    end do ! end loop over nlevs

    call destroy(bpt)

  end subroutine init_mass_stochastic

  ! call this once at the end of simulation to deallocate memory
  subroutine destroy_mass_stochastic(mla)

    type(ml_layout), intent(in   ) :: mla

    ! local
    integer :: n,nlevs,i,dm,rng
    
    type(bl_prof_timer), save :: bpt

    call build(bpt,"destroy_mass_stochastic")

    nlevs = mla%nlevel
    dm = mla%dim

    do n=1,nlevs
       do rng=1,n_rngs
          do i=1,dm
             call multifab_destroy(stoch_W_fc(n,i,rng))
          end do
       end do
    end do
    
    deallocate(stoch_W_fc)

    call destroy(bpt)

  end subroutine destroy_mass_stochastic

  subroutine fill_mass_stochastic(mla,the_bc_level)
  
    type(ml_layout), intent(in   )  :: mla
    type(bc_level) , intent(in   )  :: the_bc_level(:)

    ! Local variables
    integer :: dm,nlevs,box,i,rng
    
    type(bl_prof_timer), save :: bpt

    call build(bpt,"fill_mass_stochastic")

    nlevs = mla%nlevel
    dm    = mla%dim    
    
    ! generate and store the stochastic flux (random numbers)
    if (use_bl_rng) then
       do rng=1, n_rngs
          do i = 1,dm
             call multifab_fill_random(stoch_W_fc(:,i,rng), rng_eng=rng_eng_diffusion)
          end do
       end do
    else
       do rng=1, n_rngs
          do i = 1,dm
             call multifab_fill_random(stoch_W_fc(:,i,rng))
          end do
       end do
    end if

    ! apply boundary conditions to stochastic fluxes
    call stoch_mass_bc(mla,the_bc_level)
  
    call destroy(bpt)

  end subroutine fill_mass_stochastic

  subroutine stoch_mass_bc(mla,the_bc_level)
    
    type(ml_layout), intent(in   )  :: mla
    type(bc_level) , intent(in   )  :: the_bc_level(:)

    ! local
    integer :: n,nlevs,dm,idim,rng,i,ng_f
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: fp(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"stoch_mass_bc")

    nlevs = mla%nlevel
    dm = mla%dim
    
    ng_f = stoch_W_fc(1,1,1)%ng

    do n=1,nlevs
       do idim=1,dm
          do rng=1,n_rngs
             do i=1,nfabs(stoch_W_fc(n,idim,rng))
                fp => dataptr(stoch_W_fc(n,idim,rng),i)
                lo = lwb(get_box(stoch_W_fc(n,idim,rng),i))
                hi = upb(get_box(stoch_W_fc(n,idim,rng),i))
                select case (dm)
                case (2)
                   call stoch_mass_bc_2d(fp(:,:,1,:),ng_f,idim,lo,hi, &
                                         the_bc_level(n)%phys_bc_level_array(i,:,:))
                case (3)
                   call stoch_mass_bc_3d(fp(:,:,:,:),ng_f,idim,lo,hi, &
                                         the_bc_level(n)%phys_bc_level_array(i,:,:))
                end select
             end do
          end do
       end do
    end do

    call destroy(bpt)

  end subroutine stoch_mass_bc

  subroutine stoch_mass_bc_2d(sflux,ng_f,idim,lo,hi,phys_bc)

      integer        , intent(in   ) :: lo(:),hi(:),ng_f,idim
      integer        , intent(in   ) :: phys_bc(:,:)
      real(kind=dp_t), intent(inout) :: sflux(lo(1)-ng_f:,lo(2)-ng_f:,:)

      if (idim .eq. 1) then

         if (phys_bc(1,1) .eq. NO_SLIP_WALL .or. phys_bc(1,1) .eq. SLIP_WALL) then
            sflux(lo(1),lo(2):hi(2),:) = 0.d0
         end if

         if (phys_bc(1,2) .eq. NO_SLIP_WALL .or. phys_bc(1,2) .eq. SLIP_WALL) then
            sflux(hi(1)+1,lo(2):hi(2),:) = 0.d0
         end if

         if (phys_bc(1,1) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1),lo(2):hi(2),:) = sqrt(2.0d0)*sflux(lo(1),lo(2):hi(2),:)
         end if

         if (phys_bc(1,2) .eq. NO_SLIP_RESERVOIR) then
            sflux(hi(1)+1,lo(2):hi(2),:) = sqrt(2.0d0)*sflux(hi(1)+1,lo(2):hi(2),:)
         end if

      else

         if (phys_bc(2,1) .eq. NO_SLIP_WALL .or. phys_bc(2,1) .eq. SLIP_WALL) then
            sflux(lo(1):hi(1),lo(2),:) = 0.d0
         end if

         if (phys_bc(2,2) .eq. NO_SLIP_WALL .or. phys_bc(2,2) .eq. SLIP_WALL) then
            sflux(lo(1):hi(1),hi(2)+1,:) = 0.d0
         end if

         if (phys_bc(2,1) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1):hi(1),lo(2),:) = sqrt(2.0d0)*sflux(lo(1):hi(1),lo(2),:)
         end if

         if (phys_bc(2,2) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1):hi(1),hi(2)+1,:) = sqrt(2.0d0)*sflux(lo(1):hi(1),hi(2)+1,:)
         end if

      end if

  end subroutine stoch_mass_bc_2d

  subroutine stoch_mass_bc_3d(sflux,ng_f,idim,lo,hi,phys_bc)

      integer        , intent(in   ) :: lo(:),hi(:),ng_f,idim
      integer        , intent(in   ) :: phys_bc(:,:)
      real(kind=dp_t), intent(inout) :: sflux(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:,:)

      if (idim .eq. 1) then

         if (phys_bc(1,1) .eq. NO_SLIP_WALL .or. phys_bc(1,1) .eq. SLIP_WALL) then
            sflux(lo(1),lo(2):hi(2),lo(3):hi(3),:) = 0.d0
         end if

         if (phys_bc(1,2) .eq. NO_SLIP_WALL .or. phys_bc(1,2) .eq. SLIP_WALL) then
            sflux(hi(1)+1,lo(2):hi(2),lo(3):hi(3),:) = 0.d0
         end if

         if (phys_bc(1,1) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1),lo(2):hi(2),lo(3):hi(3),:) = sqrt(2.0d0)*sflux(lo(1),lo(2):hi(2),lo(3):hi(3),:)
         end if

         if (phys_bc(1,2) .eq. NO_SLIP_RESERVOIR) then
            sflux(hi(1)+1,lo(2):hi(2),lo(3):hi(3),:) = sqrt(2.0d0)*sflux(hi(1)+1,lo(2):hi(2),lo(3):hi(3),:)
         end if

      else if (idim .eq. 2) then

         if (phys_bc(2,1) .eq. NO_SLIP_WALL .or. phys_bc(2,1) .eq. SLIP_WALL) then
            sflux(lo(1):hi(1),lo(2),lo(3):hi(3),:) = 0.d0
         end if

         if (phys_bc(2,2) .eq. NO_SLIP_WALL .or. phys_bc(2,2) .eq. SLIP_WALL) then
            sflux(lo(1):hi(1),hi(2)+1,lo(3):hi(3),:) = 0.d0
         end if

         if (phys_bc(2,1) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1):hi(1),lo(2),lo(3):hi(3),:) = sqrt(2.0d0)*sflux(lo(1):hi(1),lo(2),lo(3):hi(3),:)
         end if

         if (phys_bc(2,2) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1):hi(1),hi(2)+1,lo(3):hi(3),:) = sqrt(2.0d0)*sflux(lo(1):hi(1),hi(2)+1,lo(3):hi(3),:)
         end if

      else

         if (phys_bc(3,1) .eq. NO_SLIP_WALL .or. phys_bc(3,1) .eq. SLIP_WALL) then
            sflux(lo(1):hi(1),lo(2):hi(2),lo(3),:) = 0.d0
         end if

         if (phys_bc(3,2) .eq. NO_SLIP_WALL .or. phys_bc(3,2) .eq. SLIP_WALL) then
            sflux(lo(1):hi(1),lo(2):hi(2),hi(3)+1,:) = 0.d0
         end if

         if (phys_bc(3,1) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1):hi(1),lo(2):hi(2),lo(3),:) = sqrt(2.0d0)*sflux(lo(1):hi(1),lo(2):hi(2),lo(3),:)
         end if

         if (phys_bc(3,2) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1):hi(1),lo(2):hi(2),hi(3)+1,:) = sqrt(2.0d0)*sflux(lo(1):hi(1),lo(2):hi(2),hi(3)+1,:)
         end if

      end if

  end subroutine stoch_mass_bc_3d

  subroutine add_init_n_fluctuations(mla,n_init,dx,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: n_init(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    integer :: n,nlevs,dm,spec,n_cell
    real(kind=dp_t) :: dn_sum, dv

    type(multifab) :: n_temp(mla%nlevel)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"add_init_n_fluctuations")

    nlevs = mla%nlevel
    dm = mla%dim

    dv = product(dx(1,1:dm))*cross_section

    ! the number of ghost cells must match variance_mfab input to multifab_fill_random
    ! the values in the ghost cells do not have to be added to n_init since we
    ! fill ghost cells for n_init afterwards
    do n=1,nlevs
       call multifab_build(n_temp(n),mla%la(n),nspecies,n_init(n)%ng)
    end do

    n_cell = multifab_volume(n_temp(1)) / nspecies

    ! create a multifab full of random numbers
    do n=1,nlevs
       
       if (use_bl_rng) then
          call multifab_fill_random(n_temp(n:n), &
                                    variance_mfab=n_init, &
                                    variance=abs(initial_variance)/dv, &  ! We do not multiply here by variance_coef_mass
                                    rng_eng=rng_eng_init)
       else
          call multifab_fill_random(n_temp(n:n), &
                                    variance_mfab=n_init, &
                                    variance=abs(initial_variance)/dv)  ! We do not multiply here by variance_coef_mass

       end if
  
       if(initial_variance<0.0d0) then
          ! Make sure this sums to zero
          do spec=1, nspecies
             dn_sum = multifab_sum_c(n_temp(n),spec,1) / dble(n_cell)
             call multifab_sub_sub_s_c(n_temp(n),spec,dn_sum,1,0)
          end do
       end if   
          
    end do

    do n=1,nlevs
       if(.false..and.sum(abs(density_weights))>0.0d0) then ! For A+B<->C tests
          ! Generate only one random perturbation and weight it: dn_k = w_k * dn_0
          do spec=1,nspecies
             call multifab_saxpy_3_cc(n_init(n),spec,density_weights(spec),n_temp(n),1,1)             
          end do
       else   
          call multifab_plus_plus_c(n_init(n),1,n_temp(n),1,nspecies,0)
       end if   
       call multifab_fill_boundary(n_init(n))
       call multifab_physbc(n_init(n),1,scal_bc_comp,nspecies, &
                            the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
    end do

    do n=1,nlevs
       call multifab_destroy(n_temp(n))
    end do    

    call destroy(bpt)

  end subroutine add_init_n_fluctuations
   
end module stochastic_n_fluxdiv_module
