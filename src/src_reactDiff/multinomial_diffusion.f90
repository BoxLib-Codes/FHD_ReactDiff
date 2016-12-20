module multinomial_diffusion_module

  use ml_layout_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use BoxLibRNGs
  use bl_rng_module
  use bl_random_module
  use probin_common_module, only: n_cells
  use probin_reactdiff_module, only: nspecies, D_Fick, cross_section, use_bl_rng

  implicit none

  private

  public :: multinomial_diffusion

contains

  ! advances n_old to n_new using multinomial diffusion
  subroutine multinomial_diffusion(mla,n_old,n_new,diff_coef_face, &
                                           dx,dt,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: n_old(:)
    type(multifab) , intent(inout) :: n_new(:)
    type(multifab) , intent(in   ) :: diff_coef_face(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    integer :: n,nlevs,i,dm,spec

    type(bl_prof_timer),save :: bpt

    nlevs = mla%nlevel
    dm = mla%dim
    
    call build(bpt,"multinomial_diffusion")

    ! set new state to zero everywhere, including ghost cells
    do n=1,nlevs
       call multifab_setval(n_new(n),0.d0,all=.true.)
    end do    

    ! copy old state into new in valid region only
    do n=1,nlevs
       call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
    end do    

    ! update with multinomial diffusion, each grid in isolation
    call multinomial_diffusion_update(mla,n_new,diff_coef_face,dx,dt,the_bc_tower)

    ! call sum_boundary to deal with grid boundaries
    do n=1,nlevs
       call multifab_sum_boundary(n_new(n),1)
    end do

    ! properly fill n_new ghost cells
    do n=1,nlevs
       call multifab_fill_boundary(n_new(n))
       call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                            the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
    end do

    call destroy(bpt)

  end subroutine multinomial_diffusion

  subroutine multinomial_diffusion_update(mla,n_new,diff_coef_face,dx,dt,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: n_new(:) ! Old state on input, new state on output in valid region, or increment in ghosts
    type(multifab) , intent(in   ) :: diff_coef_face(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    integer :: n,nlevs,i,dm,ng_n,ng_d
    real(kind=dp_t) :: dv ! Cell volume

    integer      ::  lo(mla%dim),  hi(mla%dim)
    
    real(kind=dp_t), pointer :: np(:,:,:,:)
    real(kind=dp_t), pointer :: dxp(:,:,:,:)
    real(kind=dp_t), pointer :: dyp(:,:,:,:)
    real(kind=dp_t), pointer :: dzp(:,:,:,:)

    nlevs = mla%nlevel
    dm = mla%dim

    ng_n = n_new(1)%ng
    ng_d = diff_coef_face(1,1)%ng

    dv = product(dx(1,1:dm))*cross_section

    ! cannot use OpenMP with tiling since each cell is responsible for updating
    ! cells possibly outside of its file.  OpenMP could be added at the k loop level
    ! with reduction tricks
    do n=1,nlevs
       do i=1,nfabs(n_new(n))
        np  => dataptr(n_new(n),i)
        dxp => dataptr(diff_coef_face(n,1),i)
        dyp => dataptr(diff_coef_face(n,2),i)
        lo = lwb(get_box(n_new(n),i))
        hi = upb(get_box(n_new(n),i))
        select case (dm)
        case (2)
           if(n_cells(2)==1) then ! This is really a 1D domain
              ! Note the in this case the second dimension of dxp has bounds of (0:0)
              call multinomial_diffusion_update_1d(np(:,0,1,:),ng_n, &
                                                   dxp(:,0,1,:),ng_d, &
                                                   lo(1),hi(1),dx(n,1),dt,dv)
           else
              call multinomial_diffusion_update_2d(np(:,:,1,:),ng_n, &
                                                   dxp(:,:,1,:),dyp(:,:,1,:),ng_d, &
                                                   lo,hi,dx(n,:),dt,dv)
           end if
        case (3)
           dzp => dataptr(diff_coef_face(n,3),i)
           call multinomial_diffusion_update_3d(np(:,:,:,:),ng_n, &
                                                dxp(:,:,:,:),dyp(:,:,:,:),dzp(:,:,:,:),ng_d, &
                                                lo,hi,dx(n,:),dt,dv)
        end select
      end do
    end do

  end subroutine multinomial_diffusion_update

  ! For 1D we want to be more efficient by sampling only two binomials instead of 4
  subroutine multinomial_diffusion_update_1d(n_new,ng_n,diffx,ng_d,lo,hi,dx,dt,dv)

    integer        , intent(in   ) :: lo,hi,ng_n,ng_d
    real(kind=dp_t), intent(inout) :: n_new(lo-ng_n:,:) ! Old state on input, new state on output
    real(kind=dp_t), intent(in)    :: diffx(lo-ng_d:,:)
    real(kind=dp_t), intent(in   ) :: dx,dt,dv

    ! local
    integer :: i,j,comp,n_total,n_sum,n_change
    integer, allocatable :: cell_update(:,:) ! Avoid stack overflows and put this on the heap instead
    
    integer, parameter :: n_faces=2
    integer :: fluxes(n_faces) ! Number of particles jumping out of this cell to each of the neighboring cells
    real(kind=dp_t) :: probabilities(n_faces)

    allocate(cell_update(lo-1:hi+1,nspecies))
    cell_update = 0.d0

    do comp=1,nspecies

       do i=lo,hi

          probabilities = (/diffx(i,  comp)*dt/dx**2, &
                            diffx(i+1,comp)*dt/dx**2/)

          if(sum(probabilities)>1.0_dp_t) &
             call bl_error("Explicit CFL stability limit violated for multinomial diffusion")
          if (use_bl_rng) then
             call MultinomialRNG(samples=fluxes, n_samples=n_faces, &
                                 N=max(0, nint(n_new(i,comp)*dv)), p=probabilities, &
                                 engine=rng_eng_diffusion%p)
          else
             call MultinomialRNG(samples=fluxes, n_samples=n_faces, &
                                 N=max(0, nint(n_new(i,comp)*dv)), p=probabilities)
          end if

          ! lo-x face
          cell_update(i  ,comp) = cell_update(i  ,comp) - fluxes(1)
          cell_update(i-1,comp) = cell_update(i-1,comp) + fluxes(1)

          ! hi-x face
          cell_update(i  ,comp) = cell_update(i  ,comp) - fluxes(2)
          cell_update(i+1,comp) = cell_update(i+1,comp) + fluxes(2)

       end do

    end do

    ! increment n_new for all components but remember to convert back to number densities from number of molecules
    n_new(lo-1:hi+1,1:nspecies) = n_new(lo-1:hi+1,1:nspecies) + &
                                        cell_update(lo-1:hi+1,1:nspecies) / dv
         
    deallocate(cell_update)     

  end subroutine multinomial_diffusion_update_1d

  subroutine multinomial_diffusion_update_2d(n_new,ng_n,diffx,diffy,ng_d,lo,hi,dx,dt,dv)

    integer        , intent(in   ) :: lo(:),hi(:),ng_n,ng_d
    real(kind=dp_t), intent(inout) :: n_new(lo(1)-ng_n:,lo(2)-ng_n:,:) ! Old state on input, new state on output
    real(kind=dp_t), intent(in)    :: diffx(lo(1)-ng_d:,lo(2)-ng_d:,:)
    real(kind=dp_t), intent(in)    :: diffy(lo(1)-ng_d:,lo(2)-ng_d:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt,dv

    ! local
    integer :: i,j,comp,n_total,n_sum,n_change
    integer, allocatable :: cell_update(:,:,:) ! Avoid stack overflows and put this on the heap instead
    
    integer, parameter :: n_faces=4
    integer :: fluxes(n_faces) ! Number of particles jumping out of this cell to each of the neighboring cells
    real(kind=dp_t) :: probabilities(n_faces)

    allocate(cell_update(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,nspecies))
    cell_update = 0.d0

    do comp=1,nspecies

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
          
             probabilities = (/diffx(i  ,j,comp)*dt/dx(1)**2, &
                               diffx(i+1,j,comp)*dt/dx(1)**2, &
                               diffy(i,j  ,comp)*dt/dx(2)**2, &
                               diffy(i,j+1,comp)*dt/dx(2)**2/)
                               
             !write(*,*) "species=", comp, " probability=", sum(probabilities), "D=", diffx(i  ,j,comp), " dt=", dt, " dx=", dx(1)
             if(sum(probabilities)>1.0_dp_t) &
                call bl_error("Explicit CFL stability limit violated for multinomial diffusion")

             if (use_bl_rng) then
                call MultinomialRNG(samples=fluxes, n_samples=n_faces, &
                                    N=max(0, nint(n_new(i,j,comp)*dv)), p=probabilities, &
                                    engine=rng_eng_diffusion%p)
             else
                call MultinomialRNG(samples=fluxes, n_samples=n_faces, &
                                    N=max(0, nint(n_new(i,j,comp)*dv)), p=probabilities)
             end if

             ! lo-x face
             cell_update(i  ,j,comp) = cell_update(i  ,j,comp) - fluxes(1)
             cell_update(i-1,j,comp) = cell_update(i-1,j,comp) + fluxes(1)
             
             ! hi-x face
             cell_update(i  ,j,comp) = cell_update(i  ,j,comp) - fluxes(2)
             cell_update(i+1,j,comp) = cell_update(i+1,j,comp) + fluxes(2)
             
             ! lo-y face
             cell_update(i,j  ,comp) = cell_update(i,j  ,comp) - fluxes(3)
             cell_update(i,j-1,comp) = cell_update(i,j-1,comp) + fluxes(3)

             ! hi-y face
             cell_update(i,j  ,comp) = cell_update(i,j  ,comp) - fluxes(4)
             cell_update(i,j+1,comp) = cell_update(i,j+1,comp) + fluxes(4)

          end do
       end do

    end do

    ! increment n_new for all components but remember to convert back to number densities from number of molecules
    n_new(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1:nspecies) = &
                 n_new(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1:nspecies) &
         + cell_update(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1:nspecies) / dv
         
    deallocate(cell_update)     

  end subroutine multinomial_diffusion_update_2d

  subroutine multinomial_diffusion_update_3d(n_new,ng_n,diffx,diffy,diffz,ng_d,lo,hi,dx,dt,dv)

    integer        , intent(in   ) :: lo(:),hi(:),ng_n,ng_d
    real(kind=dp_t), intent(inout) :: n_new(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:,:) ! Old state on input, new state on output
    real(kind=dp_t), intent(in)    :: diffx(lo(1)-ng_d:,lo(2)-ng_d:,lo(3)-ng_d:,:)
    real(kind=dp_t), intent(in)    :: diffy(lo(1)-ng_d:,lo(2)-ng_d:,lo(3)-ng_d:,:)
    real(kind=dp_t), intent(in)    :: diffz(lo(1)-ng_d:,lo(2)-ng_d:,lo(3)-ng_d:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt,dv

    ! local
    integer :: i,j,k,comp,n_total,n_sum,n_change
    integer, allocatable :: cell_update(:,:,:,:) ! Avoid stack overflows and put this on the heap instead
    
    integer, parameter :: n_faces=6
    integer :: fluxes(n_faces) ! Number of particles jumping out of this cell to each of the neighboring cells
    real(kind=dp_t) :: probabilities(n_faces)

    allocate(cell_update(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,nspecies))
    cell_update = 0.d0

    do comp=1,nspecies

       do k=lo(3),hi(3)
       do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          
          probabilities = (/diffx(i  ,j,k,comp)*dt/dx(1)**2, &
                            diffx(i+1,j,k,comp)*dt/dx(1)**2, &
                            diffy(i,j  ,k,comp)*dt/dx(2)**2, &
                            diffy(i,j+1,k,comp)*dt/dx(2)**2, &
                            diffz(i,j,k  ,comp)*dt/dx(3)**2, &
                            diffz(i,j,k+1,comp)*dt/dx(3)**2/)
             
          if(sum(probabilities)>1.0_dp_t) &
             call bl_error("Explicit CFL stability limit violated for multinomial diffusion")

          if (use_bl_rng) then
             call MultinomialRNG(samples=fluxes, n_samples=n_faces, &
                                 N=max(0, nint(n_new(i,j,k,comp)*dv)), p=probabilities, &
                                 engine=rng_eng_diffusion%p)
          else
             call MultinomialRNG(samples=fluxes, n_samples=n_faces, &
                                 N=max(0, nint(n_new(i,j,k,comp)*dv)), p=probabilities)
          end if

          ! lo-x face
          cell_update(i  ,j,k,comp) = cell_update(i  ,j,k,comp) - fluxes(1)
          cell_update(i-1,j,k,comp) = cell_update(i-1,j,k,comp) + fluxes(1)
          
          ! hi-x face
          cell_update(i  ,j,k,comp) = cell_update(i  ,j,k,comp) - fluxes(2)
          cell_update(i+1,j,k,comp) = cell_update(i+1,j,k,comp) + fluxes(2)

          ! lo-y face
          cell_update(i,j  ,k,comp) = cell_update(i,j  ,k,comp) - fluxes(3)
          cell_update(i,j-1,k,comp) = cell_update(i,j-1,k,comp) + fluxes(3)

          ! hi-y face
          cell_update(i,j  ,k,comp) = cell_update(i,j  ,k,comp) - fluxes(4)
          cell_update(i,j+1,k,comp) = cell_update(i,j+1,k,comp) + fluxes(4)

          ! lo-z face
          cell_update(i,j,k  ,comp) = cell_update(i,j,k  ,comp) - fluxes(5)
          cell_update(i,j,k-1,comp) = cell_update(i,j,k-1,comp) + fluxes(5)

          ! hi-z face
          cell_update(i,j,k  ,comp) = cell_update(i,j,k  ,comp) - fluxes(6)
          cell_update(i,j,k+1,comp) = cell_update(i,j,k+1,comp) + fluxes(6)

       end do
       end do
       end do

    end do

    ! increment n_new for all components but remember to convert back to number densities from number of molecules
    n_new(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1:nspecies) = &
                 n_new(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1:nspecies) &
         + cell_update(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1:nspecies) / dv
         
    deallocate(cell_update)     

  end subroutine multinomial_diffusion_update_3d

end module multinomial_diffusion_module
