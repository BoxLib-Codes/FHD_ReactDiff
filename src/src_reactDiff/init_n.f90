module init_n_module

  use bl_types
  use bl_constants_module
  use ml_layout_module
  use multifab_physbc_module
  use define_bc_module
  use bc_module
  use bl_rng_module
  use bl_random_module
  use BoxLibRNGs


  use probin_common_module, only: prob_lo, prob_hi, prob_type, initial_variance, &
                                  perturb_width, smoothing_width
  use probin_reactdiff_module, only: nspecies, n_init_in, model_file_init, &
                                     cross_section, integer_populations, use_bl_rng
  
  implicit none

  private

  public :: init_n, init_n_model

  ! prob_type codes for LowMach:
  ! 0=thermodynamic equilibrium, n=n_init_in(1,1:nspecies)
  ! 1=gaussian spreading (order of accuracy testing)
  ! 2=gradient along y, n=n_init_in(1,1:nspecies) on bottom (y=0) and n_init(2,1:nspecies) on top (y=Ly)
  ! 3=1+sin^2(pi*x)*sin^2(pi*y)*sin^2(pi*z) test problem
  ! 4=vertical stripe (thirds of domain)

contains

  subroutine init_n(mla,n_init,dx,the_bc_tower)

    ! initialize rho_i and umac in the valid region
    ! we first initialize c_i in the valid region
    ! then enforce that sum(c_i)=1 by overwriting the final concentration,
    ! and then use the EOS to compute rho_i

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: n_init(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower ), intent(in   ) :: the_bc_tower
 
    ! local variables
    integer                        :: lo(mla%dim), hi(mla%dim)
    integer                        :: i, dm, n, nlevs, ng_n
    real(kind=dp_t), pointer       :: np(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "init_n")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_n = n_init(1)%ng

    ! looping over boxes 
    do n=1,nlevs
       do i=1,nfabs(n_init(n))
          np => dataptr(n_init(n),i)
          lo = lwb(get_box(n_init(n),i))
          hi = upb(get_box(n_init(n),i))
          select case (dm)
          case (2)
             call init_n_2d(np(:,:,1,:),ng_n,lo,hi,dx(n,:))
          case (3)
             call init_n_3d(np(:,:,:,:),ng_n,lo,hi,dx(n,:))
          end select          
       end do
    end do

    call destroy(bpt)

   do n=1,nlevs
      call multifab_fill_boundary(n_init(n))
      call multifab_physbc(n_init(n),1,scal_bc_comp,nspecies, &
                           the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
   end do


  end subroutine init_n

  subroutine init_n_2d(n_init,ng_n,lo,hi,dx)

    integer          :: lo(:), hi(:), ng_n
    real(kind=dp_t)  :: n_init(lo(1)-ng_n:,lo(2)-ng_n:,:)
    real(kind=dp_t)  :: dx(:)
 
    ! local varables
    integer         :: i,j,species
    real(kind=dp_t) :: x,y,r,cen(2),L(2)
    real(kind=dp_t) :: one_fraction_domain1,one_fraction_domain2,x1,x2,stripe_ratio,rad

    L(1:2) = prob_hi(1:2)-prob_lo(1:2) ! Domain length
    
    select case (prob_type)

    case(0) 
       !============================================================
       ! Thermodynamic equilibrium
       !============================================================

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             n_init(i,j,1:nspecies) = n_init_in(1,1:nspecies)

          end do
       end do

    case(1) 
       !=============================================================
       ! Initializing from a Gaussian
       !=============================================================

       cen(1:2) = 0.6d0*prob_lo(1:2) + 0.4d0*prob_hi(1:2)

       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)

             r = sqrt((x-cen(1))**2 + (y-cen(2))**2)

             n_init(i,j,1:nspecies) = n_init_in(1,1:nspecies)*exp(-100.d0*r**2)

          end do
       end do

    case(2) 
       !=========================================================
       ! Initializing with constant gradient 
       !=========================================================

       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2) 
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1) 

             ! linear gradient in rho
             n_init(i,j,1:nspecies) = n_init_in(1,1:nspecies) + & 
                  (n_init_in(2,1:nspecies) - n_init_in(1,1:nspecies))*(y-prob_lo(2))/L(2)

          end do
       end do

    case(3) 
       !=========================================================
       ! 1+sin^2(pi*x)*sin^2(pi*y)*sin^2(pi*z) test problem
       !=========================================================

       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2) 
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1) 

             n_init(i,j,1:nspecies) = 1.d0 + sin(M_PI*x)**2 * sin(M_PI*y)**2

          end do
       end do

    case(4)
       !=================================================================
       ! vertical stripe having width = perturb_width*dx(1) 
       ! n_init = n_init_in(1,:) inside, n_init = n_init_in (2,:) outside
       ! can be discontinous or smooth depending on smoothing_width
       !=================================================================

       stripe_ratio = perturb_width*dx(1)/L(1)

       ! prob_lo(1) < one_fraction_domain2 < center < one_fraction_domain1 < prob_hi(1) 
       one_fraction_domain1=0.5d0*(1.0d0-stripe_ratio)*prob_lo(1)+0.5d0*(1.0d0+stripe_ratio)*prob_hi(1)
       one_fraction_domain2=0.5d0*(1.0d0+stripe_ratio)*prob_lo(1)+0.5d0*(1.0d0-stripe_ratio)*prob_hi(1)

       do i=lo(1),hi(1)
          x1 =(prob_lo(1) + dx(1)*(dble(i)+0.5d0) - one_fraction_domain1)
          x2 =(prob_lo(1) + dx(1)*(dble(i)+0.5d0) - one_fraction_domain2)
       
          if (smoothing_width > 0.d0) then
             ! smooth interface
             do j=lo(2),hi(2)
                n_init(i,j,1:nspecies) = n_init_in(2,1:nspecies) + &
                     0.5d0*(n_init_in(2,1:nspecies)-n_init_in(1,1:nspecies)) * &
                     (tanh(x1/(smoothing_width*dx(1))) - tanh(x2/(smoothing_width*dx(1))))
             end do
          else
             ! discontinuous interface 
             if (x2 < 0.d0 .or. x1 > 0.d0) then   ! outside the stripe
                do j=lo(2),hi(2)
                   n_init(i,j,1:nspecies) = n_init_in(2,1:nspecies)
                end do
             else                                 ! inside the stripe
                do j=lo(2),hi(2)
                   n_init(i,j,1:nspecies) = n_init_in(1,1:nspecies)
                end do
             end if
          end if
       end do

    case(5)
       !=================================================================
       ! bubble having radius = 0.5*perturb_width*dx(1) 
       ! n_init = n_init_in(1,:) inside, n_init = n_init_in (2,:) outside
       ! can be discontinous or smooth depending on smoothing_width
       !=================================================================

       rad = 0.5d0*perturb_width*dx(1) 

       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2) - 0.5d0*(prob_lo(2)+prob_hi(2))
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1) - 0.5d0*(prob_lo(1)+prob_hi(1))

             r = sqrt(x**2 + y**2)

             if (smoothing_width .eq. 0) then
                ! discontinuous interface
                if (r .lt. rad) then
                   n_init(i,j,1:nspecies) = n_init_in(1,1:nspecies)
                else
                   n_init(i,j,1:nspecies) = n_init_in(2,1:nspecies)
                end if
             else
                ! smooth interface
                n_init(i,j,1:nspecies) = n_init_in(1,1:nspecies) + &
                     (n_init_in(2,1:nspecies) - n_init_in(1,1:nspecies))* &
                     0.5d0*(1.d0 + tanh((r-rad)/(smoothing_width*dx(1))))
             end if

          end do
       end do

    case default
       
       call bl_error("init_n_2d: prob_type not supported")

    end select

    if(integer_populations) then ! Ensure that the initial number of molecules are integers
       if(initial_variance<0.0d0) then ! Distribute the particles on the box using a multinomial sampler
       
          ! To do this really correctly one should do this for the whole system, not box per box
          if (parallel_IOProcessor()) write(*,*) "Using multinomial initial distribution PER BOX, not domain"
          
          do species=1, nspecies
            call sample_integers(n_init(lo(1):hi(1),lo(2):hi(2),species), &
                     ncells=size(n_init(lo(1):hi(1),lo(2):hi(2),species)), &
                     dv=dx(1)*dx(2)*cross_section)
          end do
       
       else ! Make the number of molecules in each cell Poisson distributed with desired mean
       
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                call round_to_integers(n_init(i,j,1:nspecies), dv=dx(1)*dx(2)*cross_section)
             end do
          end do    
           
       end if    
    end if
    
  end subroutine init_n_2d

  subroutine init_n_3d(n_init,ng_n,lo,hi,dx)

    integer          :: lo(:), hi(:), ng_n
    real(kind=dp_t)  :: n_init(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:,:)
    real(kind=dp_t)  :: dx(:)
 
    ! local varables
    integer         :: i,j,k,species
    real(kind=dp_t) :: x,y,z,r,cen(3),L(3)
    real(kind=dp_t) :: one_fraction_domain1,one_fraction_domain2,x1,x2,stripe_ratio,rad

    L(1:3) = prob_hi(1:3)-prob_lo(1:3) ! Domain length
    
    select case (prob_type)

    case(0) 
       !================================================================================
       ! Thermodynamic equilibrium
       !================================================================================

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)

                n_init(i,j,k,1:nspecies) = n_init_in(1,1:nspecies)

             end do
          end do
       end do

    case(1) 
       !================================================================================
       ! Initializing from a Gaussian
       !================================================================================
       cen(1:3) = 0.6d0*prob_lo(1:3) + 0.4d0*prob_hi(1:3)

       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             do i=lo(1),hi(1)
                x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)

                r = sqrt((x-cen(1))**2 + (y-cen(2))**2 + (z-cen(3))**2)

                n_init(i,j,k,1:nspecies) = n_init_in(1,1:nspecies)*exp(-100.d0*r**2)

             end do
          end do
       end do

    case(2) 
       !========================================================
       ! Initializing with constant gradient
       !========================================================

       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3) 
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             do i=lo(1),hi(1)
                x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)

                n_init(i,j,k,1:nspecies) = n_init_in(1,1:nspecies) + &
                     (n_init_in(2,1:nspecies) - n_init_in(1,1:nspecies))*(y-prob_lo(2))/L(2)

             end do
          end do
       end do

    case(3) 
       !========================================================
       ! 1+sin^2(pi*x)*sin^2(pi*y)*sin^2(pi*z) test problem
       !========================================================

       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3) 
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             do i=lo(1),hi(1)
                x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)

                n_init(i,j,k,1:nspecies) = 1.d0 + sin(M_PI*x)**2 * sin(M_PI*y)**2 * sin(M_PI*z)**2

             end do
          end do
       end do

    case(4)
       !=================================================================
       ! vertical stripe having width = perturb_width*dx(1) 
       ! n_init = n_init_in(1,:) inside, n_init = n_init_in (2,:) outside
       ! can be discontinous or smooth depending on smoothing_width
       !=================================================================

       stripe_ratio = perturb_width*dx(1)/L(1)

       ! prob_lo(1) < one_fraction_domain2 < center < one_fraction_domain1 < prob_hi(1) 
       one_fraction_domain1=0.5d0*(1.0d0-stripe_ratio)*prob_lo(1)+0.5d0*(1.0d0+stripe_ratio)*prob_hi(1)
       one_fraction_domain2=0.5d0*(1.0d0+stripe_ratio)*prob_lo(1)+0.5d0*(1.0d0-stripe_ratio)*prob_hi(1)

       do i=lo(1),hi(1)
          x1 =(prob_lo(1) + dx(1)*(dble(i)+0.5d0) - one_fraction_domain1)
          x2 =(prob_lo(1) + dx(1)*(dble(i)+0.5d0) - one_fraction_domain2)
       
          if (smoothing_width > 0.d0) then
             ! smooth interface
             do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                n_init(i,j,k,1:nspecies) = n_init_in(2,1:nspecies) + &
                     0.5d0*(n_init_in(2,1:nspecies)-n_init_in(1,1:nspecies)) * &
                     (tanh(x1/(smoothing_width*dx(1))) - tanh(x2/(smoothing_width*dx(1))))
             end do
             end do
          else
             ! discontinuous interface 
             if (x2 < 0.d0 .or. x1 > 0.d0) then   ! outside the stripe
                do k=lo(3),hi(3)
                do j=lo(2),hi(2)
                   n_init(i,j,k,1:nspecies) = n_init_in(2,1:nspecies)
                end do
                end do
             else                                 ! inside the stripe
                do k=lo(3),hi(3)
                do j=lo(2),hi(2)
                   n_init(i,j,k,1:nspecies) = n_init_in(1,1:nspecies)
                end do
                end do
             end if
          end if
       end do

    case(5)
       !=================================================================
       ! bubble having radius = 0.5*perturb_width*dx(1) 
       ! n_init = n_init_in(1,:) inside, n_init = n_init_in (2,:) outside
       ! can be discontinous or smooth depending on smoothing_width
       !=================================================================

       rad = 0.5d0*perturb_width*dx(1) 

       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3) - 0.5d0*(prob_lo(3)+prob_hi(3))
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2) - 0.5d0*(prob_lo(2)+prob_hi(2))
             do i=lo(1),hi(1)
                x = prob_lo(1) + (dble(i)+0.5d0)*dx(1) - 0.5d0*(prob_lo(1)+prob_hi(1))

                r = sqrt(x**2 + y**2 + z**2)

                if (smoothing_width .eq. 0) then
                   ! discontinuous interface
                   if (r .lt. rad) then
                      n_init(i,j,k,1:nspecies) = n_init_in(1,1:nspecies)
                   else
                      n_init(i,j,k,1:nspecies) = n_init_in(2,1:nspecies)
                   end if
                else
                   ! smooth interface
                   n_init(i,j,k,1:nspecies) = n_init_in(1,1:nspecies) + &
                        (n_init_in(2,1:nspecies) - n_init_in(1,1:nspecies))* &
                        0.5d0*(1.d0 + tanh((r-rad)/(smoothing_width*dx(1))))
                end if
             end do
          end do
       end do

    case default

       call bl_error("init_n_3d: prob_type not supported")

    end select

    if(integer_populations) then ! Ensure that the initial number of molecules are integers
       if(initial_variance<0.0d0) then ! Distribute the particles on the box using a multinomial sampler
       
          ! To do this really correctly one should do this for the whole system, not box per box
          if (parallel_IOProcessor()) write(*,*) "Using multinomial initial distribution PER BOX, not domain"
          
          do species=1, nspecies
            call sample_integers(n_init(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),species), &
                     ncells=size(n_init(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),species)), &
                     dv=dx(1)*dx(2)*dx(3)*cross_section)
          end do
       
       else ! Make the number of molecules in each cell Poisson distributed with desired mean
       
          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   call round_to_integers(n_init(i,j,k,1:nspecies), dv=dx(1)*dx(2)*dx(3)*cross_section)
                end do
             end do
          end do  
           
       end if    
    end if

  end subroutine init_n_3d
  
  subroutine round_to_integers(n, dv)
    real(dp_t), intent(inout) :: n(nspecies)
    real(dp_t), intent(in) :: dv ! cell volume
    
    integer :: comp, nparticles
    real(dp_t) :: random(nspecies)
    
    if(initial_variance>0.0d0) then
       do comp=1, nspecies
          ! Generate the initial fluctuations using a Poisson random number generator
          ! This assumes that the distribution of initial conditions is a product Poisson measure
          if (use_bl_rng) then
             call PoissonRNG(number=nparticles, mean=n(comp)*dv, engine=rng_eng_init%p)
          else
             call PoissonRNG(number=nparticles, mean=n(comp)*dv)
          end if
          n(comp) = nparticles/dv
       end do   
    else ! Minimize fluctuations but ensure the number density is an integer
       ! It is important here to use smart rounding so that the average concentration is preserved in expectation
       ! Using nearest integer n = nint(n*dv)/dv will not work except for lots of molecules
       if (use_bl_rng) then
          call UniformRNGs(random, nspecies, engine=rng_eng_init%p)
       else
          call UniformRNGs(random, nspecies)
       end if      
        ! Round to nearest integer
       n = floor( n*dv + random ) / dv
    end if  
    
  end subroutine round_to_integers
  
  subroutine sample_integers(n, ncells, dv)
    integer, intent(in) :: ncells
    real(dp_t), intent(inout) :: n(ncells)
    real(dp_t), intent(in) :: dv ! cell volume
    
    integer, dimension(ncells) :: n_molecules ! Temporary array of integer number of molecules
    real(dp_t), dimension(ncells) :: p
    
    p=n/sum(n) ! Probability proportional to mean density
    p=p/sum(p) ! Make sure it sums to 1 to roundoff again

    if (use_bl_rng) then
       call MultinomialRNG(samples=n_molecules, n_samples=ncells, N=nint(sum(n)*dv), p=p, engine=rng_eng_init%p)
    else
       call MultinomialRNG(samples=n_molecules, n_samples=ncells, N=nint(sum(n)*dv), p=p)
    end if
    
    n = n_molecules/dv ! Convert back to number density
  
  end subroutine sample_integers
  
  subroutine init_n_model(mla,n_init,dx,the_bc_tower,input_array,comp) ! Read from a file

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: n_init(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower ), intent(in   ) :: the_bc_tower
    real(kind=dp_t), intent(in   ) :: input_array(:,:,:)
    integer        , intent(in   ) :: comp

    ! local variables
    integer                        :: lo(mla%dim), hi(mla%dim)
    integer                        :: i, dm, n, nlevs, ng_n
    real(kind=dp_t), pointer       :: np(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "init_n")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_n = n_init(1)%ng

    ! looping over boxes 
    do n=1,nlevs
       do i=1,nfabs(n_init(n))
          np => dataptr(n_init(n),i)
          lo = lwb(get_box(n_init(n),i))
          hi = upb(get_box(n_init(n),i))
          select case (dm)
          case (2)
             call init_n_model_2d(np(:,:,1,:),ng_n,lo,hi,input_array(:,:,1),comp)
          case (3)
             call init_n_model_3d(np,ng_n,lo,hi,input_array,comp)
          end select
       end do
    end do

    call destroy(bpt)

    if (comp .eq. nspecies) then
       do n=1,nlevs
          call multifab_fill_boundary(n_init(n))
          call multifab_physbc(n_init(n),1,scal_bc_comp,nspecies, &
                               the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
       end do
    end if

  end subroutine init_n_model

  subroutine init_n_model_2d(n_init,ng_n,lo,hi,input_array,comp)

    integer         :: lo(:), hi(:), ng_n, comp
    real(kind=dp_t) :: n_init(lo(1)-ng_n:,lo(2)-ng_n:,:)
    real(kind=dp_t) :: input_array(0:,0:) ! This argument is replicated so it is the whole box, not just our patch!
 
    ! local varables
    integer         :: i,j
    
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)

       if(model_file_init>0) then  
          n_init(i,j,comp) = input_array(i,j)
       else
          n_init(i,j,comp) = input_array(j,i)       
       end if   

    end do
    end do

  end subroutine init_n_model_2d

  subroutine init_n_model_3d(n_init,ng_n,lo,hi,input_array,comp)

    integer         :: lo(:), hi(:), ng_n, comp
    real(kind=dp_t) :: n_init(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:,:)
    real(kind=dp_t) :: input_array(0:,0:,0:) ! This argument is replicated so it is the whole box, not just our patch!
 
    ! local varables
    integer         :: i,j,k
    
    do k=lo(3),hi(3)
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)

       if(model_file_init>0) then  
          n_init(i,j,k,comp) = input_array(i,j,k)
       else
          n_init(i,j,k,comp) = input_array(k,j,i)       
       end if   

    end do
    end do
    end do

  end subroutine init_n_model_3d

end module init_n_module
