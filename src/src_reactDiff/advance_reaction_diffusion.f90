module advance_reaction_diffusion_module

  use ml_layout_module
  use define_bc_module
  use multifab_physbc_module
  use bc_module
  use stochastic_n_fluxdiv_module
  use diffusive_n_fluxdiv_module
  use chemical_rates_module
  use multinomial_diffusion_module
  use implicit_diffusion_module
  use probin_common_module, only: variance_coef_mass
  use probin_reactdiff_module, only: nspecies, D_Fick, temporal_integrator, &
       midpoint_stoch_flux_type, nreactions, use_Poisson_rng

  implicit none

  private

  public :: advance_reaction_diffusion

contains

  ! this solves dn/dt = div ( D grad (n)) + div (sqrt(2*variance*D*n)*W) + f(n) - g
  !  where f(n) are the chemical production rates (deterministic or stochastic)
  !  and g=ext_src (note minus sign!) is a constant (in time) *deterministic* source term.
  ! To model stochastic particle production (sources) include g in the definition of f instead
  !  or add it as a reaction 0->products

  subroutine advance_reaction_diffusion(mla,n_old,n_new,dx,dt,the_bc_tower,ext_src)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: n_old(:)
    type(multifab) , intent(inout) :: n_new(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(multifab) , intent(in   ), optional :: ext_src(:)

    ! local
    type(multifab) :: diff_fluxdiv(mla%nlevel)
    type(multifab) :: stoch_fluxdiv(mla%nlevel)
    type(multifab) :: diff_coef_face(mla%nlevel,mla%dim)
    type(multifab) :: rate1(mla%nlevel)
    type(multifab) :: rate2(mla%nlevel) 
    type(multifab) :: rhs(mla%nlevel)

    integer :: nlevs,dm,n,i,spec

    type(bl_prof_timer),save :: bpt

    real(kind=dp_t), parameter :: mattingly_lin_comb_coef(1:2) = (/-1.d0, 2.d0/)

    !!!!!!!!
    ! init !
    !!!!!!!!

    call build(bpt,"advance_reaction_diffusion")

    nlevs = mla%nlevel
    dm = mla%dim

    ! build
    do n=1,nlevs
      do i=1,dm
        call multifab_build_edge(diff_coef_face(n,i),mla%la(n),nspecies,0,i)
      end do
      call multifab_build(rate1(n),mla%la(n),nspecies,0)
    end do

    ! diffusion coefficients (for now just setting each to a different constant)
    ! if one wants a space-dependent D or state-dependent D,
    ! see multispecies code as example
    do n=1,nlevs
      do i=1,dm
        do spec=1,nspecies
          call multifab_setval_c(diff_coef_face(n,i),D_Fick(spec),spec,1,all=.true.)
        end do
      end do
    end do

    if (temporal_integrator .eq. -3) then  ! multinomial diffusion 

       ! calculate rates
       ! rates could be deterministic or stochastic depending on use_Poisson_rng
       call chemical_rates(mla,n_old,rate1,dx,dt)

       ! advance multinomial diffusion
       call multinomial_diffusion(mla,n_old,n_new,diff_coef_face,dx,dt,the_bc_tower)
      
       do n=1,nlevs
          ! add reaction contribution and external source
          call multifab_saxpy_3(n_new(n),dt,rate1(n))
          if(present(ext_src)) call multifab_saxpy_3(n_new(n),dt,ext_src(n))
          call multifab_fill_boundary(n_new(n))
          call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                               the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
       end do

       do n=1,nlevs
          do i=1,dm
             call multifab_destroy(diff_coef_face(n,i))
          end do
          call multifab_destroy(rate1(n))
       end do
       return

    end if

    do n=1,nlevs
       call multifab_build(diff_fluxdiv(n),mla%la(n),nspecies,0)
       call multifab_build(stoch_fluxdiv(n),mla%la(n),nspecies,0)
    end do

    ! compute diffusive flux divergence
    call diffusive_n_fluxdiv(mla,n_old,diff_coef_face,diff_fluxdiv,dx,the_bc_tower)

    ! compute stochastic flux divergence
    if (variance_coef_mass .gt. 0.d0) then
      ! fill random flux multifabs with new random numbers
      call fill_mass_stochastic(mla,the_bc_tower%bc_tower_array)
      ! compute stochastic flux divergence
      call stochastic_n_fluxdiv(mla,n_old,diff_coef_face,stoch_fluxdiv,dx,dt, &
                                the_bc_tower,increment_in=.false.)
    else
      do n=1,nlevs
        call multifab_setval(stoch_fluxdiv(n),0.d0,all=.true.)
      end do
    end if

    !!!!!!!!!!!!!!!!
    ! time advance !
    !!!!!!!!!!!!!!!!

    if (temporal_integrator .eq. -1) then  ! forward Euler

      ! calculate rates
      ! rates could be deterministic or stochastic depending on use_Poisson_rng
      call chemical_rates(mla,n_old,rate1,dx,dt)

      ! n_k^{n+1} = n_k^n + dt div (D_k grad n_k)^n
      !                   + dt div (sqrt(2 D_k n_k^n dt) Z) ! Gaussian noise
      !                   + 1/dV * P( f(n_k)*dt*dV )        ! Poisson noise
      !                   + dt ext_src
      do n=1,nlevs
        call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
        call multifab_saxpy_3(n_new(n),dt,diff_fluxdiv(n))
        call multifab_saxpy_3(n_new(n),dt,stoch_fluxdiv(n))
        call multifab_saxpy_3(n_new(n),dt,rate1(n))
        if(present(ext_src)) call multifab_saxpy_3(n_new(n),dt,ext_src(n))

        call multifab_fill_boundary(n_new(n))
        call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                             the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
      end do

   else if (temporal_integrator .eq. -2) then  ! explicit midpoint

      ! temporary storage for second rate
      do n=1,nlevs
         call multifab_build(rate2(n),mla%la(n),nspecies,0)
      end do

      if (use_Poisson_rng .eq. 2) then  ! explicit midpoint with SSA

         !!!!!!!!!!!!!!!
         ! predictor   !
         !!!!!!!!!!!!!!!

         ! n_k^{**} = n_k^n + (dt/2)       div (D_k grad n_k)^n
         !                  + (dt/sqrt(2)) div sqrt(2 D_k n_k^n / (dt*dV)) Z_1 ! Gaussian noise
         !                  + (dt/2)       ext_src
         do n=1,nlevs
            call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
            call multifab_saxpy_3(n_new(n),dt/2.d0,diff_fluxdiv(n))
            call multifab_saxpy_3(n_new(n),dt/sqrt(2.d0),stoch_fluxdiv(n))
            if(present(ext_src)) call multifab_saxpy_3(n_new(n),dt/2.d0,ext_src(n))
         end do
         
         ! computing rate1 = R(n^{**},dt/2) / (dt/2)
         call chemical_rates(mla,n_new,rate1,dx,dt/2.d0)

         ! n_k^* = n_k^{**} + R(n^{**},dt/2)
         do n=1,nlevs
            call multifab_saxpy_3(n_new(n),dt/2.d0,rate1(n))
            call multifab_fill_boundary(n_new(n))
            call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                                 the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
         end do

         !!!!!!!!!!!!!!!
         ! corrector   !
         !!!!!!!!!!!!!!!

         ! compute diffusive flux divergence
         call diffusive_n_fluxdiv(mla,n_new,diff_coef_face,diff_fluxdiv,dx,the_bc_tower)
         
         ! computing rate2 = R(n^*,dt/2) / (dt/2)
         call chemical_rates(mla,n_new,rate2,dx,dt/2.d0)

         ! compute stochastic flux divergence and add to the ones from the predictor stage
         if (variance_coef_mass .gt. 0.d0) then

            ! first, fill random flux multifabs with new random numbers
            call fill_mass_stochastic(mla,the_bc_tower%bc_tower_array)
            call generate_stochastic_fluxdiv_corrector()

         end if

         ! n_k^{n+1} = n_k^n + dt div (D_k grad n_k)^*
         !                   + dt div (sqrt(2 D_k n_k^n / (dt*dV)) Z_1 / sqrt(2) ) ! Gaussian noise
         !                   + dt div (sqrt(2 D_k n_k^? / (dt*dV)) Z_2 / sqrt(2) ) ! Gaussian noise
         !                   + R(n^{**},dt/2)
         !                   + R(n^{*},dt/2)
         !                   + dt ext_src
         ! where
         ! n_k^? = n_k^n               (midpoint_stoch_flux_type=1)
         !       = n_k^pred            (midpoint_stoch_flux_type=2)
         !       = 2*n_k^pred - n_k^n  (midpoint_stoch_flux_type=3)
      
         do n=1,nlevs
            call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
            call multifab_saxpy_3(n_new(n),dt,diff_fluxdiv(n))
            call multifab_saxpy_3(n_new(n),dt/sqrt(2.d0),stoch_fluxdiv(n))
            call multifab_saxpy_3(n_new(n),dt/2.d0,rate1(n))
            call multifab_saxpy_3(n_new(n),dt/2.d0,rate2(n))
            if(present(ext_src)) call multifab_saxpy_3(n_new(n),dt,ext_src(n))
         
            call multifab_fill_boundary(n_new(n))
            call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                                 the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
         end do

      else  ! explicit midpoint for det/tau/CLE

         !!!!!!!!!!!!!!!
         ! predictor   !
         !!!!!!!!!!!!!!!

         ! calculate rates from a(n_old)
         call chemical_rates(mla,n_old,rate1,dx,dt/2.d0)

         ! n_k^{n+1/2} = n_k^n + (dt/2)       div (D_k grad n_k)^n
         !                     + (dt/sqrt(2)) div sqrt(2 D_k n_k^n / (dt*dV)) Z_1 ! Gaussian noise
         !                     + 1/dV * P_1( f(n_k)*(dt/2)*dV )                   ! Poisson noise
         !                     + (dt/2)        ext_src
         do n=1,nlevs
            call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
            call multifab_saxpy_3(n_new(n),dt/2.d0,diff_fluxdiv(n))
            call multifab_saxpy_3(n_new(n),dt/sqrt(2.d0),stoch_fluxdiv(n))
            call multifab_saxpy_3(n_new(n),dt/2.d0,rate1(n))
            if(present(ext_src)) call multifab_saxpy_3(n_new(n),dt/2.d0,ext_src(n))
            
            call multifab_fill_boundary(n_new(n))
            call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                                 the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
         end do

         !!!!!!!!!!!!!!!
         ! corrector   !
         !!!!!!!!!!!!!!!

         ! Here we do not write this in the form that Mattingly et al do
         !  where we just continue the second half of the time step from where we left
         ! Rather, we compute terms at the midpoint and then add contributions from both 
         ! halves of the time step to n_old
         ! This works simpler with diffusion but we have to store both rates1 and rates2
         
         ! compute diffusive flux divergence
         call diffusive_n_fluxdiv(mla,n_new,diff_coef_face,diff_fluxdiv,dx,the_bc_tower)

         ! calculate rates from 2*a(n_pred)-a(n_old)
         call chemical_rates(mla,n_old,rate2,dx,dt/2.d0,n_new,mattingly_lin_comb_coef)

         ! compute stochastic flux divergence and add to the ones from the predictor stage
         if (variance_coef_mass .gt. 0.d0) then

            ! first, fill random flux multifabs with new random numbers
            call fill_mass_stochastic(mla,the_bc_tower%bc_tower_array)
            call generate_stochastic_fluxdiv_corrector()

         end if

         ! n_k^{n+1} = n_k^n + dt div (D_k grad n_k)^{n+1/2}
         !                   + dt div (sqrt(2 D_k n_k^n / (dt*dV)) Z_1 / sqrt(2) ) ! Gaussian noise
         !                   + dt div (sqrt(2 D_k n_k^? / (dt*dV)) Z_2 / sqrt(2) ) ! Gaussian noise
         !                   + 1/dV * P_1( f(n_k)*(dt/2)*dV )                        ! Poisson noise
         !                   + 1/dV * P_2( (2*f(n_k^pred)-f(n_k))*(dt/2)*dV )        ! Poisson noise
         !                   + dt ext_src
         ! where
         ! n_k^? = n_k^n               (midpoint_stoch_flux_type=1)
         !       = n_k^pred            (midpoint_stoch_flux_type=2)
         !       = 2*n_k^pred - n_k^n  (midpoint_stoch_flux_type=3)
         
         do n=1,nlevs
            call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
            call multifab_saxpy_3(n_new(n),dt,diff_fluxdiv(n))
            call multifab_saxpy_3(n_new(n),dt/sqrt(2.d0),stoch_fluxdiv(n))
            call multifab_saxpy_3(n_new(n),dt/2.d0,rate1(n))
            call multifab_saxpy_3(n_new(n),dt/2.d0,rate2(n))
            if(present(ext_src)) call multifab_saxpy_3(n_new(n),dt,ext_src(n))
         
            call multifab_fill_boundary(n_new(n))
            call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                                 the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
         end do

      end if  ! end if/else (use_Poisson_rng)
     
      do n=1,nlevs
         call multifab_destroy(rate2(n))
      end do
     
   else if (temporal_integrator .eq. -4) then  ! implicit midpoint

      if (use_Poisson_rng .eq. 2) then  ! implicit midpoint with SSA

         ! backward Euler predictor to half-time
         ! n_k^* = n_k^n + (dt/2)       div (D_k grad n_k)^{n+1/2}
         !               + (dt/sqrt(2)) div sqrt(2 D_k n_k^n / (dt*dV)) Z_1 ! Gaussian noise
         !               + (dt/2)       ext_src
         !
         ! in delta form
         !
         ! (I - div (dt/2) D_k grad) delta n_k =   (dt/2)       div (D_k grad n_k^n)
         !                                       + (dt/sqrt(2)) div (sqrt(2 D_k n_k^n / (dt*dV)) Z_1
         !                                       + (dt/2) ext_src
         
         do n=1,nlevs
            call multifab_build(rhs(n),mla%la(n),nspecies,0)
         end do

         do n=1,nlevs
            call multifab_setval(rhs(n),0.d0)
            call multifab_saxpy_3(rhs(n),dt/2.d0,diff_fluxdiv(n))
            call multifab_saxpy_3(rhs(n),dt/sqrt(2.d0),stoch_fluxdiv(n))
            if(present(ext_src)) call multifab_saxpy_3(rhs(n),dt/2.d0,ext_src(n))
         end do

         call implicit_diffusion(mla,n_old,n_new,rhs,diff_coef_face,dx,dt,the_bc_tower)

         ! corrector

         ! compute R(n^*,dt) / dt
         call chemical_rates(mla,n_new,rate1,dx,dt)

         ! compute stochastic flux divergence and add to the ones from the predictor stage
         if (variance_coef_mass .gt. 0.d0) then
            
            ! first, fill random flux multifabs with new random numbers
            call fill_mass_stochastic(mla,the_bc_tower%bc_tower_array)

            ! compute n on faces to use in the stochastic flux in the corrector
            ! three possibilities
            call generate_stochastic_fluxdiv_corrector()

         end if

         ! Crank-Nicolson
         ! n_k^{n+1} = n_k^n + (dt/2) div (D_k grad n_k)^n
         !                   + (dt/2) div (D_k grad n_k)^{n+1}
         !                   + dt div (sqrt(2 D_k n_k^n / (dt*dV)) Z_1 / sqrt(2) ) ! Gaussian noise
         !                   + dt div (sqrt(2 D_k n_k^? / (dt*dV)) Z_2 / sqrt(2) ) ! Gaussian noise
         !                   + R(n^*,dt)
         !                   + dt ext_src
         !
         ! in delta form
         !
         ! (I - div (dt/2) D_k grad) delta n_k =   dt div (D_k grad n_k^n)
         !                   + dt div (sqrt(2 D_k n_k^n / (dt*dV)) Z_1 / sqrt(2) ) ! Gaussian noise
         !                   + dt div (sqrt(2 D_k n_k^? / (dt*dV)) Z_2 / sqrt(2) ) ! Gaussian noise
         !                   + R(n^*,dt)
         !                   + dt ext_src

         do n=1,nlevs
            call multifab_setval(rhs(n),0.d0)
            call multifab_saxpy_3(rhs(n),dt,diff_fluxdiv(n))
            call multifab_saxpy_3(rhs(n),dt/sqrt(2.d0),stoch_fluxdiv(n))
            call multifab_saxpy_3(rhs(n),dt,rate1(n))
            if(present(ext_src)) call multifab_saxpy_3(rhs(n),dt,ext_src(n))
         end do
      
         call implicit_diffusion(mla,n_old,n_new,rhs,diff_coef_face,dx,dt,the_bc_tower)

         do n=1,nlevs
            call multifab_destroy(rhs(n))
         end do

      else  ! implicit midpoint for det/tau/CLE

         ! backward Euler predictor to half-time
         ! n_k^{n+1/2} = n_k^n + (dt/2)       div (D_k grad n_k)^{n+1/2}
         !                     + (dt/sqrt(2)) div sqrt(2 D_k n_k^n / (dt*dV)) Z_1 ! Gaussian noise
         !                     + 1/dV * P_1( f(n_k)*(dt/2)*dV )                   ! Poisson noise
         !                     + (dt/2)       ext_src
         !
         ! in delta form
         !
         ! (I - div (dt/2) D_k grad) delta n_k =   (dt/2)       div (D_k grad n_k^n)
         !                                       + (dt/sqrt(2)) div (sqrt(2 D_k n_k^n / (dt*dV)) Z_1
         !                                       + 1/dV * P_1( f(n_k)*(dt/2)*dV )
         !                                       + (dt/2) ext_src
         
         do n=1,nlevs
            call multifab_build(rhs(n),mla%la(n),nspecies,0)
            call multifab_build(rate2(n),mla%la(n),nspecies,0)
         end do

         ! calculate rates
         ! rates could be deterministic or stochastic depending on use_Poisson_rng
         call chemical_rates(mla,n_old,rate1,dx,dt/2.d0)

         do n=1,nlevs
            call multifab_setval(rhs(n),0.d0)
            call multifab_saxpy_3(rhs(n),dt/2.d0,diff_fluxdiv(n))
            call multifab_saxpy_3(rhs(n),dt/sqrt(2.d0),stoch_fluxdiv(n))
            call multifab_saxpy_3(rhs(n),dt/2.d0,rate1(n))
            if(present(ext_src)) call multifab_saxpy_3(rhs(n),dt/2.d0,ext_src(n))
         end do

         call implicit_diffusion(mla,n_old,n_new,rhs,diff_coef_face,dx,dt,the_bc_tower)

         ! corrector

         ! calculate rates from 2*a(n_pred)-a(n_old)
         call chemical_rates(mla,n_old,rate2,dx,dt/2.d0,n_new,mattingly_lin_comb_coef)

         ! compute stochastic flux divergence and add to the ones from the predictor stage
         if (variance_coef_mass .gt. 0.d0) then
            
            ! first, fill random flux multifabs with new random numbers
            call fill_mass_stochastic(mla,the_bc_tower%bc_tower_array)

            ! compute n on faces to use in the stochastic flux in the corrector
            ! three possibilities
            call generate_stochastic_fluxdiv_corrector()

         end if

         ! Crank-Nicolson
         ! n_k^{n+1} = n_k^n + (dt/2) div (D_k grad n_k)^n
         !                   + (dt/2) div (D_k grad n_k)^{n+1}
         !                   + dt div (sqrt(2 D_k n_k^n / (dt*dV)) Z_1 / sqrt(2) ) ! Gaussian noise
         !                   + dt div (sqrt(2 D_k n_k^? / (dt*dV)) Z_2 / sqrt(2) ) ! Gaussian noise
         !                   + 1/dV * P_1( f(n_k)*(dt/2)*dV )                        ! Poisson noise
         !                   + 1/dV * P_2( (2*f(n_k^pred)-f(n_k))*(dt/2)*dV )        ! Poisson noise
         !                   + dt ext_src
         !
         ! in delta form
         !
         ! (I - div (dt/2) D_k grad) delta n_k =   dt div (D_k grad n_k^n)
         !                   + dt div (sqrt(2 D_k n_k^n / (dt*dV)) Z_1 / sqrt(2) ) ! Gaussian noise
         !                   + dt div (sqrt(2 D_k n_k^? / (dt*dV)) Z_2 / sqrt(2) ) ! Gaussian noise
         !                   + 1/dV * P_1( f(n_k)*(dt/2)*dV )                        ! Poisson noise
         !                   + 1/dV * P_2( (2*f(n_k^pred)-f(n_k))*(dt/2)*dV )        ! Poisson noise
         !                   + dt ext_src

         do n=1,nlevs
            call multifab_setval(rhs(n),0.d0)
            call multifab_saxpy_3(rhs(n),dt,diff_fluxdiv(n))
            call multifab_saxpy_3(rhs(n),dt/sqrt(2.d0),stoch_fluxdiv(n))
            call multifab_saxpy_3(rhs(n),dt/2.d0,rate1(n))
            call multifab_saxpy_3(rhs(n),dt/2.d0,rate2(n))
            if(present(ext_src)) call multifab_saxpy_3(rhs(n),dt,ext_src(n))
         end do
      
         call implicit_diffusion(mla,n_old,n_new,rhs,diff_coef_face,dx,dt,the_bc_tower)

         do n=1,nlevs
            call multifab_destroy(rhs(n))
            call multifab_destroy(rate2(n))
         end do

      end if ! end if/else (use_Poisson_rng)

   else
      call bl_error("advance_reaction_diffusion: invalid temporal_integrator")
   end if

    !!!!!!!!!!!
    ! destroy !
    !!!!!!!!!!!

    do n=1,nlevs
      call multifab_destroy(diff_fluxdiv(n))
      call multifab_destroy(stoch_fluxdiv(n))
      do i=1,dm
        call multifab_destroy(diff_coef_face(n,i))
      end do
      call multifab_destroy(rate1(n))
    end do

    call destroy(bpt)
    
  contains
  
      subroutine generate_stochastic_fluxdiv_corrector()
         ! compute n on faces to use in the stochastic flux in the corrector
         ! three possibilities: old, mid, or 2*mid-old
         select case (midpoint_stoch_flux_type)
         case (1)
            ! use n_old
            call stochastic_n_fluxdiv(mla,n_old,diff_coef_face,stoch_fluxdiv,dx,dt, &
                                      the_bc_tower,increment_in=.true.)
         case (2)
            ! use n_pred 
            call stochastic_n_fluxdiv(mla,n_new,diff_coef_face,stoch_fluxdiv,dx,dt, &
                                      the_bc_tower,increment_in=.true.)
         case (3)
            ! compute n_new=2*n_pred-n_old
            ! here we use n_new as temporary storage since it will be overwritten shortly
            do n=1,nlevs
               call multifab_mult_mult_s_c(n_new(n),1,2.d0,nspecies,n_new(n)%ng)
               call multifab_sub_sub_c(n_new(n),1,n_old(n),1,nspecies,n_new(n)%ng)
            end do
            ! use n_new=2*n_pred-n_old
            call stochastic_n_fluxdiv(mla,n_new,diff_coef_face,stoch_fluxdiv,dx,dt, &
                                      the_bc_tower,increment_in=.true.)
         case default
            call bl_error("advance_reaction_diffusion: invalid midpoint_stoch_flux_type")
         end select
      end subroutine generate_stochastic_fluxdiv_corrector    

  end subroutine advance_reaction_diffusion

end module advance_reaction_diffusion_module
