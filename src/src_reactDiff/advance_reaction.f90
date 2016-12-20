module advance_reaction_module

  use ml_layout_module
  use define_bc_module
  use multifab_physbc_module
  use bc_module
  use chemical_rates_module
  use probin_reactdiff_module, only: nspecies, nreactions, reaction_type, inhomogeneous_bc_fix

  implicit none

  private

  public :: advance_reaction

  ! here we use Mattingly's predictor-corrector with theta=0.5d0 (for rection_type=1).
  ! with these parameters this is actually equivalent to a traditional midpoint scheme.
  real(kind=dp_t), parameter :: theta = 0.5d0
  real(kind=dp_t), parameter :: alpha1 = 2.d0
  real(kind=dp_t), parameter :: alpha2 = 1.d0

contains

  ! this solves dn/dt = f(n) - g (note the minus sign for g)
  ! where f(n) are the chemical production rates (deterministic or stochastic)
  ! and g=ext_src is a constant (in time) *deterministic* source term.
  ! to model stochastic particle production (sources) include g in the definition of f instead.
  ! or add it as a reaction 0->products

  subroutine advance_reaction(mla,n_old,n_new,dx,dt,the_bc_tower,ext_src)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: n_old(:)
    type(multifab) , intent(inout) :: n_new(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(multifab) , intent(in   ) :: ext_src(:)

    ! local
    type(multifab) :: rate(mla%nlevel)

    integer :: nlevs, dm, n

    type(bl_prof_timer),save :: bpt

    real(kind=dp_t) :: mattingly_lin_comb_coef(1:2)  ! only used for reaction_type=1

    !!!!!!!!
    ! init !
    !!!!!!!!

    nlevs = mla%nlevel
    dm = mla%dim

    ! if there are no reactions to process, copy n_old to n_new,
    ! account for ext_src and return
    if(nreactions<1) then
       do n=1,nlevs
          call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,n_new(n)%ng)
          call multifab_saxpy_3(n_new(n),-dt,ext_src(n))
       end do
       return
    end if
    
    call build(bpt,"advance_reaction")

    do n=1,nlevs
       call multifab_build(rate(n),mla%la(n),nspecies,0)
    end do

    !!!!!!!!!!!!!!!!!!
    ! advancing time !
    !!!!!!!!!!!!!!!!!!

    if (reaction_type .eq. 0) then ! first-order det/tau-leaping/CLE, or SSA

      ! calculate rates
      ! rates could be deterministic or stochastic depending on use_Poisson_rng
      call chemical_rates(mla,n_old,rate,dx,dt)

      ! update
      do n=1,nlevs
        call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
        call multifab_saxpy_3(n_new(n),dt,rate(n))
        call multifab_saxpy_3(n_new(n),-dt,ext_src(n))  ! note the negative sign

        call multifab_fill_boundary(n_new(n))
        call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                             the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
      end do

    else if (reaction_type .eq. 1) then  ! second-order det/tau-leaping/CLE

      !!!!!!!!!!!!!!!
      ! predictor   !
      !!!!!!!!!!!!!!!

      ! calculate rates from a(n_old)
      call chemical_rates(mla,n_old,rate,dx,theta*dt)

      do n=1,nlevs
        call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
        call multifab_saxpy_3(n_new(n),theta*dt,rate(n))
        call multifab_saxpy_3(n_new(n),-theta*dt,ext_src(n))  ! note the negative sign

        call multifab_fill_boundary(n_new(n))
        call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                             the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
      end do

      !!!!!!!!!!!!!!!
      ! corrector   !
      !!!!!!!!!!!!!!!
      
      ! Here we write this in the form that Mattingly et al do
      !  where we just continue the second half of the time step from where we left

      mattingly_lin_comb_coef(1) = -alpha2
      mattingly_lin_comb_coef(2) = alpha1 

      ! calculate rates from 2*a(n_pred)-a(n_old)
      call chemical_rates(mla,n_old,rate,dx,(1.d0-theta)*dt,n_new,mattingly_lin_comb_coef)

      ! update
      do n=1,nlevs
        call multifab_saxpy_3(n_new(n),(1.d0-theta)*dt,rate(n))
        ! note the negative sign
        ! also note that ext_src does not change in the time interval (t,t+dt) 
        call multifab_saxpy_3(n_new(n),-dt*(1.d0-theta)*(alpha1-alpha2),ext_src(n))

        call multifab_fill_boundary(n_new(n))
        call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                             the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
      end do

    else

      call bl_error("advance_reaction: invalid reaction_type")

    end if

    !!!!!!!!!!!
    ! destroy !
    !!!!!!!!!!!
    do n=1,nlevs
        call multifab_destroy(rate(n))
    end do

    call destroy(bpt)

  end subroutine advance_reaction

end module advance_reaction_module
