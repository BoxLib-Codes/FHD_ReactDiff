module chemical_rates_module

  use ml_layout_module
  use BoxLibRNGs
  use bl_rng_module
  use bl_random_module
  use compute_reaction_rates_module
  use probin_reactdiff_module, only: nspecies, nreactions, stoichiometric_factors, &
                                     temporal_integrator, &
                                     use_Poisson_rng, cross_section, use_bl_rng

  implicit none

  private

  public :: chemical_rates, deterministic_chemical_rates, n_rejections
  
  ! keep track of the number of negative values for linearly combined average reaction rates,
  !  which are to be rejected and set as zero. 
  !  (see the case of lin_comb_avg_react_rate=.true. in chemical_rates_cell)
  integer*8 :: n_rejections = 0 

contains

  ! this routine returns chemical (production) rates with units (number density)/(time) 
  !  for a given state n_cc.
  ! note that chemical *production* rate (per species) is the sum over reactions of
  !  chemical *reaction* rates (per reaction) times with the stoichiometric coefficients.
  ! if optional arguments n_interm and lin_comb_coef_in are given,
  !  chemical rates are obtained by sampling reaction rates from the following average reaction rates:
  !  [lin_comb_coef(1)*f(n_cc)+lin_comb_coef(2)*f(n_interm)]^+
  !  where f(nn) is the average reaction rate of reaction r for state nn.
  ! hence, for the second stage of Mattingly's midpoint scheme (with theta=1/2), lin_comb_coef_in=(-1,2).

  subroutine chemical_rates(mla,n_cc,chem_rate,dx,dt,n_interm,lin_comb_coef_in)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: n_cc(:)
    type(multifab) , intent(inout) :: chem_rate(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:), dt
    type(multifab) , intent(in   ), optional :: n_interm(:)
    real(kind=dp_t), intent(in   ), optional :: lin_comb_coef_in(:)

    ! local
    integer         :: nlevs, dm, n, i
    real(kind=dp_t) :: dv

    type(mfiter)    :: mfi
    type(box)       :: tilebox
    integer         :: tlo(mla%dim), thi(mla%dim)
    integer         ::  lo(mla%dim),  hi(mla%dim)
    integer         :: ng_n, ng_r, ng_i

    real(kind=dp_t), pointer  :: np(:,:,:,:)  ! "n" for n_cc
    real(kind=dp_t), pointer  :: rp(:,:,:,:)  ! "r" for chem_rate 
    real(kind=dp_t), pointer  :: ip(:,:,:,:)  ! "i" for n_interm 

    type(bl_prof_timer), save :: bpt
    
    logical         :: lin_comb_avg_react_rate  ! whether to use linearly combined average reaction rates
    real(kind=dp_t) :: lin_comb_coef(1:2)       ! coefficients for the linear combination 


    nlevs = mla%nlevel
    dm = mla%dim
    
    ! if there are no reactions to process, set rates zero and return. 
    if (nreactions<1) then
      do n=1,nlevs
        call setval(chem_rate(n),0.d0,all=.true.)
      end do
      return
    end if   

    ! otherwise, complete the remaining part
    call build(bpt,"chemical_rates")

    lin_comb_avg_react_rate = .false.
    lin_comb_coef(1) = 1.d0  ! if lin_comb_avg_react_rate=.false.,
    lin_comb_coef(2) = 0.d0  ! lin_comb_coef won't be used.

    if (present(lin_comb_coef_in)) then
      if (present(n_interm)) then
        lin_comb_avg_react_rate = .true.
        lin_comb_coef(1:2) = lin_comb_coef_in(1:2)
      else
        call bl_error("chemical_rates: n_interm missing")
      end if
    else
      if (present(n_interm)) then
        call bl_error("chemical_rates: lin_comb_coef_in missing")
      end if
    end if

    dv = product(dx(1,1:dm))*cross_section


    !!!!! omp tiling

    ng_n = n_cc(1)%ng
    ng_r = chem_rate(1)%ng

    if (lin_comb_avg_react_rate) then
       ng_i = n_interm(1)%ng
    else
       ng_i = n_cc(1)%ng  ! won't be used
    end if

    !$omp parallel private(n,i,mfi,tilebox,tlo,thi) &
    !$omp private(np,rp,ip,lo,hi)

    do n=1,nlevs
      call mfiter_build(mfi, n_cc(n), tiling=.true.)

      do while (more_tile(mfi))
        i = get_fab_index(mfi)

        tilebox = get_tilebox(mfi)
        tlo = lwb(tilebox)
        thi = upb(tilebox)

        np => dataptr(n_cc(n),i)
        rp => dataptr(chem_rate(n),i)

        if (lin_comb_avg_react_rate) then
          ip => dataptr(n_interm(n),i)
        else
          ip => dataptr(n_cc(n),i)  ! won't be used
        end if

        lo = lwb(get_box(n_cc(n),i))
        hi = upb(get_box(n_cc(n),i))

        select case (dm)
        case (2)
          call chemical_rates_2d(np(:,:,1,:),ng_n,rp(:,:,1,:),ng_r,lo,hi,tlo,thi,dv,dt,  &
                                 lin_comb_avg_react_rate,ip(:,:,1,:),ng_i,lin_comb_coef)
        case (3)
          call chemical_rates_3d(np(:,:,:,:),ng_n,rp(:,:,:,:),ng_r,lo,hi,tlo,thi,dv,dt,  &
                                 lin_comb_avg_react_rate,ip(:,:,:,:),ng_i,lin_comb_coef)
        end select
      end do
    end do
    !$omp end parallel

    call destroy(bpt)

  end subroutine chemical_rates 
 
  
  ! this routine returns deterministic chemical (production) rates
  !  irrespective of the value of use_Poisson_rng.
  subroutine deterministic_chemical_rates(mla,n_cc,chem_rate,dx,dt)

    type(ml_layout), intent(in   ) :: mla
    type(multifab ), intent(in   ) :: n_cc(:)
    type(multifab ), intent(inout) :: chem_rate(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:), dt

    ! local
    integer :: use_Poisson_rng_sav

    ! before calling chemical_rates, change the value of use_Poisson_rng as -1
    !  so that deterministic rates are computed.
    ! then, revert the value of use_Poisson_rng.
    use_Poisson_rng_sav = use_Poisson_rng
    use_Poisson_rng = -1
    call chemical_rates(mla,n_cc,chem_rate,dx,dt)
    use_Poisson_rng = use_Poisson_rng_sav

  end subroutine deterministic_chemical_rates


  subroutine chemical_rates_2d(n_cc,ng_n,chem_rate,ng_r,glo,ghi,tlo,thi,dv,dt,      &
                               lin_comb_avg_react_rate,n_interm,ng_i,lin_comb_coef)

    integer        , intent(in   ) :: glo(:), ghi(:), tlo(:), thi(:), ng_n, ng_r, ng_i
    real(kind=dp_t), intent(in   ) :: n_cc     (glo(1)-ng_n:,glo(2)-ng_n:,:)
    real(kind=dp_t), intent(inout) :: chem_rate(glo(1)-ng_r:,glo(2)-ng_r:,:)
    real(kind=dp_t), intent(in   ) :: n_interm (glo(1)-ng_i:,glo(2)-ng_i:,:)
    real(kind=dp_t), intent(in   ) :: dv, dt
    logical        , intent(in   ) :: lin_comb_avg_react_rate 
    real(kind=dp_t), intent(in   ) :: lin_comb_coef(:)
    
    ! local
    integer :: i, j

    do j=tlo(2),thi(2)
    do i=tlo(1),thi(1)
      call chemical_rates_cell(n_cc(i,j,1:nspecies),chem_rate(i,j,1:nspecies),dv,dt,           &
                               lin_comb_avg_react_rate,n_interm(i,j,1:nspecies),lin_comb_coef)
    end do
    end do

  end subroutine chemical_rates_2d


  subroutine chemical_rates_3d(n_cc,ng_n,chem_rate,ng_r,glo,ghi,tlo,thi,dv,dt,      &
                               lin_comb_avg_react_rate,n_interm,ng_i,lin_comb_coef)

    integer        , intent(in   ) :: glo(:), ghi(:), tlo(:), thi(:), ng_n, ng_r, ng_i
    real(kind=dp_t), intent(in   ) :: n_cc     (glo(1)-ng_n:,glo(2)-ng_n:,glo(3)-ng_n:,:)
    real(kind=dp_t), intent(inout) :: chem_rate(glo(1)-ng_r:,glo(2)-ng_r:,glo(3)-ng_r:,:)
    real(kind=dp_t), intent(in   ) :: n_interm (glo(1)-ng_i:,glo(2)-ng_i:,glo(3)-ng_i:,:)
    real(kind=dp_t), intent(in   ) :: dv, dt
    logical,         intent(in   ) :: lin_comb_avg_react_rate 
    real(kind=dp_t), intent(in   ) :: lin_comb_coef(:)
    
    ! local
    integer :: i, j, k

    do k=tlo(3),thi(3)
    do j=tlo(2),thi(2)
    do i=tlo(1),thi(1)
      call chemical_rates_cell(n_cc(i,j,k,1:nspecies),chem_rate(i,j,k,1:nspecies),dv,dt,         &
                               lin_comb_avg_react_rate,n_interm(i,j,k,1:nspecies),lin_comb_coef)
    end do
    end do
    end do

  end subroutine chemical_rates_3d


  subroutine chemical_rates_cell(n_cc,chem_rate,dv,dt,lin_comb_avg_react_rate,n_interm,lin_comb_coef)

    real(kind=dp_t), intent(in   ) :: n_cc(:)
    real(kind=dp_t), intent(inout) :: chem_rate(:)
    real(kind=dp_t), intent(in   ) :: dv, dt
    logical        , intent(in   ) :: lin_comb_avg_react_rate 
    real(kind=dp_t), intent(in   ) :: n_interm(:)
    real(kind=dp_t), intent(in   ) :: lin_comb_coef(:)

    ! local
    integer         :: reaction
    real(kind=dp_t) :: avg_react_rate       (1:nreactions)
    real(kind=dp_t) :: avg_react_rate_interm(1:nreactions)

    ! avg_num_reaction and num_reactions are used in sample_num_reactions
    real(kind=dp_t) :: avg_num_reactions    (1:nreactions)
    real(kind=dp_t) :: num_reactions        (1:nreactions)

    real(kind=dp_t) :: n_tmp(1:nspecies)

    if (use_Poisson_rng .eq. 2) then

       call advance_reaction_SSA_cell(n_cc,n_tmp,dv,dt)
       chem_rate(1:nspecies) = (n_tmp(1:nspecies) - n_cc(1:nspecies)) / dt

    else

       ! obtain average reaction rates
       if (lin_comb_avg_react_rate) then
          ! calculate linearly combined average reaction rates
          call compute_reaction_rates(n_cc    (1:nspecies),avg_react_rate       ,dv)
          call compute_reaction_rates(n_interm(1:nspecies),avg_react_rate_interm,dv)
          avg_react_rate = lin_comb_coef(1)*avg_react_rate + lin_comb_coef(2)*avg_react_rate_interm

          ! check whether the resulting rates are negative - if so, we set them zero below
          do reaction=1,nreactions 
             if (avg_react_rate(reaction) .lt. 0.d0) then 
                n_rejections = n_rejections+1
             end if
          end do
       else
          ! calculate the average reaction rates from n_cc
          call compute_reaction_rates(n_cc(1:nspecies),avg_react_rate,dv)
       end if

       ! convert each average reaction rates into the average number of reactions
       ! zero out any negative reactions 
       avg_num_reactions = max(0.d0,avg_react_rate*dv*dt)

       ! calculate chemical rates 
       chem_rate(1:nspecies) = 0.d0

       do reaction=1,nreactions
          ! for each reaction, sample the number of reaction
          ! sample_num_reactions determines num_reactions(reaction) by using avg_num_reactions(reaction)
          call sample_num_reactions(reaction)

          ! convert this into chemical rates  
          chem_rate(1:nspecies) = chem_rate(1:nspecies) + num_reactions(reaction)/dv/dt *                  &
               (stoichiometric_factors(1:nspecies,2,reaction)-stoichiometric_factors(1:nspecies,1,reaction))
       end do

    end if

  contains

    ! compute num_reactions from avg_num_reactions
    ! note that avg_num_reactions(:) and num_reactions(:) are defined in chemical_rates_cell
    subroutine sample_num_reactions(comp) ! auxilliary routine (should be inlined by compiler)

      integer, intent(in) :: comp

      ! local
      integer :: tmp

      ! for a given reaction, compute how many reactions will happen
      !  by sampling a Poisson (tau leaping) or Gaussian (CLE) random number
      if (avg_num_reactions(comp) .gt. 0.d0) then
         select case(use_Poisson_rng)           
         case(1)
            ! need a Poisson random number for tau leaping
            if (use_bl_rng) then
               call PoissonRNG(number=tmp, mean=avg_num_reactions(comp), engine=rng_eng_reaction%p)
            else
               call PoissonRNG(number=tmp, mean=avg_num_reactions(comp))
            end if
            num_reactions(comp) = tmp ! convert to real
         case(0)
            ! need a Gaussian random number for CLE
            if (use_bl_rng) then
               call NormalRNG(num_reactions(comp), rng_eng_reaction%p)
            else
               call NormalRNG(num_reactions(comp))
            end if
            num_reactions(comp) = avg_num_reactions(comp) + sqrt(avg_num_reactions(comp))*num_reactions(comp)
         case(-1)
            ! do deterministic chemistry   
            num_reactions(comp) = avg_num_reactions(comp)
         case default    
            call bl_error("chemical_rates: invalid use_Poisson_rng")
         end select
      else
         num_reactions(comp) = 0
      end if
      
    end subroutine sample_num_reactions

  end subroutine chemical_rates_cell

  subroutine advance_reaction_SSA_cell(n_old,n_new,dv,dt)

    real(kind=dp_t), intent(in   ) :: n_old(:)
    real(kind=dp_t), intent(inout) :: n_new(:)
    real(kind=dp_t), intent(in   ) :: dv, dt

    real(kind=dp_t) :: avg_react_rate(1:nreactions)
    real(kind=dp_t) :: rTotal, rr, rSum, tau, t_local

    integer :: spec,reaction, which_reaction
    integer :: n_steps_SSA
    
    ! copy old state into new
    n_new = n_old

    t_local = 0.d0
    n_steps_SSA = 0

    EventLoop: do

       ! compute reaction rates in units (# reactions) / (unit time) / (unit volume)
       call compute_reaction_rates(n_new(1:nspecies),avg_react_rate,dv)

       ! compute reaction rates in units (# reactions) / (unit time)
       avg_react_rate = max(0.0d0,avg_react_rate*dv)

       ! sum the reaction rates
       rTotal = sum(avg_react_rate(1:nreactions))

       ! if rTotal is zero (=no reaction), then exit the EventLoop
       if (rTotal .eq. 0.d0) exit EventLoop

       ! generate pseudorandom number in interval [0,1).
       if (use_bl_rng) then
          call UniformRNG(rr, rng_eng_reaction%p)
       else
          call UniformRNG(rr)
       end if
       ! tau is how long until the next reaction occurs
       tau = -log(1-rr)/rTotal
       t_local = t_local + tau;

       if (t_local .gt. dt) exit EventLoop

       ! Select the next reaction according to relative rates
       if (use_bl_rng) then
          call UniformRNG(rr, rng_eng_reaction%p)
       else
          call UniformRNG(rr)
       end if
       rr = rr*rTotal
       rSum = 0
       FindReaction: do reaction=1,nreactions
          rSum = rSum + avg_react_rate(reaction)
          which_reaction = reaction
          if( rSum >= rr ) exit FindReaction
       end do FindReaction

       ! update number densities for this reaction
       do spec=1,nspecies
          n_new(spec) = n_new(spec) + &
               (stoichiometric_factors(spec,2,which_reaction)-stoichiometric_factors(spec,1,which_reaction)) / dv
       end do
       
       n_steps_SSA = n_steps_SSA+1

    end do EventLoop

  end subroutine advance_reaction_SSA_cell

end module chemical_rates_module
