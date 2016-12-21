module probin_reactdiff_module

  use bl_types
  use bl_space
  use probin_common_module, only: dim_in
 
  implicit none

  integer, parameter :: max_species=10
  integer, parameter :: max_reactions=20

  ! Problem description
  !----------------------
  integer, save         :: nspecies = 2             ! number of species
  integer, save         :: nreactions = 1           ! number of reactions
  
  ! Control of algorithm
  !----------------------
  integer, save :: temporal_integrator = 0          ! 0=D + R (first-order splitting)
                                                    ! 1=(1/2)R + D + (1/2)R (Strang option 1)
                                                    ! 2=(1/2)D + R + (1/2)D (Strang option 2)
                                                    ! -1=unsplit forward Euler
                                                    ! -2=unsplit explicit midpoint 
                                                    ! -3=unsplit multinomial diffusion
                                                    ! -4=unsplit implicit midpoint
  integer, save :: diffusion_type = 0               ! only used for splitting schemes (temporal_integrator>=0)
                                                    ! 0=explicit trapezoidal predictor/corrector
                                                    ! 1=Crank-Nicolson semi-implicit
                                                    ! 2=explicit midpoint
                                                    ! 3=multinomial diffusion
                                                    ! 4=forward Euler  
  integer, save :: reaction_type = 0                ! only used for splitting schemes (temporal_integrator>=0)
                                                    ! 0=first-order (deterministic, tau leaping, CLE, or SSA)
                                                    ! 1=second-order (determinisitc, tau leaping, or CLE only)
  integer, save :: use_Poisson_rng = 1              ! how to calculate chemical production rates
                                                    ! 2=SSA
                                                    ! 1=do tau leaping (Poisson increments)
                                                    ! 0=do CLE (Gaussian increments)
                                                    ! -1=do deterministic chemistry
  integer, save :: midpoint_stoch_flux_type = 1     ! only used for midpoint diffusion schemes (split as well as unsplit)
                                                    ! corrector formulation of noise
                                                    ! 1 = K(nold) * W1 + K(nold)         * W2
                                                    ! 2 = K(nold) * W1 + K(npred)        * W2
                                                    ! 3 = K(nold) * W1 + K(2*npred-nold) * W2

  logical, save :: inhomogeneous_bc_fix = .false.   ! use the Einkemmer boundary condition fix (split schemes only)
  integer, save :: avg_type = 1                     ! how to compute n on faces for stochastic weighting
                                                    ! 1=arithmetic (with C0-Heaviside), 2=geometric, 3=harmonic
                                                    ! 10=arithmetic average with discontinuous Heaviside function
                                                    ! 11=arithmetic average with C1-smoothed Heaviside function
                                                    ! 12=arithmetic average with C2-smoothed Heaviside function

  logical, save :: use_bl_rng = .false.             ! if true, use F_BaseLib/bl_random RNGs
                                                    ! if false, use HydroGrid RNGs

  ! Random number seeds for each physical process for use_bl_rng=T
  ! for positive value, the value is assigned as seed value
  ! for 0, a positive value is randomly chosen
  ! if -1 (only for restart), RNGs status is restored from checkpoint data
  integer, save :: seed_diffusion = 1
  integer, save :: seed_reaction = 1
  integer, save :: seed_init = 1

  ! Initial and boundary conditions
  !----------------------
  real(kind=dp_t), save :: n_init_in(2,max_species) = 1.d0     ! initial values to be used in init_n.f90
  real(kind=dp_t), save :: n_bc(3,2,max_species) = 0.d0        ! n_i boundary conditions (dir,lohi,species)

  integer, save            :: model_file_init = 0              ! initialize from model files:
                                                               ! 0=no, 1=usual order (Fortran), -1=transpose order (C)
  character(len=128), save :: model_file(max_species)          ! one model file for each species
  
  logical, save :: integer_populations=.false.                 ! initialize with all number of molecules strictly integer

  ! Diffusion     
  !----------------------                          
  real(kind=dp_t), save :: D_Fick(max_species) = 1.d0          ! Fickian diffusion coeffs
  integer, save         :: diffusion_stencil_order = 1         ! diffusion boundary stencil order
  integer, save         :: mg_verbose = 0                      ! implicit diffusion solve verbosity
  integer, save         :: cg_verbose = 0                      ! implicit diffusion solve bottom solver verbosity
  real(kind=dp_t), save :: implicit_diffusion_rel_eps = 1.d-10 ! relative eps for implicit diffusion solve
  real(kind=dp_t), save :: implicit_diffusion_abs_eps = -1.d0  ! absolute eps for implicit diffusion solve
  
  ! Chemical reactions
  !----------------------
  real(kind=dp_t), save :: cross_section = 1.d0                ! in 2D, thickness of cell
                                                               ! in general, dv = product(dx(1,1:dm))*cross_section

  ! whether to compute chemical rates using classical LMA or integer-based one
  logical, save         :: include_discrete_LMA_correction = .true. 

  ! LMA chemical reaction rate for each reaction (assuming Law of Mass holds)
  real(kind=dp_t), save :: rate_const(max_reactions) = 0.0d0, rate_multiplier=1.0d0

  ! stoichiometric factors for each reaction (species,LHS(1)/RHS(2),reaction)
  ! Example: For N1 + 2*N2 -> N3 use
  ! stoichiometric_factors(1:3,1,1) = 1 2 0
  ! stoichiometric_factors(1:3,2,1) = 0 0 1
  integer, save         :: stoichiometric_factors(max_species,2,max_reactions) = 0 
  
  ! Controlling output:
  integer, save :: n_steps_write_avg = 0 ! If non-zero, its absolute value tells how many steps before writing total densites
                                         ! If positive, it writes average number densities in the system
                                         ! If negative, it writes the total number of molecules in the system
  ! interval to write histograms
  integer, save :: hist_int = 0


  namelist /probin_reactdiff/ nspecies, nreactions
  namelist /probin_reactdiff/ temporal_integrator, diffusion_type, midpoint_stoch_flux_type
  namelist /probin_reactdiff/ reaction_type, use_Poisson_rng, avg_type
  namelist /probin_reactdiff/ use_bl_rng, seed_diffusion, seed_reaction, seed_init
  namelist /probin_reactdiff/ inhomogeneous_bc_fix, n_init_in, n_bc, model_file_init, model_file, integer_populations
  namelist /probin_reactdiff/ D_Fick, diffusion_stencil_order, mg_verbose, cg_verbose
  namelist /probin_reactdiff/ implicit_diffusion_rel_eps, implicit_diffusion_abs_eps
  namelist /probin_reactdiff/ cross_section, include_discrete_LMA_correction
  namelist /probin_reactdiff/ rate_const, rate_multiplier, stoichiometric_factors
  namelist /probin_reactdiff/ n_steps_write_avg, hist_int

contains

  subroutine probin_reactdiff_init()

    use f2kcli
    use parallel
    use bl_IO_module
    use bl_prof_module
    use bl_error_module
    use bl_constants_module
    use cluster_module
    
    integer            :: narg, farg
    character(len=128) :: fname
    integer            :: un
    logical            :: lexist,need_inputs
    
    narg = command_argument_count()

    ! You can put default values here if you want, but we have specified them above 
    ! in the variable declaration
 
    ! read from input file 
    need_inputs = .true.
    farg = 1
    if ( need_inputs .AND. narg >= 1 ) then
       call get_command_argument(farg, value = fname)
       inquire(file = fname, exist = lexist )
       if ( lexist ) then
          farg = farg + 1
          un = unit_new()
          open(unit=un, file = fname, status = 'old', action = 'read')
          read(unit=un, nml = probin_reactdiff)
          close(unit=un)
          need_inputs = .false.
       end if
    end if
    
    ! also can be read in from the command line by appending 
    do while ( farg <= narg )
       call get_command_argument(farg, value = fname)
       select case (fname)

       case ('--nspecies')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) nspecies

       case ('--nreactions')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) nreactions

       case ('--temporal_integrator')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) temporal_integrator 

       case ('--diffusion_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) diffusion_type

       case ('--midpoint_stoch_flux_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) midpoint_stoch_flux_type

       case ('--reaction_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) reaction_type

       case ('--use_Poisson_rng')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) use_Poisson_rng

       case ('--inhomogeneous_bc_fix')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) inhomogeneous_bc_fix

       case ('--integer_populations')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) integer_populations 

       case ('--avg_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) avg_type

       case ('--use_bl_rng')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) use_bl_rng

       case ('--seed_diffusion')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) seed_diffusion

       case ('--seed_reaction')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) seed_reaction

       case ('--seed_init')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) seed_init

       case ('--model_file_init')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) model_file_init

       case ('--diffusion_stencil_order')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) diffusion_stencil_order

       case ('--mg_verbose')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mg_verbose

       case ('--cg_verbose')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) cg_verbose

       case ('--implicit_diffusion_rel_eps')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) implicit_diffusion_rel_eps

       case ('--implicit_diffusion_abs_eps')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) implicit_diffusion_abs_eps

       case ('--cross_section')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) cross_section

       case ('--rate_multiplier')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) rate_multiplier

       case ('--include_discrete_LMA_correction')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) include_discrete_LMA_correction

       case ('--n_steps_write_avg')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) n_steps_write_avg 

       case ('--hist_int')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) hist_int

       case ('--D_Fick_1')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) D_Fick(1)

       case ('--D_Fick_2')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) D_Fick(2)

       case ('--D_Fick_3')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) D_Fick(3)

       case ('--no_diffusion')
          D_Fick(1:max_species) = 0.d0

       case default
          if (parallel_IOProcessor() ) then
             print*,'probin_reactdiff: command-line input ',trim(fname),' not read'
          end if

       end select

       farg = farg + 1
    end do
    
    ! check that nspecies<=max_species, otherwise abort with error message
    if(nspecies.gt.max_species) then 
       call bl_error(" nspecies greater than max_species - Aborting")
    end if
    
    ! check that nreactions<=max_reactions, otherwise abort with error message
    if(nreactions.gt.max_reactions) then 
       call bl_error(" nreactions greater than max_reactions - Aborting")
    end if

    if (inhomogeneous_bc_fix .and. temporal_integrator .lt. 0) then
       call bl_error("inhomogeneous_bc_fix only appropriate for split schemes")
    end if

    if (temporal_integrator .ge. 0 .and. reaction_type .ne. 0) then
       if (use_Poisson_rng .eq. 2) then
          call bl_error("SSA (use_Poisson_rng=2) requires reaction_type=0 for split schemes")
       end if
    end if
    
  end subroutine probin_reactdiff_init

end module probin_reactdiff_module
