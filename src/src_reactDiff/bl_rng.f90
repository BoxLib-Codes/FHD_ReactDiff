module bl_rng_module

  use bl_types
  use bl_random_module
  use parallel
  use bl_error_module
  use probin_common_module, only: restart
  use probin_reactdiff_module, only: temporal_integrator, diffusion_type, &
                                     use_Poisson_rng, seed_diffusion, seed_reaction, &
                                     seed_init, integer_populations

  implicit none

  private

  public :: rng_init, rng_destroy, &
            rng_eng_diffusion, &
            rng_eng_reaction, &
            rng_eng_init

  ! randon number engines
  type(bl_rng_engine)      , save :: rng_eng_diffusion
  type(bl_rng_engine)      , save :: rng_eng_reaction
  type(bl_rng_engine)      , save :: rng_eng_init

contains

  subroutine rng_init()

    if (seed_diffusion .eq. -1 .and. restart .lt. 0) then
       call bl_error("seed_diffusion = -1 requires restart")
    end if

    if (seed_reaction .eq. -1 .and. restart .lt. 0) then
       call bl_error("seed_reaction = -1 requires restart")
    end if

    !!!!!!!!!!!!!!!!!!!!!!
    ! build engines
    !!!!!!!!!!!!!!!!!!!!!!

    if (seed_diffusion .ne. -1) then
       call bl_rng_build_engine(rng_eng_diffusion, seed_diffusion)
    end if
    
    if (seed_reaction .ne. -1) then
       call bl_rng_build_engine(rng_eng_reaction, seed_reaction)
    end if

    call bl_rng_build_engine(rng_eng_init, seed_init)

  end subroutine rng_init

  subroutine rng_destroy()

    call bl_rng_destroy_engine(rng_eng_diffusion)
    call bl_rng_destroy_engine(rng_eng_reaction)
    call bl_rng_destroy_engine(rng_eng_init)

  end subroutine rng_destroy

end module bl_rng_module
