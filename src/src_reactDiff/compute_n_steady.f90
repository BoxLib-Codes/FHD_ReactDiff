module compute_n_steady_module

  use ml_layout_module
  use multifab_physbc_module
  use define_bc_module
  use bc_module
  use bndry_reg_module
  use ml_solve_module
  use probin_reactdiff_module, only: nspecies, mg_verbose, cg_verbose, &
       implicit_diffusion_rel_eps, implicit_diffusion_abs_eps, D_Fick

  implicit none

  private

  public :: compute_n_steady

contains

  subroutine compute_n_steady(mla,n_steady,dx,dt,the_bc_tower)

    ! compute div D_k grad n_steady_k = 0 with inhomogeneous boundary conditions
    
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: n_steady(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! for multigrid solver; (alpha - div beta grad) phi = rhs
    type(multifab) :: alpha(mla%nlevel)
    type(multifab) :: rhs(mla%nlevel)
    type(multifab) :: phi(mla%nlevel)
    type(multifab) :: beta(mla%nlevel,mla%dim)

    ! for diffusion multigrid - not used but needs to be passed in
    type(bndry_reg) :: fine_flx(mla%nlevel)

    integer :: i,dm,n,nlevs,spec

    type(bl_prof_timer),save :: bpt

    nlevs = mla%nlevel
    dm = mla%dim

    ! if we are periodic in all directions, set n_steady=0 and return
    if (all(mla%pmask(1:dm))) then
       do n=1,nlevs
          call multifab_setval(n_steady(n),0.d0,all=.true.)
       end do
       return
    end if

    call build(bpt,"compute_n_steady")

    ! for multigrid solver; (alpha - div beta grad) phi = rhs
    do n=1,nlevs
       call multifab_build(alpha(n),mla%la(n),1,0)
       call multifab_build(rhs(n)  ,mla%la(n),1,0)
       call multifab_build(phi(n)  ,mla%la(n),1,1)
       do i=1,dm
          call multifab_build_edge(beta(n,i),mla%la(n),1,0,i)
       end do
    end do

    ! stores beta*grad phi/dx_fine on coarse-fine interfaces
    ! this gets computed inside of ml_cc_solve
    ! we pass it back out because some algorithms (like projection methods) 
    ! use this information
    do n = 1,nlevs
       call bndry_reg_build(fine_flx(n),mla%la(n),ml_layout_get_pd(mla,n))
    end do

    do n=1,nlevs
       call multifab_setval(alpha(n),0.d0,all=.true.)
       call multifab_setval(rhs(n),0.d0,all=.true.)
    end do

    do spec=1,nspecies


       ! beta = D_k
       do n=1,nlevs
          do i=1,dm
             call multifab_setval(beta(n,i),D_Fick(spec),all=.true.)
          end do
       end do

       ! initial guess for phi
       ! zero on the interior, ghost cells filled with Dirichlet value
       do n=1,nlevs
          call multifab_setval(phi(n),0.d0,all=.true.)
          call multifab_physbc(phi(n),1,scal_bc_comp+spec-1,1, &
                               the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
       end do

       ! solve the Poisson equation
       call ml_cc_solve(mla,rhs,phi,fine_flx,alpha,beta,dx, &
                        the_bc_tower,scal_bc_comp+spec-1, &
                        stencil_order=1, &
                        verbose=mg_verbose, &
                        cg_verbose=cg_verbose, &
                        eps=implicit_diffusion_rel_eps, &
                        abs_eps=implicit_diffusion_abs_eps)

       ! copy solution into n_new
       do n=1,nlevs
          call multifab_copy_c(n_steady(n),spec,phi(n),1,1,0)
       end do

    end do

    do n=1,nlevs
       call multifab_destroy(alpha(n))
       call multifab_destroy(rhs(n))
       call multifab_destroy(phi(n))
       do i=1,dm
          call multifab_destroy(beta(n,i))
       end do
    end do

    do n = 1,nlevs
       call bndry_reg_destroy(fine_flx(n))
    end do

    call destroy(bpt)

  end subroutine compute_n_steady

end module compute_n_steady_module
