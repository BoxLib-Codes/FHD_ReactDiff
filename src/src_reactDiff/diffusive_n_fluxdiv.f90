module diffusive_n_fluxdiv_module

  use define_bc_module
  use bc_module
  use ml_layout_module
  use cc_applyop_module
  use probin_common_module, only: fixed_dt
  use probin_reactdiff_module, only: nspecies, diffusion_stencil_order

  implicit none

  private

  public :: diffusive_n_fluxdiv

contains

  subroutine diffusive_n_fluxdiv(mla,n_cc,diff_coef_face,diff_fluxdiv,dx,the_bc_tower)

    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(in   )  :: n_cc(:)
    type(multifab) , intent(in   )  :: diff_coef_face(:,:)
    type(multifab) , intent(inout)  :: diff_fluxdiv(:)
    real(kind=dp_t), intent(in   )  :: dx(:,:)
    type(bc_tower) , intent(in   )  :: the_bc_tower

    ! local variables
    integer i,dm,n,nlevs,spec

    type(multifab) :: alpha(mla%nlevel)
    type(multifab) :: res(mla%nlevel)
    type(multifab) :: phi(mla%nlevel)
    type(multifab) :: beta(mla%nlevel,mla%dim)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "diffusive_n_fluxdiv")

    nlevs = mla%nlevel  ! number of levels 
    dm    = mla%dim     ! dimensionality

    ! single cell case set diffusive mass fluxdiv to zero and return
    if ((multifab_volume(n_cc(1))/nspecies)<=1) then
       do n=1,nlevs
          call multifab_setval(diff_fluxdiv(n),0.d0,all=.true.)
       end do
       return
    end if

    do n=1,nlevs
       ! cc_applyop compute res = (alpha - div beta grad) phi
       call multifab_build(alpha(n),mla%la(n),1,0)
       call multifab_build(res(n)  ,mla%la(n),1,0)
       call multifab_build(phi(n)  ,mla%la(n),1,1) 
       do i=1,dm
          call multifab_build_edge(beta(n,i),mla%la(n),1,0,i)
       end do
    end do

    do spec=1,nspecies
       ! cc_applyop compute res = (alpha - div beta grad) phi
       ! phi = n_cc
       ! alpha = 0
       ! beta = -D_k
       do n=1,nlevs
          call multifab_copy_c(phi(n),1,n_cc(n),spec,1,1)
          call multifab_setval(alpha(n),0.d0,all=.true.)
          do i=1,dm
             call multifab_copy_c(beta(n,i),1,diff_coef_face(n,i),spec,1,0)
             call multifab_mult_mult_s(beta(n,i),-1.d0)
          end do
       end do

       call cc_applyop(mla,res,phi,alpha,beta,dx,the_bc_tower,scal_bc_comp+spec-1, &
                       stencil_order_in=diffusion_stencil_order)

       ! copy solution into diff_fluxdiv
       do n=1,nlevs
          call multifab_copy_c(diff_fluxdiv(n),spec,res(n),1,1,0)
       end do
    end do

    do n=1,nlevs
       call multifab_destroy(alpha(n))
       call multifab_destroy(res(n))
       call multifab_destroy(phi(n))
       do i=1,dm
          call multifab_destroy(beta(n,i))
       end do
    end do

    call destroy(bpt)

  end subroutine diffusive_n_fluxdiv

end module diffusive_n_fluxdiv_module
