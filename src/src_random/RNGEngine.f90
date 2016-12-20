module RNGEngine
   ! The RNG egines are written in C and C++ and part of BoxLib
   ! This is just a collection of interfaces to the C codes callable from Fortran
   use iso_c_binding
   implicit none
   public
   
   integer, parameter, private :: dp=kind(0.0d0), sp=kind(0.0)

   ! These are written in C and part of BoxLib
   interface SeedRNG_C
      subroutine srandgen(seed) bind(c)
         ! Seed the generator (should be called at the beginning)
         ! Use the shell script $RANDOM to generate a truly random seed
         integer, value :: seed
      end subroutine
   end interface

! This is not used at present and since it requires a C++ compiler I leave it out here   
!   interface SeedParallelRNG_C
!      subroutine SeedParallelRNG(seed) bind(c, name="SeedParallelRNG")
!         ! Seed the generator on each processor (should be called at the beginning)
!         ! This requires MPI to be initialized and running and is compiled with C++
!         integer, value :: seed
!      end subroutine
!   end interface   

   ! Donev: We have observed some difficulties calling this from Fortran
   ! It either has to do with unsigned versus signed integers or passing by value
   ! Appears to be a compiler bug in gfortran as it works in ifort.
   ! For now, best avoid this routine and use the floating-point RNG to compute integers
   interface UniformInteger
      ! void genrandint (unsigned long int *r, unsigned long int n)
      subroutine genrandint(number, range) bind(c)
         ! Returns an integer uniformly distributed in the range [1,range]
         import
         integer, intent(out) :: number
         integer, intent(in), value :: range
         !integer, intent(in) :: range ! Donev: Temporary change to pass by address here
      end subroutine   
   end interface

   ! Returns pseudorandom number in interval [0,1).
   interface UniformRNG_C
      subroutine genrand(number) bind(c) ! Double precision
         import
         real(dp), intent(out) :: number
      end subroutine   
      subroutine genrand_sp(number) bind(c) ! Single-precision
         import
         real(sp), intent(out) :: number
      end subroutine   
      subroutine hg_genrand(number, engine) bind(c)
        import
        real(dp), intent(out) :: number
        type(c_ptr), value :: engine
      end subroutine hg_genrand
      subroutine hg_genrand_sp(number, engine) bind(c)
        import
        real(sp), intent(out) :: number
        type(c_ptr), value :: engine
      end subroutine hg_genrand_sp
   end interface
   
   ! Donev: Consider using the Fortran version in RNG.f90 instead of NormalRNG_C from the code   
   interface NormalRNG_C
      subroutine genrandn(number) bind(c)
         ! Returns a normally-distributed number with mean 0 and variance 1
         import
         real(dp), intent(out) :: number
      end subroutine 
      subroutine hg_genrandn(number, engine) bind(c)
        ! Returns a normally-distributed number with mean 0 and variance 1
        import
        real(dp), intent(out) :: number
        type(c_ptr), value :: engine
      end subroutine hg_genrandn
   end interface

end module
