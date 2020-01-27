module dns_libis
  use,intrinsic :: iso_c_binding
  implicit none

  interface
    subroutine libis_send(udata, nx32, ny, nz32) bind(C, name="libis_send")
      use,intrinsic :: iso_c_binding
      implicit none
      real(c_double) :: udata(nx32, ny, nz32)
      integer(c_intptr_t) :: nx32, ny, nz32
    end subroutine libis_send
  end interface

  interface
    subroutine libis_initialize(udata, nx32, ny, nz32, Xyst, Xysz, Xzst, Xzsz, LMFx, LMFy, LMFz) bind(C, name="libis_initialize")
      use,intrinsic :: iso_c_binding
      implicit none
      integer(c_intptr_t) :: nx32, ny, nz32
      integer(c_intptr_t) :: Xyst, Xysz, Xzst, Xzsz
      real(c_double) :: LMFx, LMFy, LMFz
      real(c_double) :: udata(nx32, Xysz, Xzsz)
    end subroutine libis_initialize
  end interface

  interface
    subroutine libis_finalize() bind(C, name="libis_finalize")
      use,intrinsic :: iso_c_binding
      implicit none
    end subroutine libis_finalize
  end interface

end module dns_libis

