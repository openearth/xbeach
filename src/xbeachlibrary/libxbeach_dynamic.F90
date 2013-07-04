module libxbeach_dynamic

  use libxbeach_module
  use iso_c_binding

  implicit none

contains

  integer(c_int) function xbeach_init() bind(C, name="init")
    !DEC$ ATTRIBUTES DLLEXPORT::xbeach_init

    xbeach_init = init()

  end function xbeach_init

  integer(c_int) function xbeach_outputext() bind(C, name="outputext")
    !DEC$ ATTRIBUTES DLLEXPORT::xbeach_outputext

    xbeach_outputext = outputext()

  end function xbeach_outputext

  integer(c_int) function xbeach_executestep() bind(C, name="executestep")
    !DEC$ ATTRIBUTES DLLEXPORT::xbeach_executestep

    xbeach_executestep = executestep()

  end function xbeach_executestep

  integer(c_int) function xbeach_finalize() bind(C, name="finalize")
    !DEC$ ATTRIBUTES DLLEXPORT::xbeach_finalize

    xbeach_finalize = finalize()

  end function xbeach_finalize

end module libxbeach_dynamic
