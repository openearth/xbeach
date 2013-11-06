module libxbeach_dynamic

  use libxbeach_module
  use iso_c_binding
  use logging_module

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
  
  subroutine assignlogdelegate(fPtr) bind(C,name="assignlogdelegate")
    !DEC$ ATTRIBUTES DLLEXPORT::assignlogdelegate
    use logging_module

    procedure(distributeloginterface) :: fPtr
    
    call assignlogdelegate_internal(fPtr)
    
  end subroutine assignlogdelegate
  
  integer(c_int) function writetolog()  bind(C, name="writetolog")
    !DEC$ ATTRIBUTES DLLEXPORT::writetolog
    use logging_module
    
    integer :: i, tmp

    writetolog = -1
    
    do i = 1, 10
        tmp = i
        call distributelog(tmp,"test iets anders",16)
    end do
    
    writetolog = 0
  end function writetolog
  
    end module libxbeach_dynamic