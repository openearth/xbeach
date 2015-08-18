

  subroutine set_var(c_var_name, xptr) bind(C, name="set_var")
    !DEC$ ATTRIBUTES DLLEXPORT :: set_var
    ! Return a pointer to the variable
    use iso_c_binding, only: c_double, c_char, c_loc, c_f_pointer

    character(kind=c_char), intent(in) :: c_var_name(*)
    type(c_ptr), value, intent(in) :: xptr

    real(c_double), pointer :: x_0d_double_ptr
    real(c_double), pointer :: x_1d_double_ptr(:)
    real(c_double), pointer :: x_2d_double_ptr(:,:)
    real(c_double), pointer :: x_3d_double_ptr(:,:,:)
    real(c_double), pointer :: x_4d_double_ptr(:,:,:,:)
    integer(c_int), pointer :: x_0d_int_ptr
    integer(c_int), pointer :: x_1d_int_ptr(:)
    integer(c_int), pointer :: x_2d_int_ptr(:,:)
    integer(c_int), pointer :: x_3d_int_ptr(:,:,:)
    real(c_float), pointer  :: x_1d_float_ptr(:)
    real(c_float), pointer  :: x_2d_float_ptr(:,:)
    real(c_float), pointer  :: x_3d_float_ptr(:,:,:)
    ! The fortran name of the attribute name
    character(len=strlen(c_var_name)) :: var_name

    ! Store the name
    ! The real array
    integer :: index
    type(arraytype) :: array

    var_name = char_array_to_string(c_var_name)

    index =  chartoindex(var_name)
    if (index .eq. -1) return
    call indextos(s,index,array)

    select case(var_name)
    case("x")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%x))
!       write(*,*) "setting x to", x_2d_double_ptr
       s%x(:,:) = x_2d_double_ptr
    case("y")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%y))
!       write(*,*) "setting y to", x_2d_double_ptr
       s%y(:,:) = x_2d_double_ptr
    case("xz")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%xz))
!       write(*,*) "setting xz to", x_2d_double_ptr
       s%xz(:,:) = x_2d_double_ptr
    case("yz")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%yz))
!       write(*,*) "setting yz to", x_2d_double_ptr
       s%yz(:,:) = x_2d_double_ptr
    case("xu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%xu))
!       write(*,*) "setting xu to", x_2d_double_ptr
       s%xu(:,:) = x_2d_double_ptr
    case("yu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%yu))
!       write(*,*) "setting yu to", x_2d_double_ptr
       s%yu(:,:) = x_2d_double_ptr
    case("xv")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%xv))
!       write(*,*) "setting xv to", x_2d_double_ptr
       s%xv(:,:) = x_2d_double_ptr
    case("yv")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%yv))
!       write(*,*) "setting yv to", x_2d_double_ptr
       s%yv(:,:) = x_2d_double_ptr
    case("dsu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dsu))
!       write(*,*) "setting dsu to", x_2d_double_ptr
       s%dsu(:,:) = x_2d_double_ptr
    case("dsv")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dsv))
!       write(*,*) "setting dsv to", x_2d_double_ptr
       s%dsv(:,:) = x_2d_double_ptr
    case("dsz")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dsz))
!       write(*,*) "setting dsz to", x_2d_double_ptr
       s%dsz(:,:) = x_2d_double_ptr
    case("dsc")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dsc))
!       write(*,*) "setting dsc to", x_2d_double_ptr
       s%dsc(:,:) = x_2d_double_ptr
    case("dnu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dnu))
!       write(*,*) "setting dnu to", x_2d_double_ptr
       s%dnu(:,:) = x_2d_double_ptr
    case("dnv")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dnv))
!       write(*,*) "setting dnv to", x_2d_double_ptr
       s%dnv(:,:) = x_2d_double_ptr
    case("dnz")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dnz))
!       write(*,*) "setting dnz to", x_2d_double_ptr
       s%dnz(:,:) = x_2d_double_ptr
    case("dnc")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dnc))
!       write(*,*) "setting dnc to", x_2d_double_ptr
       s%dnc(:,:) = x_2d_double_ptr
    case("dsdnui")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dsdnui))
!       write(*,*) "setting dsdnui to", x_2d_double_ptr
       s%dsdnui(:,:) = x_2d_double_ptr
    case("dsdnvi")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dsdnvi))
!       write(*,*) "setting dsdnvi to", x_2d_double_ptr
       s%dsdnvi(:,:) = x_2d_double_ptr
    case("dsdnzi")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dsdnzi))
!       write(*,*) "setting dsdnzi to", x_2d_double_ptr
       s%dsdnzi(:,:) = x_2d_double_ptr
    case("alfaz")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%alfaz))
!       write(*,*) "setting alfaz to", x_2d_double_ptr
       s%alfaz(:,:) = x_2d_double_ptr
    case("alfau")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%alfau))
!       write(*,*) "setting alfau to", x_2d_double_ptr
       s%alfau(:,:) = x_2d_double_ptr
    case("alfav")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%alfav))
!       write(*,*) "setting alfav to", x_2d_double_ptr
       s%alfav(:,:) = x_2d_double_ptr
    case("sdist")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%sdist))
!       write(*,*) "setting sdist to", x_2d_double_ptr
       s%sdist(:,:) = x_2d_double_ptr
    case("ndist")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%ndist))
!       write(*,*) "setting ndist to", x_2d_double_ptr
       s%ndist(:,:) = x_2d_double_ptr
    case("dx")
       call c_f_pointer(xptr, x_0d_double_ptr, shape(s%dx))
!       write(*,*) "setting dx to", x_0d_double_ptr
       s%dx = x_0d_double_ptr
    case("dy")
       call c_f_pointer(xptr, x_0d_double_ptr, shape(s%dy))
!       write(*,*) "setting dy to", x_0d_double_ptr
       s%dy = x_0d_double_ptr
    case("xori")
       call c_f_pointer(xptr, x_0d_double_ptr, shape(s%xori))
!       write(*,*) "setting xori to", x_0d_double_ptr
       s%xori = x_0d_double_ptr
    case("yori")
       call c_f_pointer(xptr, x_0d_double_ptr, shape(s%yori))
!       write(*,*) "setting yori to", x_0d_double_ptr
       s%yori = x_0d_double_ptr
    case("alfa")
       call c_f_pointer(xptr, x_0d_double_ptr, shape(s%alfa))
!       write(*,*) "setting alfa to", x_0d_double_ptr
       s%alfa = x_0d_double_ptr
    case("posdwn")
       call c_f_pointer(xptr, x_0d_double_ptr, shape(s%posdwn))
!       write(*,*) "setting posdwn to", x_0d_double_ptr
       s%posdwn = x_0d_double_ptr
    case("thetamin")
       call c_f_pointer(xptr, x_0d_double_ptr, shape(s%thetamin))
!       write(*,*) "setting thetamin to", x_0d_double_ptr
       s%thetamin = x_0d_double_ptr
    case("thetamax")
       call c_f_pointer(xptr, x_0d_double_ptr, shape(s%thetamax))
!       write(*,*) "setting thetamax to", x_0d_double_ptr
       s%thetamax = x_0d_double_ptr
    case("nx")
       call c_f_pointer(xptr, x_0d_int_ptr, shape(s%nx))
!       write(*,*) "setting nx to", x_0d_int_ptr
       s%nx = x_0d_int_ptr
    case("ny")
       call c_f_pointer(xptr, x_0d_int_ptr, shape(s%ny))
!       write(*,*) "setting ny to", x_0d_int_ptr
       s%ny = x_0d_int_ptr
    case("zs01")
       call c_f_pointer(xptr, x_0d_double_ptr, shape(s%zs01))
!       write(*,*) "setting zs01 to", x_0d_double_ptr
       s%zs01 = x_0d_double_ptr
    case("zs02")
       call c_f_pointer(xptr, x_0d_double_ptr, shape(s%zs02))
!       write(*,*) "setting zs02 to", x_0d_double_ptr
       s%zs02 = x_0d_double_ptr
    case("zs03")
       call c_f_pointer(xptr, x_0d_double_ptr, shape(s%zs03))
!       write(*,*) "setting zs03 to", x_0d_double_ptr
       s%zs03 = x_0d_double_ptr
    case("zs04")
       call c_f_pointer(xptr, x_0d_double_ptr, shape(s%zs04))
!       write(*,*) "setting zs04 to", x_0d_double_ptr
       s%zs04 = x_0d_double_ptr
    case("xyzs01")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%xyzs01))
!       write(*,*) "setting xyzs01 to", x_1d_double_ptr
       s%xyzs01(:) = x_1d_double_ptr
    case("xyzs02")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%xyzs02))
!       write(*,*) "setting xyzs02 to", x_1d_double_ptr
       s%xyzs02(:) = x_1d_double_ptr
    case("xyzs03")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%xyzs03))
!       write(*,*) "setting xyzs03 to", x_1d_double_ptr
       s%xyzs03(:) = x_1d_double_ptr
    case("xyzs04")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%xyzs04))
!       write(*,*) "setting xyzs04 to", x_1d_double_ptr
       s%xyzs04(:) = x_1d_double_ptr
    case("tidelen")
       call c_f_pointer(xptr, x_0d_int_ptr, shape(s%tidelen))
!       write(*,*) "setting tidelen to", x_0d_int_ptr
       s%tidelen = x_0d_int_ptr
    case("windlen")
       call c_f_pointer(xptr, x_0d_int_ptr, shape(s%windlen))
!       write(*,*) "setting windlen to", x_0d_int_ptr
       s%windlen = x_0d_int_ptr
    case("zb")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%zb))
!       write(*,*) "setting zb to", x_2d_double_ptr
       s%zb(:,:) = x_2d_double_ptr
    case("zb0")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%zb0))
!       write(*,*) "setting zb0 to", x_2d_double_ptr
       s%zb0(:,:) = x_2d_double_ptr
    case("theta")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%theta))
!       write(*,*) "setting theta to", x_1d_double_ptr
       s%theta(:) = x_1d_double_ptr
    case("theta_s")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%theta_s))
!       write(*,*) "setting theta_s to", x_1d_double_ptr
       s%theta_s(:) = x_1d_double_ptr
    case("ntheta")
       call c_f_pointer(xptr, x_0d_int_ptr, shape(s%ntheta))
!       write(*,*) "setting ntheta to", x_0d_int_ptr
       s%ntheta = x_0d_int_ptr
    case("ntheta_s")
       call c_f_pointer(xptr, x_0d_int_ptr, shape(s%ntheta_s))
!       write(*,*) "setting ntheta_s to", x_0d_int_ptr
       s%ntheta_s = x_0d_int_ptr
    case("dtheta")
       call c_f_pointer(xptr, x_0d_double_ptr, shape(s%dtheta))
!       write(*,*) "setting dtheta to", x_0d_double_ptr
       s%dtheta = x_0d_double_ptr
    case("dtheta_s")
       call c_f_pointer(xptr, x_0d_double_ptr, shape(s%dtheta_s))
!       write(*,*) "setting dtheta_s to", x_0d_double_ptr
       s%dtheta_s = x_0d_double_ptr
    case("theta0")
       call c_f_pointer(xptr, x_0d_double_ptr, shape(s%theta0))
!       write(*,*) "setting theta0 to", x_0d_double_ptr
       s%theta0 = x_0d_double_ptr
    case("thetamean")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%thetamean))
!       write(*,*) "setting thetamean to", x_2d_double_ptr
       s%thetamean(:,:) = x_2d_double_ptr
    case("Fx")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Fx))
!       write(*,*) "setting Fx to", x_2d_double_ptr
       s%Fx(:,:) = x_2d_double_ptr
    case("Fy")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Fy))
!       write(*,*) "setting Fy to", x_2d_double_ptr
       s%Fy(:,:) = x_2d_double_ptr
    case("Sxy")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Sxy))
!       write(*,*) "setting Sxy to", x_2d_double_ptr
       s%Sxy(:,:) = x_2d_double_ptr
    case("Syy")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Syy))
!       write(*,*) "setting Syy to", x_2d_double_ptr
       s%Syy(:,:) = x_2d_double_ptr
    case("Sxx")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Sxx))
!       write(*,*) "setting Sxx to", x_2d_double_ptr
       s%Sxx(:,:) = x_2d_double_ptr
    case("n")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%n))
!       write(*,*) "setting n to", x_2d_double_ptr
       s%n(:,:) = x_2d_double_ptr
    case("H")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%H))
!       write(*,*) "setting H to", x_2d_double_ptr
       s%H(:,:) = x_2d_double_ptr
    case("cgx")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%cgx))
!       write(*,*) "setting cgx to", x_3d_double_ptr
       s%cgx(:,:,:) = x_3d_double_ptr
    case("cgy")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%cgy))
!       write(*,*) "setting cgy to", x_3d_double_ptr
       s%cgy(:,:,:) = x_3d_double_ptr
    case("cx")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%cx))
!       write(*,*) "setting cx to", x_3d_double_ptr
       s%cx(:,:,:) = x_3d_double_ptr
    case("cy")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%cy))
!       write(*,*) "setting cy to", x_3d_double_ptr
       s%cy(:,:,:) = x_3d_double_ptr
    case("ctheta")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%ctheta))
!       write(*,*) "setting ctheta to", x_3d_double_ptr
       s%ctheta(:,:,:) = x_3d_double_ptr
    case("ee")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%ee))
!       write(*,*) "setting ee to", x_3d_double_ptr
       s%ee(:,:,:) = x_3d_double_ptr
    case("thet")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%thet))
!       write(*,*) "setting thet to", x_3d_double_ptr
       s%thet(:,:,:) = x_3d_double_ptr
    case("costh")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%costh))
!       write(*,*) "setting costh to", x_3d_double_ptr
       s%costh(:,:,:) = x_3d_double_ptr
    case("sinth")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%sinth))
!       write(*,*) "setting sinth to", x_3d_double_ptr
       s%sinth(:,:,:) = x_3d_double_ptr
    case("sigt")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%sigt))
!       write(*,*) "setting sigt to", x_3d_double_ptr
       s%sigt(:,:,:) = x_3d_double_ptr
    case("rr")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%rr))
!       write(*,*) "setting rr to", x_3d_double_ptr
       s%rr(:,:,:) = x_3d_double_ptr
    case("cgx_s")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%cgx_s))
!       write(*,*) "setting cgx_s to", x_3d_double_ptr
       s%cgx_s(:,:,:) = x_3d_double_ptr
    case("cgy_s")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%cgy_s))
!       write(*,*) "setting cgy_s to", x_3d_double_ptr
       s%cgy_s(:,:,:) = x_3d_double_ptr
    case("ctheta_s")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%ctheta_s))
!       write(*,*) "setting ctheta_s to", x_3d_double_ptr
       s%ctheta_s(:,:,:) = x_3d_double_ptr
    case("ee_s")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%ee_s))
!       write(*,*) "setting ee_s to", x_3d_double_ptr
       s%ee_s(:,:,:) = x_3d_double_ptr
    case("thet_s")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%thet_s))
!       write(*,*) "setting thet_s to", x_3d_double_ptr
       s%thet_s(:,:,:) = x_3d_double_ptr
    case("costh_s")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%costh_s))
!       write(*,*) "setting costh_s to", x_3d_double_ptr
       s%costh_s(:,:,:) = x_3d_double_ptr
    case("sinth_s")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%sinth_s))
!       write(*,*) "setting sinth_s to", x_3d_double_ptr
       s%sinth_s(:,:,:) = x_3d_double_ptr
    case("k")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%k))
!       write(*,*) "setting k to", x_2d_double_ptr
       s%k(:,:) = x_2d_double_ptr
    case("c")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%c))
!       write(*,*) "setting c to", x_2d_double_ptr
       s%c(:,:) = x_2d_double_ptr
    case("cg")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%cg))
!       write(*,*) "setting cg to", x_2d_double_ptr
       s%cg(:,:) = x_2d_double_ptr
    case("sigm")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%sigm))
!       write(*,*) "setting sigm to", x_2d_double_ptr
       s%sigm(:,:) = x_2d_double_ptr
    case("wm")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%wm))
!       write(*,*) "setting wm to", x_2d_double_ptr
       s%wm(:,:) = x_2d_double_ptr
    case("hh")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%hh))
!       write(*,*) "setting hh to", x_2d_double_ptr
       s%hh(:,:) = x_2d_double_ptr
    case("zs")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%zs))
!       write(*,*) "setting zs to", x_2d_double_ptr
       s%zs(:,:) = x_2d_double_ptr
    case("zs0")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%zs0))
!       write(*,*) "setting zs0 to", x_2d_double_ptr
       s%zs0(:,:) = x_2d_double_ptr
    case("tideinpt")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%tideinpt))
!       write(*,*) "setting tideinpt to", x_1d_double_ptr
       s%tideinpt(:) = x_1d_double_ptr
    case("tideinpz")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%tideinpz))
!       write(*,*) "setting tideinpz to", x_2d_double_ptr
       s%tideinpz(:,:) = x_2d_double_ptr
    case("windinpt")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%windinpt))
!       write(*,*) "setting windinpt to", x_1d_double_ptr
       s%windinpt(:) = x_1d_double_ptr
    case("windvelts")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%windvelts))
!       write(*,*) "setting windvelts to", x_1d_double_ptr
       s%windvelts(:) = x_1d_double_ptr
    case("winddirts")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%winddirts))
!       write(*,*) "setting winddirts to", x_1d_double_ptr
       s%winddirts(:) = x_1d_double_ptr
    case("windxts")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%windxts))
!       write(*,*) "setting windxts to", x_1d_double_ptr
       s%windxts(:) = x_1d_double_ptr
    case("windyts")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%windyts))
!       write(*,*) "setting windyts to", x_1d_double_ptr
       s%windyts(:) = x_1d_double_ptr
    case("windsu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%windsu))
!       write(*,*) "setting windsu to", x_2d_double_ptr
       s%windsu(:,:) = x_2d_double_ptr
    case("windnv")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%windnv))
!       write(*,*) "setting windnv to", x_2d_double_ptr
       s%windnv(:,:) = x_2d_double_ptr
    case("dzsdt")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dzsdt))
!       write(*,*) "setting dzsdt to", x_2d_double_ptr
       s%dzsdt(:,:) = x_2d_double_ptr
    case("dzsdx")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dzsdx))
!       write(*,*) "setting dzsdx to", x_2d_double_ptr
       s%dzsdx(:,:) = x_2d_double_ptr
    case("dzsdy")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dzsdy))
!       write(*,*) "setting dzsdy to", x_2d_double_ptr
       s%dzsdy(:,:) = x_2d_double_ptr
    case("dzbdx")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dzbdx))
!       write(*,*) "setting dzbdx to", x_2d_double_ptr
       s%dzbdx(:,:) = x_2d_double_ptr
    case("dzbdy")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dzbdy))
!       write(*,*) "setting dzbdy to", x_2d_double_ptr
       s%dzbdy(:,:) = x_2d_double_ptr
    case("dzbdt")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dzbdt))
!       write(*,*) "setting dzbdt to", x_2d_double_ptr
       s%dzbdt(:,:) = x_2d_double_ptr
    case("dzbnow")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dzbnow))
!       write(*,*) "setting dzbnow to", x_2d_double_ptr
       s%dzbnow(:,:) = x_2d_double_ptr
    case("uu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%uu))
!       write(*,*) "setting uu to", x_2d_double_ptr
       s%uu(:,:) = x_2d_double_ptr
    case("vu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vu))
!       write(*,*) "setting vu to", x_2d_double_ptr
       s%vu(:,:) = x_2d_double_ptr
    case("uv")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%uv))
!       write(*,*) "setting uv to", x_2d_double_ptr
       s%uv(:,:) = x_2d_double_ptr
    case("vv")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vv))
!       write(*,*) "setting vv to", x_2d_double_ptr
       s%vv(:,:) = x_2d_double_ptr
    case("qx")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%qx))
!       write(*,*) "setting qx to", x_2d_double_ptr
       s%qx(:,:) = x_2d_double_ptr
    case("qy")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%qy))
!       write(*,*) "setting qy to", x_2d_double_ptr
       s%qy(:,:) = x_2d_double_ptr
    case("sedero")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%sedero))
!       write(*,*) "setting sedero to", x_2d_double_ptr
       s%sedero(:,:) = x_2d_double_ptr
    case("dcbdx")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dcbdx))
!       write(*,*) "setting dcbdx to", x_2d_double_ptr
       s%dcbdx(:,:) = x_2d_double_ptr
    case("dcbdy")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dcbdy))
!       write(*,*) "setting dcbdy to", x_2d_double_ptr
       s%dcbdy(:,:) = x_2d_double_ptr
    case("dcsdx")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dcsdx))
!       write(*,*) "setting dcsdx to", x_2d_double_ptr
       s%dcsdx(:,:) = x_2d_double_ptr
    case("dcsdy")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dcsdy))
!       write(*,*) "setting dcsdy to", x_2d_double_ptr
       s%dcsdy(:,:) = x_2d_double_ptr
    case("ui")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%ui))
!       write(*,*) "setting ui to", x_2d_double_ptr
       s%ui(:,:) = x_2d_double_ptr
    case("vi")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vi))
!       write(*,*) "setting vi to", x_2d_double_ptr
       s%vi(:,:) = x_2d_double_ptr
    case("E")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%E))
!       write(*,*) "setting E to", x_2d_double_ptr
       s%E(:,:) = x_2d_double_ptr
    case("R")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%R))
!       write(*,*) "setting R to", x_2d_double_ptr
       s%R(:,:) = x_2d_double_ptr
    case("urms")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%urms))
!       write(*,*) "setting urms to", x_2d_double_ptr
       s%urms(:,:) = x_2d_double_ptr
    case("D")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%D))
!       write(*,*) "setting D to", x_2d_double_ptr
       s%D(:,:) = x_2d_double_ptr
    case("Qb")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Qb))
!       write(*,*) "setting Qb to", x_2d_double_ptr
       s%Qb(:,:) = x_2d_double_ptr
    case("ust")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%ust))
!       write(*,*) "setting ust to", x_2d_double_ptr
       s%ust(:,:) = x_2d_double_ptr
    case("ueu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%ueu))
!       write(*,*) "setting ueu to", x_2d_double_ptr
       s%ueu(:,:) = x_2d_double_ptr
    case("vev")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vev))
!       write(*,*) "setting vev to", x_2d_double_ptr
       s%vev(:,:) = x_2d_double_ptr
    case("vmagu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vmagu))
!       write(*,*) "setting vmagu to", x_2d_double_ptr
       s%vmagu(:,:) = x_2d_double_ptr
    case("vmageu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vmageu))
!       write(*,*) "setting vmageu to", x_2d_double_ptr
       s%vmageu(:,:) = x_2d_double_ptr
    case("vmagv")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vmagv))
!       write(*,*) "setting vmagv to", x_2d_double_ptr
       s%vmagv(:,:) = x_2d_double_ptr
    case("vmagev")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vmagev))
!       write(*,*) "setting vmagev to", x_2d_double_ptr
       s%vmagev(:,:) = x_2d_double_ptr
    case("u")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%u))
!       write(*,*) "setting u to", x_2d_double_ptr
       s%u(:,:) = x_2d_double_ptr
    case("v")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%v))
!       write(*,*) "setting v to", x_2d_double_ptr
       s%v(:,:) = x_2d_double_ptr
    case("ue")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%ue))
!       write(*,*) "setting ue to", x_2d_double_ptr
       s%ue(:,:) = x_2d_double_ptr
    case("ve")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%ve))
!       write(*,*) "setting ve to", x_2d_double_ptr
       s%ve(:,:) = x_2d_double_ptr
    case("hold")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%hold))
!       write(*,*) "setting hold to", x_2d_double_ptr
       s%hold(:,:) = x_2d_double_ptr
    case("wetu")
       call c_f_pointer(xptr, x_2d_int_ptr, shape(s%wetu))
!       write(*,*) "setting wetu to", x_2d_int_ptr
       s%wetu(:,:) = x_2d_int_ptr
    case("wetv")
       call c_f_pointer(xptr, x_2d_int_ptr, shape(s%wetv))
!       write(*,*) "setting wetv to", x_2d_int_ptr
       s%wetv(:,:) = x_2d_int_ptr
    case("wetz")
       call c_f_pointer(xptr, x_2d_int_ptr, shape(s%wetz))
!       write(*,*) "setting wetz to", x_2d_int_ptr
       s%wetz(:,:) = x_2d_int_ptr
    case("hu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%hu))
!       write(*,*) "setting hu to", x_2d_double_ptr
       s%hu(:,:) = x_2d_double_ptr
    case("hv")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%hv))
!       write(*,*) "setting hv to", x_2d_double_ptr
       s%hv(:,:) = x_2d_double_ptr
    case("hum")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%hum))
!       write(*,*) "setting hum to", x_2d_double_ptr
       s%hum(:,:) = x_2d_double_ptr
    case("hvm")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%hvm))
!       write(*,*) "setting hvm to", x_2d_double_ptr
       s%hvm(:,:) = x_2d_double_ptr
    case("vmag")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vmag))
!       write(*,*) "setting vmag to", x_2d_double_ptr
       s%vmag(:,:) = x_2d_double_ptr
    case("ccg")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%ccg))
!       write(*,*) "setting ccg to", x_3d_double_ptr
       s%ccg(:,:,:) = x_3d_double_ptr
    case("uwf")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%uwf))
!       write(*,*) "setting uwf to", x_2d_double_ptr
       s%uwf(:,:) = x_2d_double_ptr
    case("vwf")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vwf))
!       write(*,*) "setting vwf to", x_2d_double_ptr
       s%vwf(:,:) = x_2d_double_ptr
    case("ustr")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%ustr))
!       write(*,*) "setting ustr to", x_2d_double_ptr
       s%ustr(:,:) = x_2d_double_ptr
    case("usd")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%usd))
!       write(*,*) "setting usd to", x_2d_double_ptr
       s%usd(:,:) = x_2d_double_ptr
    case("bi")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%bi))
!       write(*,*) "setting bi to", x_1d_double_ptr
       s%bi(:) = x_1d_double_ptr
    case("DR")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%DR))
!       write(*,*) "setting DR to", x_2d_double_ptr
       s%DR(:,:) = x_2d_double_ptr
    case("umean")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%umean))
!       write(*,*) "setting umean to", x_2d_double_ptr
       s%umean(:,:) = x_2d_double_ptr
    case("vmean")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vmean))
!       write(*,*) "setting vmean to", x_2d_double_ptr
       s%vmean(:,:) = x_2d_double_ptr
    case("ur")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%ur))
!       write(*,*) "setting ur to", x_2d_double_ptr
       s%ur(:,:) = x_2d_double_ptr
    case("vardx")
       call c_f_pointer(xptr, x_0d_int_ptr, shape(s%vardx))
!       write(*,*) "setting vardx to", x_0d_int_ptr
       s%vardx = x_0d_int_ptr
    case("D15")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%D15))
!       write(*,*) "setting D15 to", x_1d_double_ptr
       s%D15(:) = x_1d_double_ptr
    case("D50")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%D50))
!       write(*,*) "setting D50 to", x_1d_double_ptr
       s%D50(:) = x_1d_double_ptr
    case("D90")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%D90))
!       write(*,*) "setting D90 to", x_1d_double_ptr
       s%D90(:) = x_1d_double_ptr
    case("sedcal")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%sedcal))
!       write(*,*) "setting sedcal to", x_1d_double_ptr
       s%sedcal(:) = x_1d_double_ptr
    case("ucrcal")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%ucrcal))
!       write(*,*) "setting ucrcal to", x_1d_double_ptr
       s%ucrcal(:) = x_1d_double_ptr
    case("Tsg")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%Tsg))
!       write(*,*) "setting Tsg to", x_3d_double_ptr
       s%Tsg(:,:,:) = x_3d_double_ptr
    case("Susg")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%Susg))
!       write(*,*) "setting Susg to", x_3d_double_ptr
       s%Susg(:,:,:) = x_3d_double_ptr
    case("Svsg")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%Svsg))
!       write(*,*) "setting Svsg to", x_3d_double_ptr
       s%Svsg(:,:,:) = x_3d_double_ptr
    case("Subg")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%Subg))
!       write(*,*) "setting Subg to", x_3d_double_ptr
       s%Subg(:,:,:) = x_3d_double_ptr
    case("Svbg")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%Svbg))
!       write(*,*) "setting Svbg to", x_3d_double_ptr
       s%Svbg(:,:,:) = x_3d_double_ptr
    case("ceqbg")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%ceqbg))
!       write(*,*) "setting ceqbg to", x_3d_double_ptr
       s%ceqbg(:,:,:) = x_3d_double_ptr
    case("ceqsg")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%ceqsg))
!       write(*,*) "setting ceqsg to", x_3d_double_ptr
       s%ceqsg(:,:,:) = x_3d_double_ptr
    case("ua")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%ua))
!       write(*,*) "setting ua to", x_2d_double_ptr
       s%ua(:,:) = x_2d_double_ptr
    case("BR")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%BR))
!       write(*,*) "setting BR to", x_2d_double_ptr
       s%BR(:,:) = x_2d_double_ptr
    case("kb")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%kb))
!       write(*,*) "setting kb to", x_2d_double_ptr
       s%kb(:,:) = x_2d_double_ptr
    case("Tbore")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Tbore))
!       write(*,*) "setting Tbore to", x_2d_double_ptr
       s%Tbore(:,:) = x_2d_double_ptr
    case("dzav")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dzav))
!       write(*,*) "setting dzav to", x_2d_double_ptr
       s%dzav(:,:) = x_2d_double_ptr
    case("maxzs")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%maxzs))
!       write(*,*) "setting maxzs to", x_2d_double_ptr
       s%maxzs(:,:) = x_2d_double_ptr
    case("minzs")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%minzs))
!       write(*,*) "setting minzs to", x_2d_double_ptr
       s%minzs(:,:) = x_2d_double_ptr
    case("L1")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%L1))
!       write(*,*) "setting L1 to", x_2d_double_ptr
       s%L1(:,:) = x_2d_double_ptr
    case("Sk")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Sk))
!       write(*,*) "setting Sk to", x_2d_double_ptr
       s%Sk(:,:) = x_2d_double_ptr
    case("As")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%As))
!       write(*,*) "setting As to", x_2d_double_ptr
       s%As(:,:) = x_2d_double_ptr
    case("gwhead")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%gwhead))
!       write(*,*) "setting gwhead to", x_2d_double_ptr
       s%gwhead(:,:) = x_2d_double_ptr
    case("gwheadb")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%gwheadb))
!       write(*,*) "setting gwheadb to", x_2d_double_ptr
       s%gwheadb(:,:) = x_2d_double_ptr
    case("gwlevel")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%gwlevel))
!       write(*,*) "setting gwlevel to", x_2d_double_ptr
       s%gwlevel(:,:) = x_2d_double_ptr
    case("gwheight")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%gwheight))
!       write(*,*) "setting gwheight to", x_2d_double_ptr
       s%gwheight(:,:) = x_2d_double_ptr
    case("gwbottom")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%gwbottom))
!       write(*,*) "setting gwbottom to", x_2d_double_ptr
       s%gwbottom(:,:) = x_2d_double_ptr
    case("gwu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%gwu))
!       write(*,*) "setting gwu to", x_2d_double_ptr
       s%gwu(:,:) = x_2d_double_ptr
    case("gwv")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%gwv))
!       write(*,*) "setting gwv to", x_2d_double_ptr
       s%gwv(:,:) = x_2d_double_ptr
    case("gwqx")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%gwqx))
!       write(*,*) "setting gwqx to", x_2d_double_ptr
       s%gwqx(:,:) = x_2d_double_ptr
    case("gwqy")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%gwqy))
!       write(*,*) "setting gwqy to", x_2d_double_ptr
       s%gwqy(:,:) = x_2d_double_ptr
    case("gww")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%gww))
!       write(*,*) "setting gww to", x_2d_double_ptr
       s%gww(:,:) = x_2d_double_ptr
    case("gwcurv")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%gwcurv))
!       write(*,*) "setting gwcurv to", x_2d_double_ptr
       s%gwcurv(:,:) = x_2d_double_ptr
    case("dinfil")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dinfil))
!       write(*,*) "setting dinfil to", x_2d_double_ptr
       s%dinfil(:,:) = x_2d_double_ptr
    case("infil")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%infil))
!       write(*,*) "setting infil to", x_2d_double_ptr
       s%infil(:,:) = x_2d_double_ptr
    case("gw0back")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%gw0back))
!       write(*,*) "setting gw0back to", x_2d_double_ptr
       s%gw0back(:,:) = x_2d_double_ptr
    case("kturb")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%kturb))
!       write(*,*) "setting kturb to", x_2d_double_ptr
       s%kturb(:,:) = x_2d_double_ptr
    case("ero")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%ero))
!       write(*,*) "setting ero to", x_3d_double_ptr
       s%ero(:,:,:) = x_3d_double_ptr
    case("depo_im")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%depo_im))
!       write(*,*) "setting depo_im to", x_3d_double_ptr
       s%depo_im(:,:,:) = x_3d_double_ptr
    case("depo_ex")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%depo_ex))
!       write(*,*) "setting depo_ex to", x_3d_double_ptr
       s%depo_ex(:,:,:) = x_3d_double_ptr
    case("nd")
       call c_f_pointer(xptr, x_2d_int_ptr, shape(s%nd))
!       write(*,*) "setting nd to", x_2d_int_ptr
       s%nd(:,:) = x_2d_int_ptr
    case("pbbed")
       call c_f_pointer(xptr, x_4d_double_ptr, shape(s%pbbed))
!       write(*,*) "setting pbbed to", x_4d_double_ptr
       s%pbbed(:,:,:,:) = x_4d_double_ptr
    case("dzbed")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%dzbed))
!       write(*,*) "setting dzbed to", x_3d_double_ptr
       s%dzbed(:,:,:) = x_3d_double_ptr
    case("z0bed")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%z0bed))
!       write(*,*) "setting z0bed to", x_2d_double_ptr
       s%z0bed(:,:) = x_2d_double_ptr
    case("ureps")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%ureps))
!       write(*,*) "setting ureps to", x_2d_double_ptr
       s%ureps(:,:) = x_2d_double_ptr
    case("vreps")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vreps))
!       write(*,*) "setting vreps to", x_2d_double_ptr
       s%vreps(:,:) = x_2d_double_ptr
    case("urepb")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%urepb))
!       write(*,*) "setting urepb to", x_2d_double_ptr
       s%urepb(:,:) = x_2d_double_ptr
    case("vrepb")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vrepb))
!       write(*,*) "setting vrepb to", x_2d_double_ptr
       s%vrepb(:,:) = x_2d_double_ptr
    case("umwci")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%umwci))
!       write(*,*) "setting umwci to", x_2d_double_ptr
       s%umwci(:,:) = x_2d_double_ptr
    case("vmwci")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vmwci))
!       write(*,*) "setting vmwci to", x_2d_double_ptr
       s%vmwci(:,:) = x_2d_double_ptr
    case("rolthick")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%rolthick))
!       write(*,*) "setting rolthick to", x_2d_double_ptr
       s%rolthick(:,:) = x_2d_double_ptr
    case("zswci")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%zswci))
!       write(*,*) "setting zswci to", x_2d_double_ptr
       s%zswci(:,:) = x_2d_double_ptr
    case("pres")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%pres))
!       write(*,*) "setting pres to", x_2d_double_ptr
       s%pres(:,:) = x_2d_double_ptr
    case("dU")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dU))
!       write(*,*) "setting dU to", x_2d_double_ptr
       s%dU(:,:) = x_2d_double_ptr
    case("dV")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dV))
!       write(*,*) "setting dV to", x_2d_double_ptr
       s%dV(:,:) = x_2d_double_ptr
    case("wb")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%wb))
!       write(*,*) "setting wb to", x_2d_double_ptr
       s%wb(:,:) = x_2d_double_ptr
    case("ws")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%ws))
!       write(*,*) "setting ws to", x_2d_double_ptr
       s%ws(:,:) = x_2d_double_ptr
    case("wscrit")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%wscrit))
!       write(*,*) "setting wscrit to", x_2d_double_ptr
       s%wscrit(:,:) = x_2d_double_ptr
    case("bedfriccoef")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%bedfriccoef))
!       write(*,*) "setting bedfriccoef to", x_2d_double_ptr
       s%bedfriccoef(:,:) = x_2d_double_ptr
    case("taubx")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%taubx))
!       write(*,*) "setting taubx to", x_2d_double_ptr
       s%taubx(:,:) = x_2d_double_ptr
    case("tauby")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%tauby))
!       write(*,*) "setting tauby to", x_2d_double_ptr
       s%tauby(:,:) = x_2d_double_ptr
    case("Df")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Df))
!       write(*,*) "setting Df to", x_2d_double_ptr
       s%Df(:,:) = x_2d_double_ptr
    case("Dp")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Dp))
!       write(*,*) "setting Dp to", x_2d_double_ptr
       s%Dp(:,:) = x_2d_double_ptr
    case("Sutot")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Sutot))
!       write(*,*) "setting Sutot to", x_2d_double_ptr
       s%Sutot(:,:) = x_2d_double_ptr
    case("Svtot")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Svtot))
!       write(*,*) "setting Svtot to", x_2d_double_ptr
       s%Svtot(:,:) = x_2d_double_ptr
    case("cctot")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%cctot))
!       write(*,*) "setting cctot to", x_2d_double_ptr
       s%cctot(:,:) = x_2d_double_ptr
    case("wi")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%wi))
!       write(*,*) "setting wi to", x_2d_double_ptr
       s%wi(:,:) = x_2d_double_ptr
    case("dUi")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dUi))
!       write(*,*) "setting dUi to", x_2d_double_ptr
       s%dUi(:,:) = x_2d_double_ptr
    case("zi")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%zi))
!       write(*,*) "setting zi to", x_2d_double_ptr
       s%zi(:,:) = x_2d_double_ptr
    case("nuh")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%nuh))
!       write(*,*) "setting nuh to", x_2d_double_ptr
       s%nuh(:,:) = x_2d_double_ptr
    case("cf")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%cf))
!       write(*,*) "setting cf to", x_2d_double_ptr
       s%cf(:,:) = x_2d_double_ptr
    case("cfu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%cfu))
!       write(*,*) "setting cfu to", x_2d_double_ptr
       s%cfu(:,:) = x_2d_double_ptr
    case("cfv")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%cfv))
!       write(*,*) "setting cfv to", x_2d_double_ptr
       s%cfv(:,:) = x_2d_double_ptr
    case("D50top")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%D50top))
!       write(*,*) "setting D50top to", x_2d_double_ptr
       s%D50top(:,:) = x_2d_double_ptr
    case("D90top")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%D90top))
!       write(*,*) "setting D90top to", x_2d_double_ptr
       s%D90top(:,:) = x_2d_double_ptr
    case("structdepth")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%structdepth))
!       write(*,*) "setting structdepth to", x_2d_double_ptr
       s%structdepth(:,:) = x_2d_double_ptr
    case("zs0fac")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%zs0fac))
!       write(*,*) "setting zs0fac to", x_3d_double_ptr
       s%zs0fac(:,:,:) = x_3d_double_ptr
    case("tdisch")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%tdisch))
!       write(*,*) "setting tdisch to", x_1d_double_ptr
       s%tdisch(:) = x_1d_double_ptr
    case("pdisch")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%pdisch))
!       write(*,*) "setting pdisch to", x_2d_double_ptr
       s%pdisch(:,:) = x_2d_double_ptr
    case("pntdisch")
       call c_f_pointer(xptr, x_1d_int_ptr, shape(s%pntdisch))
!       write(*,*) "setting pntdisch to", x_1d_int_ptr
       s%pntdisch(:) = x_1d_int_ptr
    case("qdisch")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%qdisch))
!       write(*,*) "setting qdisch to", x_2d_double_ptr
       s%qdisch(:,:) = x_2d_double_ptr
    case("idrift")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%idrift))
!       write(*,*) "setting idrift to", x_1d_double_ptr
       s%idrift(:) = x_1d_double_ptr
    case("jdrift")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%jdrift))
!       write(*,*) "setting jdrift to", x_1d_double_ptr
       s%jdrift(:) = x_1d_double_ptr
    case("tdriftb")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%tdriftb))
!       write(*,*) "setting tdriftb to", x_1d_double_ptr
       s%tdriftb(:) = x_1d_double_ptr
    case("tdrifte")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%tdrifte))
!       write(*,*) "setting tdrifte to", x_1d_double_ptr
       s%tdrifte(:) = x_1d_double_ptr
    case("runup")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%runup))
!       write(*,*) "setting runup to", x_1d_double_ptr
       s%runup(:) = x_1d_double_ptr
    case("Hrunup")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%Hrunup))
!       write(*,*) "setting Hrunup to", x_1d_double_ptr
       s%Hrunup(:) = x_1d_double_ptr
    case("xHrunup")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%xHrunup))
!       write(*,*) "setting xHrunup to", x_1d_double_ptr
       s%xHrunup(:) = x_1d_double_ptr
    case("istruct")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%istruct))
!       write(*,*) "setting istruct to", x_1d_double_ptr
       s%istruct(:) = x_1d_double_ptr
    case("iwl")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%iwl))
!       write(*,*) "setting iwl to", x_1d_double_ptr
       s%iwl(:) = x_1d_double_ptr
    case("strucslope")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%strucslope))
!       write(*,*) "setting strucslope to", x_1d_double_ptr
       s%strucslope(:) = x_1d_double_ptr
    case("Dc")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Dc))
!       write(*,*) "setting Dc to", x_2d_double_ptr
       s%Dc(:,:) = x_2d_double_ptr
    case("ph")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%ph))
!       write(*,*) "setting ph to", x_2d_double_ptr
       s%ph(:,:) = x_2d_double_ptr
    case("newstatbc")
       call c_f_pointer(xptr, x_0d_int_ptr, shape(s%newstatbc))
!       write(*,*) "setting newstatbc to", x_0d_int_ptr
       s%newstatbc = x_0d_int_ptr
    case("dobs")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dobs))
!       write(*,*) "setting dobs to", x_2d_double_ptr
       s%dobs(:,:) = x_2d_double_ptr
    case("sig2prior")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%sig2prior))
!       write(*,*) "setting sig2prior to", x_2d_double_ptr
       s%sig2prior(:,:) = x_2d_double_ptr
    case("zbobs")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%zbobs))
!       write(*,*) "setting zbobs to", x_2d_double_ptr
       s%zbobs(:,:) = x_2d_double_ptr
    case("shobs")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%shobs))
!       write(*,*) "setting shobs to", x_2d_double_ptr
       s%shobs(:,:) = x_2d_double_ptr
    case("bwalpha")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%bwalpha))
!       write(*,*) "setting bwalpha to", x_2d_double_ptr
       s%bwalpha(:,:) = x_2d_double_ptr
    case("dcmdo")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dcmdo))
!       write(*,*) "setting dcmdo to", x_2d_double_ptr
       s%dcmdo(:,:) = x_2d_double_ptr
    case("dassim")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dassim))
!       write(*,*) "setting dassim to", x_2d_double_ptr
       s%dassim(:,:) = x_2d_double_ptr
    case("cobs")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%cobs))
!       write(*,*) "setting cobs to", x_2d_double_ptr
       s%cobs(:,:) = x_2d_double_ptr
    case("shipxCG")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%shipxCG))
!       write(*,*) "setting shipxCG to", x_1d_double_ptr
       s%shipxCG(:) = x_1d_double_ptr
    case("shipyCG")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%shipyCG))
!       write(*,*) "setting shipyCG to", x_1d_double_ptr
       s%shipyCG(:) = x_1d_double_ptr
    case("shipzCG")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%shipzCG))
!       write(*,*) "setting shipzCG to", x_1d_double_ptr
       s%shipzCG(:) = x_1d_double_ptr
    case("shipFx")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%shipFx))
!       write(*,*) "setting shipFx to", x_1d_double_ptr
       s%shipFx(:) = x_1d_double_ptr
    case("shipFy")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%shipFy))
!       write(*,*) "setting shipFy to", x_1d_double_ptr
       s%shipFy(:) = x_1d_double_ptr
    case("shipFz")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%shipFz))
!       write(*,*) "setting shipFz to", x_1d_double_ptr
       s%shipFz(:) = x_1d_double_ptr
    case("shipMx")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%shipMx))
!       write(*,*) "setting shipMx to", x_1d_double_ptr
       s%shipMx(:) = x_1d_double_ptr
    case("shipMy")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%shipMy))
!       write(*,*) "setting shipMy to", x_1d_double_ptr
       s%shipMy(:) = x_1d_double_ptr
    case("shipMz")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%shipMz))
!       write(*,*) "setting shipMz to", x_1d_double_ptr
       s%shipMz(:) = x_1d_double_ptr
    case("shipphi")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%shipphi))
!       write(*,*) "setting shipphi to", x_1d_double_ptr
       s%shipphi(:) = x_1d_double_ptr
    case("shipchi")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%shipchi))
!       write(*,*) "setting shipchi to", x_1d_double_ptr
       s%shipchi(:) = x_1d_double_ptr
    case("shippsi")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%shippsi))
!       write(*,*) "setting shippsi to", x_1d_double_ptr
       s%shippsi(:) = x_1d_double_ptr
    case("Dveg")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Dveg))
!       write(*,*) "setting Dveg to", x_2d_double_ptr
       s%Dveg(:,:) = x_2d_double_ptr
    case("Fvegu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Fvegu))
!       write(*,*) "setting Fvegu to", x_2d_double_ptr
       s%Fvegu(:,:) = x_2d_double_ptr
    case("Fvegv")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Fvegv))
!       write(*,*) "setting Fvegv to", x_2d_double_ptr
       s%Fvegv(:,:) = x_2d_double_ptr
    case("ududx")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%ududx))
!       write(*,*) "setting ududx to", x_2d_double_ptr
       s%ududx(:,:) = x_2d_double_ptr
    case("vdvdy")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vdvdy))
!       write(*,*) "setting vdvdy to", x_2d_double_ptr
       s%vdvdy(:,:) = x_2d_double_ptr
    case("udvdx")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%udvdx))
!       write(*,*) "setting udvdx to", x_2d_double_ptr
       s%udvdx(:,:) = x_2d_double_ptr
    case("vdudy")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vdudy))
!       write(*,*) "setting vdudy to", x_2d_double_ptr
       s%vdudy(:,:) = x_2d_double_ptr
    case("viscu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%viscu))
!       write(*,*) "setting viscu to", x_2d_double_ptr
       s%viscu(:,:) = x_2d_double_ptr
    case("viscv")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%viscv))
!       write(*,*) "setting viscv to", x_2d_double_ptr
       s%viscv(:,:) = x_2d_double_ptr
    case("setbathy")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%setbathy))
!       write(*,*) "setting setbathy to", x_3d_double_ptr
       s%setbathy(:,:,:) = x_3d_double_ptr
    case("tsetbathy")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%tsetbathy))
!       write(*,*) "setting tsetbathy to", x_1d_double_ptr
       s%tsetbathy(:) = x_1d_double_ptr
    case("breaking")
       call c_f_pointer(xptr, x_2d_int_ptr, shape(s%breaking))
!       write(*,*) "setting breaking to", x_2d_int_ptr
       s%breaking(:,:) = x_2d_int_ptr
    case default
       write(*,*) 'Setting unknown variable', var_name
    end select

  end subroutine set_var
