

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
       s%x(:,:) = x_2d_double_ptr
    case("y")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%y))
       s%y(:,:) = x_2d_double_ptr
    case("xz")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%xz))
       s%xz(:,:) = x_2d_double_ptr
    case("yz")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%yz))
       s%yz(:,:) = x_2d_double_ptr
    case("xu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%xu))
       s%xu(:,:) = x_2d_double_ptr
    case("yu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%yu))
       s%yu(:,:) = x_2d_double_ptr
    case("xv")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%xv))
       s%xv(:,:) = x_2d_double_ptr
    case("yv")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%yv))
       s%yv(:,:) = x_2d_double_ptr
    case("dsu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dsu))
       s%dsu(:,:) = x_2d_double_ptr
    case("dsv")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dsv))
       s%dsv(:,:) = x_2d_double_ptr
    case("dsz")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dsz))
       s%dsz(:,:) = x_2d_double_ptr
    case("dsc")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dsc))
       s%dsc(:,:) = x_2d_double_ptr
    case("dnu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dnu))
       s%dnu(:,:) = x_2d_double_ptr
    case("dnv")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dnv))
       s%dnv(:,:) = x_2d_double_ptr
    case("dnz")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dnz))
       s%dnz(:,:) = x_2d_double_ptr
    case("dnc")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dnc))
       s%dnc(:,:) = x_2d_double_ptr
    case("dsdnui")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dsdnui))
       s%dsdnui(:,:) = x_2d_double_ptr
    case("dsdnvi")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dsdnvi))
       s%dsdnvi(:,:) = x_2d_double_ptr
    case("dsdnzi")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dsdnzi))
       s%dsdnzi(:,:) = x_2d_double_ptr
    case("alfaz")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%alfaz))
       s%alfaz(:,:) = x_2d_double_ptr
    case("alfau")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%alfau))
       s%alfau(:,:) = x_2d_double_ptr
    case("alfav")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%alfav))
       s%alfav(:,:) = x_2d_double_ptr
    case("sdist")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%sdist))
       s%sdist(:,:) = x_2d_double_ptr
    case("ndist")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%ndist))
       s%ndist(:,:) = x_2d_double_ptr
    case("dx")
       call c_f_pointer(xptr, x_0d_double_ptr) ! skip shape in case of rank 0
       s%dx = x_0d_double_ptr
    case("dy")
       call c_f_pointer(xptr, x_0d_double_ptr) ! skip shape in case of rank 0
       s%dy = x_0d_double_ptr
    case("xori")
       call c_f_pointer(xptr, x_0d_double_ptr) ! skip shape in case of rank 0
       s%xori = x_0d_double_ptr
    case("yori")
       call c_f_pointer(xptr, x_0d_double_ptr) ! skip shape in case of rank 0
       s%yori = x_0d_double_ptr
    case("alfa")
       call c_f_pointer(xptr, x_0d_double_ptr) ! skip shape in case of rank 0
       s%alfa = x_0d_double_ptr
    case("posdwn")
       call c_f_pointer(xptr, x_0d_double_ptr) ! skip shape in case of rank 0
       s%posdwn = x_0d_double_ptr
    case("thetamin")
       call c_f_pointer(xptr, x_0d_double_ptr) ! skip shape in case of rank 0
       s%thetamin = x_0d_double_ptr
    case("thetamax")
       call c_f_pointer(xptr, x_0d_double_ptr) ! skip shape in case of rank 0
       s%thetamax = x_0d_double_ptr
    case("nx")
       call c_f_pointer(xptr, x_0d_int_ptr) ! skip shape in case of rank 0
       s%nx = x_0d_int_ptr
    case("ny")
       call c_f_pointer(xptr, x_0d_int_ptr) ! skip shape in case of rank 0
       s%ny = x_0d_int_ptr
    case("zs01")
       call c_f_pointer(xptr, x_0d_double_ptr) ! skip shape in case of rank 0
       s%zs01 = x_0d_double_ptr
    case("zs02")
       call c_f_pointer(xptr, x_0d_double_ptr) ! skip shape in case of rank 0
       s%zs02 = x_0d_double_ptr
    case("zs03")
       call c_f_pointer(xptr, x_0d_double_ptr) ! skip shape in case of rank 0
       s%zs03 = x_0d_double_ptr
    case("zs04")
       call c_f_pointer(xptr, x_0d_double_ptr) ! skip shape in case of rank 0
       s%zs04 = x_0d_double_ptr
    case("xyzs01")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%xyzs01))
       s%xyzs01(:) = x_1d_double_ptr
    case("xyzs02")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%xyzs02))
       s%xyzs02(:) = x_1d_double_ptr
    case("xyzs03")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%xyzs03))
       s%xyzs03(:) = x_1d_double_ptr
    case("xyzs04")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%xyzs04))
       s%xyzs04(:) = x_1d_double_ptr
    case("tidelen")
       call c_f_pointer(xptr, x_0d_int_ptr) ! skip shape in case of rank 0
       s%tidelen = x_0d_int_ptr
    case("windlen")
       call c_f_pointer(xptr, x_0d_int_ptr) ! skip shape in case of rank 0
       s%windlen = x_0d_int_ptr
    case("zb")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%zb))
       s%zb(:,:) = x_2d_double_ptr
    case("zb0")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%zb0))
       s%zb0(:,:) = x_2d_double_ptr
    case("theta")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%theta))
       s%theta(:) = x_1d_double_ptr
    case("theta_s")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%theta_s))
       s%theta_s(:) = x_1d_double_ptr
    case("ntheta")
       call c_f_pointer(xptr, x_0d_int_ptr) ! skip shape in case of rank 0
       s%ntheta = x_0d_int_ptr
    case("ntheta_s")
       call c_f_pointer(xptr, x_0d_int_ptr) ! skip shape in case of rank 0
       s%ntheta_s = x_0d_int_ptr
    case("dtheta")
       call c_f_pointer(xptr, x_0d_double_ptr) ! skip shape in case of rank 0
       s%dtheta = x_0d_double_ptr
    case("dtheta_s")
       call c_f_pointer(xptr, x_0d_double_ptr) ! skip shape in case of rank 0
       s%dtheta_s = x_0d_double_ptr
    case("theta0")
       call c_f_pointer(xptr, x_0d_double_ptr) ! skip shape in case of rank 0
       s%theta0 = x_0d_double_ptr
    case("thetamean")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%thetamean))
       s%thetamean(:,:) = x_2d_double_ptr
    case("Fx")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Fx))
       s%Fx(:,:) = x_2d_double_ptr
    case("Fy")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Fy))
       s%Fy(:,:) = x_2d_double_ptr
    case("Sxy")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Sxy))
       s%Sxy(:,:) = x_2d_double_ptr
    case("Syy")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Syy))
       s%Syy(:,:) = x_2d_double_ptr
    case("Sxx")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Sxx))
       s%Sxx(:,:) = x_2d_double_ptr
    case("n")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%n))
       s%n(:,:) = x_2d_double_ptr
    case("H")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%H))
       s%H(:,:) = x_2d_double_ptr
    case("cgx")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%cgx))
       s%cgx(:,:,:) = x_3d_double_ptr
    case("cgy")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%cgy))
       s%cgy(:,:,:) = x_3d_double_ptr
    case("cx")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%cx))
       s%cx(:,:,:) = x_3d_double_ptr
    case("cy")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%cy))
       s%cy(:,:,:) = x_3d_double_ptr
    case("ctheta")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%ctheta))
       s%ctheta(:,:,:) = x_3d_double_ptr
    case("ee")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%ee))
       s%ee(:,:,:) = x_3d_double_ptr
    case("thet")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%thet))
       s%thet(:,:,:) = x_3d_double_ptr
    case("costh")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%costh))
       s%costh(:,:,:) = x_3d_double_ptr
    case("sinth")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%sinth))
       s%sinth(:,:,:) = x_3d_double_ptr
    case("sigt")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%sigt))
       s%sigt(:,:,:) = x_3d_double_ptr
    case("rr")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%rr))
       s%rr(:,:,:) = x_3d_double_ptr
    case("cgx_s")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%cgx_s))
       s%cgx_s(:,:,:) = x_3d_double_ptr
    case("cgy_s")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%cgy_s))
       s%cgy_s(:,:,:) = x_3d_double_ptr
    case("ctheta_s")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%ctheta_s))
       s%ctheta_s(:,:,:) = x_3d_double_ptr
    case("ee_s")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%ee_s))
       s%ee_s(:,:,:) = x_3d_double_ptr
    case("thet_s")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%thet_s))
       s%thet_s(:,:,:) = x_3d_double_ptr
    case("costh_s")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%costh_s))
       s%costh_s(:,:,:) = x_3d_double_ptr
    case("sinth_s")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%sinth_s))
       s%sinth_s(:,:,:) = x_3d_double_ptr
    case("k")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%k))
       s%k(:,:) = x_2d_double_ptr
    case("c")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%c))
       s%c(:,:) = x_2d_double_ptr
    case("cg")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%cg))
       s%cg(:,:) = x_2d_double_ptr
    case("sigm")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%sigm))
       s%sigm(:,:) = x_2d_double_ptr
    case("wm")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%wm))
       s%wm(:,:) = x_2d_double_ptr
    case("hh")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%hh))
       s%hh(:,:) = x_2d_double_ptr
    case("zs")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%zs))
       s%zs(:,:) = x_2d_double_ptr
    case("zs0")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%zs0))
       s%zs0(:,:) = x_2d_double_ptr
    case("tideinpt")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%tideinpt))
       s%tideinpt(:) = x_1d_double_ptr
    case("tideinpz")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%tideinpz))
       s%tideinpz(:,:) = x_2d_double_ptr
    case("windinpt")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%windinpt))
       s%windinpt(:) = x_1d_double_ptr
    case("windvelts")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%windvelts))
       s%windvelts(:) = x_1d_double_ptr
    case("winddirts")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%winddirts))
       s%winddirts(:) = x_1d_double_ptr
    case("windxts")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%windxts))
       s%windxts(:) = x_1d_double_ptr
    case("windyts")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%windyts))
       s%windyts(:) = x_1d_double_ptr
    case("windsu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%windsu))
       s%windsu(:,:) = x_2d_double_ptr
    case("windnv")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%windnv))
       s%windnv(:,:) = x_2d_double_ptr
    case("dzsdt")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dzsdt))
       s%dzsdt(:,:) = x_2d_double_ptr
    case("dzsdx")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dzsdx))
       s%dzsdx(:,:) = x_2d_double_ptr
    case("dzsdy")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dzsdy))
       s%dzsdy(:,:) = x_2d_double_ptr
    case("dzbdx")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dzbdx))
       s%dzbdx(:,:) = x_2d_double_ptr
    case("dzbdy")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dzbdy))
       s%dzbdy(:,:) = x_2d_double_ptr
    case("dzbdt")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dzbdt))
       s%dzbdt(:,:) = x_2d_double_ptr
    case("dzbnow")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dzbnow))
       s%dzbnow(:,:) = x_2d_double_ptr
    case("uu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%uu))
       s%uu(:,:) = x_2d_double_ptr
    case("vu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vu))
       s%vu(:,:) = x_2d_double_ptr
    case("uv")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%uv))
       s%uv(:,:) = x_2d_double_ptr
    case("vv")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vv))
       s%vv(:,:) = x_2d_double_ptr
    case("qx")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%qx))
       s%qx(:,:) = x_2d_double_ptr
    case("qy")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%qy))
       s%qy(:,:) = x_2d_double_ptr
    case("sedero")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%sedero))
       s%sedero(:,:) = x_2d_double_ptr
    case("dcbdx")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dcbdx))
       s%dcbdx(:,:) = x_2d_double_ptr
    case("dcbdy")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dcbdy))
       s%dcbdy(:,:) = x_2d_double_ptr
    case("dcsdx")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dcsdx))
       s%dcsdx(:,:) = x_2d_double_ptr
    case("dcsdy")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dcsdy))
       s%dcsdy(:,:) = x_2d_double_ptr
    case("ui")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%ui))
       s%ui(:,:) = x_2d_double_ptr
    case("vi")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vi))
       s%vi(:,:) = x_2d_double_ptr
    case("E")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%E))
       s%E(:,:) = x_2d_double_ptr
    case("R")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%R))
       s%R(:,:) = x_2d_double_ptr
    case("urms")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%urms))
       s%urms(:,:) = x_2d_double_ptr
    case("D")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%D))
       s%D(:,:) = x_2d_double_ptr
    case("Qb")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Qb))
       s%Qb(:,:) = x_2d_double_ptr
    case("ust")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%ust))
       s%ust(:,:) = x_2d_double_ptr
    case("ueu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%ueu))
       s%ueu(:,:) = x_2d_double_ptr
    case("vev")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vev))
       s%vev(:,:) = x_2d_double_ptr
    case("vmagu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vmagu))
       s%vmagu(:,:) = x_2d_double_ptr
    case("vmageu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vmageu))
       s%vmageu(:,:) = x_2d_double_ptr
    case("vmagv")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vmagv))
       s%vmagv(:,:) = x_2d_double_ptr
    case("vmagev")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vmagev))
       s%vmagev(:,:) = x_2d_double_ptr
    case("u")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%u))
       s%u(:,:) = x_2d_double_ptr
    case("v")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%v))
       s%v(:,:) = x_2d_double_ptr
    case("ue")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%ue))
       s%ue(:,:) = x_2d_double_ptr
    case("ve")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%ve))
       s%ve(:,:) = x_2d_double_ptr
    case("hold")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%hold))
       s%hold(:,:) = x_2d_double_ptr
    case("wetu")
       call c_f_pointer(xptr, x_2d_int_ptr, shape(s%wetu))
       s%wetu(:,:) = x_2d_int_ptr
    case("wetv")
       call c_f_pointer(xptr, x_2d_int_ptr, shape(s%wetv))
       s%wetv(:,:) = x_2d_int_ptr
    case("wetz")
       call c_f_pointer(xptr, x_2d_int_ptr, shape(s%wetz))
       s%wetz(:,:) = x_2d_int_ptr
    case("wete")
       call c_f_pointer(xptr, x_2d_int_ptr, shape(s%wete))
       s%wete(:,:) = x_2d_int_ptr
    case("hu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%hu))
       s%hu(:,:) = x_2d_double_ptr
    case("hv")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%hv))
       s%hv(:,:) = x_2d_double_ptr
    case("hum")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%hum))
       s%hum(:,:) = x_2d_double_ptr
    case("hvm")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%hvm))
       s%hvm(:,:) = x_2d_double_ptr
    case("vmag")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vmag))
       s%vmag(:,:) = x_2d_double_ptr
    case("ccg")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%ccg))
       s%ccg(:,:,:) = x_3d_double_ptr
    case("uwf")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%uwf))
       s%uwf(:,:) = x_2d_double_ptr
    case("vwf")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vwf))
       s%vwf(:,:) = x_2d_double_ptr
    case("ustr")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%ustr))
       s%ustr(:,:) = x_2d_double_ptr
    case("usd")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%usd))
       s%usd(:,:) = x_2d_double_ptr
    case("bi")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%bi))
       s%bi(:) = x_1d_double_ptr
    case("DR")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%DR))
       s%DR(:,:) = x_2d_double_ptr
    case("umean")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%umean))
       s%umean(:,:) = x_2d_double_ptr
    case("vmean")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vmean))
       s%vmean(:,:) = x_2d_double_ptr
    case("ur")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%ur))
       s%ur(:,:) = x_2d_double_ptr
    case("vardx")
       call c_f_pointer(xptr, x_0d_int_ptr) ! skip shape in case of rank 0
       s%vardx = x_0d_int_ptr
    case("D15")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%D15))
       s%D15(:) = x_1d_double_ptr
    case("D50")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%D50))
       s%D50(:) = x_1d_double_ptr
    case("D90")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%D90))
       s%D90(:) = x_1d_double_ptr
    case("sedcal")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%sedcal))
       s%sedcal(:) = x_1d_double_ptr
    case("ucrcal")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%ucrcal))
       s%ucrcal(:) = x_1d_double_ptr
    case("Tsg")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%Tsg))
       s%Tsg(:,:,:) = x_3d_double_ptr
    case("Susg")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%Susg))
       s%Susg(:,:,:) = x_3d_double_ptr
    case("Svsg")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%Svsg))
       s%Svsg(:,:,:) = x_3d_double_ptr
    case("Subg")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%Subg))
       s%Subg(:,:,:) = x_3d_double_ptr
    case("Svbg")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%Svbg))
       s%Svbg(:,:,:) = x_3d_double_ptr
    case("ceqbg")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%ceqbg))
       s%ceqbg(:,:,:) = x_3d_double_ptr
    case("ceqsg")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%ceqsg))
       s%ceqsg(:,:,:) = x_3d_double_ptr
    case("ua")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%ua))
       s%ua(:,:) = x_2d_double_ptr
    case("BR")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%BR))
       s%BR(:,:) = x_2d_double_ptr
    case("kb")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%kb))
       s%kb(:,:) = x_2d_double_ptr
    case("Tbore")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Tbore))
       s%Tbore(:,:) = x_2d_double_ptr
    case("dzav")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dzav))
       s%dzav(:,:) = x_2d_double_ptr
    case("maxzs")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%maxzs))
       s%maxzs(:,:) = x_2d_double_ptr
    case("minzs")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%minzs))
       s%minzs(:,:) = x_2d_double_ptr
    case("L1")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%L1))
       s%L1(:,:) = x_2d_double_ptr
    case("Sk")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Sk))
       s%Sk(:,:) = x_2d_double_ptr
    case("As")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%As))
       s%As(:,:) = x_2d_double_ptr
    case("gwhead")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%gwhead))
       s%gwhead(:,:) = x_2d_double_ptr
    case("gwheadb")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%gwheadb))
       s%gwheadb(:,:) = x_2d_double_ptr
    case("gwlevel")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%gwlevel))
       s%gwlevel(:,:) = x_2d_double_ptr
    case("gwheight")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%gwheight))
       s%gwheight(:,:) = x_2d_double_ptr
    case("gwbottom")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%gwbottom))
       s%gwbottom(:,:) = x_2d_double_ptr
    case("gwu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%gwu))
       s%gwu(:,:) = x_2d_double_ptr
    case("gwv")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%gwv))
       s%gwv(:,:) = x_2d_double_ptr
    case("gwqx")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%gwqx))
       s%gwqx(:,:) = x_2d_double_ptr
    case("gwqy")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%gwqy))
       s%gwqy(:,:) = x_2d_double_ptr
    case("gww")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%gww))
       s%gww(:,:) = x_2d_double_ptr
    case("gwcurv")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%gwcurv))
       s%gwcurv(:,:) = x_2d_double_ptr
    case("dinfil")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dinfil))
       s%dinfil(:,:) = x_2d_double_ptr
    case("infil")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%infil))
       s%infil(:,:) = x_2d_double_ptr
    case("gw0back")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%gw0back))
       s%gw0back(:,:) = x_2d_double_ptr
    case("kturb")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%kturb))
       s%kturb(:,:) = x_2d_double_ptr
    case("ero")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%ero))
       s%ero(:,:,:) = x_3d_double_ptr
    case("depo_im")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%depo_im))
       s%depo_im(:,:,:) = x_3d_double_ptr
    case("depo_ex")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%depo_ex))
       s%depo_ex(:,:,:) = x_3d_double_ptr
    case("nd")
       call c_f_pointer(xptr, x_2d_int_ptr, shape(s%nd))
       s%nd(:,:) = x_2d_int_ptr
    case("pbbed")
       call c_f_pointer(xptr, x_4d_double_ptr, shape(s%pbbed))
       s%pbbed(:,:,:,:) = x_4d_double_ptr
    case("dzbed")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%dzbed))
       s%dzbed(:,:,:) = x_3d_double_ptr
    case("z0bed")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%z0bed))
       s%z0bed(:,:) = x_2d_double_ptr
    case("ureps")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%ureps))
       s%ureps(:,:) = x_2d_double_ptr
    case("vreps")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vreps))
       s%vreps(:,:) = x_2d_double_ptr
    case("urepb")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%urepb))
       s%urepb(:,:) = x_2d_double_ptr
    case("vrepb")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vrepb))
       s%vrepb(:,:) = x_2d_double_ptr
    case("umwci")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%umwci))
       s%umwci(:,:) = x_2d_double_ptr
    case("vmwci")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vmwci))
       s%vmwci(:,:) = x_2d_double_ptr
    case("rolthick")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%rolthick))
       s%rolthick(:,:) = x_2d_double_ptr
    case("zswci")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%zswci))
       s%zswci(:,:) = x_2d_double_ptr
    case("pres")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%pres))
       s%pres(:,:) = x_2d_double_ptr
    case("dU")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dU))
       s%dU(:,:) = x_2d_double_ptr
    case("dV")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dV))
       s%dV(:,:) = x_2d_double_ptr
    case("wb")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%wb))
       s%wb(:,:) = x_2d_double_ptr
    case("ws")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%ws))
       s%ws(:,:) = x_2d_double_ptr
    case("wscrit")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%wscrit))
       s%wscrit(:,:) = x_2d_double_ptr
    case("bedfriccoef")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%bedfriccoef))
       s%bedfriccoef(:,:) = x_2d_double_ptr
    case("taubx")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%taubx))
       s%taubx(:,:) = x_2d_double_ptr
    case("tauby")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%tauby))
       s%tauby(:,:) = x_2d_double_ptr
    case("Df")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Df))
       s%Df(:,:) = x_2d_double_ptr
    case("Dp")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Dp))
       s%Dp(:,:) = x_2d_double_ptr
    case("Sutot")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Sutot))
       s%Sutot(:,:) = x_2d_double_ptr
    case("Svtot")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Svtot))
       s%Svtot(:,:) = x_2d_double_ptr
    case("cctot")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%cctot))
       s%cctot(:,:) = x_2d_double_ptr
    case("wi")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%wi))
       s%wi(:,:) = x_2d_double_ptr
    case("dUi")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dUi))
       s%dUi(:,:) = x_2d_double_ptr
    case("zi")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%zi))
       s%zi(:,:) = x_2d_double_ptr
    case("nuh")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%nuh))
       s%nuh(:,:) = x_2d_double_ptr
    case("cf")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%cf))
       s%cf(:,:) = x_2d_double_ptr
    case("cfu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%cfu))
       s%cfu(:,:) = x_2d_double_ptr
    case("cfv")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%cfv))
       s%cfv(:,:) = x_2d_double_ptr
    case("D50top")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%D50top))
       s%D50top(:,:) = x_2d_double_ptr
    case("D90top")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%D90top))
       s%D90top(:,:) = x_2d_double_ptr
    case("structdepth")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%structdepth))
       s%structdepth(:,:) = x_2d_double_ptr
    case("zs0fac")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%zs0fac))
       s%zs0fac(:,:,:) = x_3d_double_ptr
    case("tdisch")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%tdisch))
       s%tdisch(:) = x_1d_double_ptr
    case("pdisch")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%pdisch))
       s%pdisch(:,:) = x_2d_double_ptr
    case("pntdisch")
       call c_f_pointer(xptr, x_1d_int_ptr, shape(s%pntdisch))
       s%pntdisch(:) = x_1d_int_ptr
    case("qdisch")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%qdisch))
       s%qdisch(:,:) = x_2d_double_ptr
    case("idrift")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%idrift))
       s%idrift(:) = x_1d_double_ptr
    case("jdrift")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%jdrift))
       s%jdrift(:) = x_1d_double_ptr
    case("tdriftb")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%tdriftb))
       s%tdriftb(:) = x_1d_double_ptr
    case("tdrifte")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%tdrifte))
       s%tdrifte(:) = x_1d_double_ptr
    case("runup")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%runup))
       s%runup(:) = x_1d_double_ptr
    case("Hrunup")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%Hrunup))
       s%Hrunup(:) = x_1d_double_ptr
    case("xHrunup")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%xHrunup))
       s%xHrunup(:) = x_1d_double_ptr
    case("istruct")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%istruct))
       s%istruct(:) = x_1d_double_ptr
    case("iwl")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%iwl))
       s%iwl(:) = x_1d_double_ptr
    case("strucslope")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%strucslope))
       s%strucslope(:) = x_1d_double_ptr
    case("Dc")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Dc))
       s%Dc(:,:) = x_2d_double_ptr
    case("ph")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%ph))
       s%ph(:,:) = x_2d_double_ptr
    case("newstatbc")
       call c_f_pointer(xptr, x_0d_int_ptr) ! skip shape in case of rank 0
       s%newstatbc = x_0d_int_ptr
    case("dobs")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dobs))
       s%dobs(:,:) = x_2d_double_ptr
    case("sig2prior")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%sig2prior))
       s%sig2prior(:,:) = x_2d_double_ptr
    case("zbobs")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%zbobs))
       s%zbobs(:,:) = x_2d_double_ptr
    case("shobs")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%shobs))
       s%shobs(:,:) = x_2d_double_ptr
    case("bwalpha")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%bwalpha))
       s%bwalpha(:,:) = x_2d_double_ptr
    case("dcmdo")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dcmdo))
       s%dcmdo(:,:) = x_2d_double_ptr
    case("dassim")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%dassim))
       s%dassim(:,:) = x_2d_double_ptr
    case("cobs")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%cobs))
       s%cobs(:,:) = x_2d_double_ptr
    case("shipxCG")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%shipxCG))
       s%shipxCG(:) = x_1d_double_ptr
    case("shipyCG")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%shipyCG))
       s%shipyCG(:) = x_1d_double_ptr
    case("shipzCG")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%shipzCG))
       s%shipzCG(:) = x_1d_double_ptr
    case("shipFx")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%shipFx))
       s%shipFx(:) = x_1d_double_ptr
    case("shipFy")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%shipFy))
       s%shipFy(:) = x_1d_double_ptr
    case("shipFz")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%shipFz))
       s%shipFz(:) = x_1d_double_ptr
    case("shipMx")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%shipMx))
       s%shipMx(:) = x_1d_double_ptr
    case("shipMy")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%shipMy))
       s%shipMy(:) = x_1d_double_ptr
    case("shipMz")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%shipMz))
       s%shipMz(:) = x_1d_double_ptr
    case("shipphi")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%shipphi))
       s%shipphi(:) = x_1d_double_ptr
    case("shipchi")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%shipchi))
       s%shipchi(:) = x_1d_double_ptr
    case("shippsi")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%shippsi))
       s%shippsi(:) = x_1d_double_ptr
    case("Dveg")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Dveg))
       s%Dveg(:,:) = x_2d_double_ptr
    case("Fvegu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Fvegu))
       s%Fvegu(:,:) = x_2d_double_ptr
    case("Fvegv")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%Fvegv))
       s%Fvegv(:,:) = x_2d_double_ptr
    case("ududx")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%ududx))
       s%ududx(:,:) = x_2d_double_ptr
    case("vdvdy")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vdvdy))
       s%vdvdy(:,:) = x_2d_double_ptr
    case("udvdx")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%udvdx))
       s%udvdx(:,:) = x_2d_double_ptr
    case("vdudy")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%vdudy))
       s%vdudy(:,:) = x_2d_double_ptr
    case("viscu")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%viscu))
       s%viscu(:,:) = x_2d_double_ptr
    case("viscv")
       call c_f_pointer(xptr, x_2d_double_ptr, shape(s%viscv))
       s%viscv(:,:) = x_2d_double_ptr
    case("setbathy")
       call c_f_pointer(xptr, x_3d_double_ptr, shape(s%setbathy))
       s%setbathy(:,:,:) = x_3d_double_ptr
    case("tsetbathy")
       call c_f_pointer(xptr, x_1d_double_ptr, shape(s%tsetbathy))
       s%tsetbathy(:) = x_1d_double_ptr
    case("breaking")
       call c_f_pointer(xptr, x_2d_int_ptr, shape(s%breaking))
       s%breaking(:,:) = x_2d_int_ptr
    case default
       write(*,*) 'Setting unknown variable', var_name
    end select

  end subroutine set_var
