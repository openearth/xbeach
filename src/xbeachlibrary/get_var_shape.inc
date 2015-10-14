subroutine get_var_shape(c_var_name, shape) bind(C, name="get_var_shape")
  !DEC$ ATTRIBUTES DLLEXPORT :: get_var_shape

  character(kind=c_char), intent(in) :: c_var_name(*)
  integer(c_int), intent(inout) :: shape(MAXDIMS)

  character(len=strlen(c_var_name)) :: var_name

  var_name = char_array_to_string(c_var_name)
  shape = (/0, 0, 0, 0, 0, 0/)

  select case(var_name)

    case("x")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("y")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("xz")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("yz")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("xu")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("yu")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("xv")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("yv")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("dsu")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("dsv")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("dsz")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("dsc")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("dnu")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("dnv")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("dnz")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("dnc")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("dsdnui")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("dsdnvi")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("dsdnzi")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("alfaz")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("alfau")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("alfav")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("sdist")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("ndist")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("dx")
       shape(1) = 0
    case("dy")
       shape(1) = 0
    case("xori")
       shape(1) = 0
    case("yori")
       shape(1) = 0
    case("alfa")
       shape(1) = 0
    case("posdwn")
       shape(1) = 0
    case("thetamin")
       shape(1) = 0
    case("thetamax")
       shape(1) = 0
    case("nx")
       shape(1) = 0
    case("ny")
       shape(1) = 0
    case("zs01")
       shape(1) = 0
    case("zs02")
       shape(1) = 0
    case("zs03")
       shape(1) = 0
    case("zs04")
       shape(1) = 0
    case("xyzs01")
       ! return in c memory order
       shape(1) = 2
    case("xyzs02")
       ! return in c memory order
       shape(1) = 2
    case("xyzs03")
       ! return in c memory order
       shape(1) = 2
    case("xyzs04")
       ! return in c memory order
       shape(1) = 2
    case("tidelen")
       shape(1) = 0
    case("windlen")
       shape(1) = 0
    case("zb")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("zb0")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("theta")
       ! return in c memory order
       shape(1) = s%ntheta
    case("theta_s")
       ! return in c memory order
       shape(1) = s%ntheta_s
    case("ntheta")
       shape(1) = 0
    case("ntheta_s")
       shape(1) = 0
    case("dtheta")
       shape(1) = 0
    case("dtheta_s")
       shape(1) = 0
    case("theta0")
       shape(1) = 0
    case("thetamean")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("Fx")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("Fy")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("Sxy")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("Syy")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("Sxx")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("n")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("H")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("cgx")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = s%ntheta
    case("cgy")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = s%ntheta
    case("cx")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = s%ntheta
    case("cy")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = s%ntheta
    case("ctheta")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = s%ntheta
    case("ee")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = s%ntheta
    case("thet")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = s%ntheta
    case("costh")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = s%ntheta
    case("sinth")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = s%ntheta
    case("sigt")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = s%ntheta
    case("rr")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = s%ntheta
    case("cgx_s")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = s%ntheta_s
    case("cgy_s")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = s%ntheta_s
    case("ctheta_s")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = s%ntheta_s
    case("ee_s")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = s%ntheta_s
    case("thet_s")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = s%ntheta_s
    case("costh_s")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = s%ntheta_s
    case("sinth_s")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = s%ntheta_s
    case("k")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("c")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("cg")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("sigm")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("wm")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("hh")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("zs")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("zs0")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("tideinpt")
       ! return in c memory order
       shape(1) = s%tidelen
    case("tideinpz")
       ! return in c memory order
       shape(2) = s%tidelen
       ! return in c memory order
       shape(1) = par%tideloc
    case("windinpt")
       ! return in c memory order
       shape(1) = s%windlen
    case("windvelts")
       ! return in c memory order
       shape(1) = s%windlen
    case("winddirts")
       ! return in c memory order
       shape(1) = s%windlen
    case("windxts")
       ! return in c memory order
       shape(1) = s%windlen
    case("windyts")
       ! return in c memory order
       shape(1) = s%windlen
    case("windsu")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("windnv")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("dzsdt")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("dzsdx")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("dzsdy")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("dzbdx")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("dzbdy")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("dzbdt")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("dzbnow")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("uu")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("vu")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("uv")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("vv")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("qx")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("qy")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("sedero")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("dcbdx")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("dcbdy")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("dcsdx")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("dcsdy")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("ui")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("vi")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("E")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("R")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("urms")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("D")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("Qb")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("ust")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("ueu")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("vev")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("vmagu")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("vmageu")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("vmagv")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("vmagev")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("u")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("v")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("ue")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("ve")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("hold")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("wetu")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("wetv")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("wetz")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("hu")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("hv")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("hum")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("hvm")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("vmag")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("ccg")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = par%ngd
    case("uwf")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("vwf")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("ustr")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("usd")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("bi")
       ! return in c memory order
       shape(1) = s%ny+1
    case("DR")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("umean")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("vmean")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("ur")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("vardx")
       shape(1) = 0
    case("D15")
       ! return in c memory order
       shape(1) = par%ngd
    case("D50")
       ! return in c memory order
       shape(1) = par%ngd
    case("D90")
       ! return in c memory order
       shape(1) = par%ngd
    case("sedcal")
       ! return in c memory order
       shape(1) = par%ngd
    case("ucrcal")
       ! return in c memory order
       shape(1) = par%ngd
    case("Tsg")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = par%ngd
    case("Susg")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = par%ngd
    case("Svsg")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = par%ngd
    case("Subg")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = par%ngd
    case("Svbg")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = par%ngd
    case("ceqbg")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = par%ngd
    case("ceqsg")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = par%ngd
    case("ua")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("BR")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("kb")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("Tbore")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("dzav")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("maxzs")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("minzs")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("L1")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("Sk")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("As")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("gwhead")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("gwheadb")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("gwlevel")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("gwheight")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("gwbottom")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("gwu")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("gwv")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("gwqx")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("gwqy")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("gww")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("gwcurv")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("dinfil")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("infil")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("gw0back")
       ! return in c memory order
       shape(2) = 2
       ! return in c memory order
       shape(1) = s%ny+1
    case("kturb")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("ero")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = par%ngd
    case("depo_im")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = par%ngd
    case("depo_ex")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = par%ngd
    case("nd")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("pbbed")
       ! return in c memory order
       shape(4) = s%nx+1
       ! return in c memory order
       shape(3) = s%ny+1
       ! return in c memory order
       shape(2) = max(par%nd,2)
       ! return in c memory order
       shape(1) = par%ngd
    case("dzbed")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = max(par%nd,2)
    case("z0bed")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("ureps")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("vreps")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("urepb")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("vrepb")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("umwci")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("vmwci")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("rolthick")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("zswci")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("pres")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("dU")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("dV")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("wb")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("ws")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("wscrit")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("bedfriccoef")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("taubx")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("tauby")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("Df")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("Dp")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("Sutot")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("Svtot")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("cctot")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("wi")
       ! return in c memory order
       shape(2) = 2
       ! return in c memory order
       shape(1) = s%ny+1
    case("dUi")
       ! return in c memory order
       shape(2) = 2
       ! return in c memory order
       shape(1) = s%ny+1
    case("zi")
       ! return in c memory order
       shape(2) = 2
       ! return in c memory order
       shape(1) = s%ny+1
    case("nuh")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("cf")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("cfu")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("cfv")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("D50top")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("D90top")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("structdepth")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("zs0fac")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = 2
    case("tdisch")
       ! return in c memory order
       shape(1) = par%ntdischarge
    case("pdisch")
       ! return in c memory order
       shape(2) = par%ndischarge
       ! return in c memory order
       shape(1) = 4
    case("pntdisch")
       ! return in c memory order
       shape(1) = par%ndischarge
    case("qdisch")
       ! return in c memory order
       shape(2) = par%ntdischarge
       ! return in c memory order
       shape(1) = par%ndischarge
    case("idrift")
       ! return in c memory order
       shape(1) = par%ndrifter
    case("jdrift")
       ! return in c memory order
       shape(1) = par%ndrifter
    case("tdriftb")
       ! return in c memory order
       shape(1) = par%ndrifter
    case("tdrifte")
       ! return in c memory order
       shape(1) = par%ndrifter
    case("runup")
       ! return in c memory order
       shape(1) = s%ny+1
    case("Hrunup")
       ! return in c memory order
       shape(1) = s%ny+1
    case("xHrunup")
       ! return in c memory order
       shape(1) = s%ny+1
    case("istruct")
       ! return in c memory order
       shape(1) = s%ny+1
    case("iwl")
       ! return in c memory order
       shape(1) = s%ny+1
    case("strucslope")
       ! return in c memory order
       shape(1) = s%ny+1
    case("Dc")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("ph")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("newstatbc")
       shape(1) = 0
    case("dobs")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("sig2prior")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("zbobs")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("shobs")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("bwalpha")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("dcmdo")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("dassim")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("cobs")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("shipxCG")
       ! return in c memory order
       shape(1) = par%nship
    case("shipyCG")
       ! return in c memory order
       shape(1) = par%nship
    case("shipzCG")
       ! return in c memory order
       shape(1) = par%nship
    case("shipFx")
       ! return in c memory order
       shape(1) = par%nship
    case("shipFy")
       ! return in c memory order
       shape(1) = par%nship
    case("shipFz")
       ! return in c memory order
       shape(1) = par%nship
    case("shipMx")
       ! return in c memory order
       shape(1) = par%nship
    case("shipMy")
       ! return in c memory order
       shape(1) = par%nship
    case("shipMz")
       ! return in c memory order
       shape(1) = par%nship
    case("shipphi")
       ! return in c memory order
       shape(1) = par%nship
    case("shipchi")
       ! return in c memory order
       shape(1) = par%nship
    case("shippsi")
       ! return in c memory order
       shape(1) = par%nship
    case("vegtype")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("Cdrag")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("Dveg")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("Fvegu")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("Fvegv")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("ududx")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("vdvdy")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("udvdx")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("vdudy")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("viscu")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("viscv")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
    case("setbathy")
       ! return in c memory order
       shape(3) = s%nx+1
       ! return in c memory order
       shape(2) = s%ny+1
       ! return in c memory order
       shape(1) = par%nsetbathy
    case("tsetbathy")
       ! return in c memory order
       shape(1) = par%nsetbathy
    case("breaking")
       ! return in c memory order
       shape(2) = s%nx+1
       ! return in c memory order
       shape(1) = s%ny+1
  end select

end subroutine get_var_shape


