
subroutine get_var(c_var_name, x) bind(C, name="get_var")
  !DEC$ ATTRIBUTES DLLEXPORT :: get_var

  ! Return a pointer to the variable

  character(kind=c_char), intent(in) :: c_var_name(*)
  type(c_ptr), intent(inout) :: x

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
   case ('x')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%x(lbound(s%x, 1),lbound(s%x, 2)))
   case ('y')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%y(lbound(s%y, 1),lbound(s%y, 2)))
   case ('xz')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%xz(lbound(s%xz, 1),lbound(s%xz, 2)))
   case ('yz')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%yz(lbound(s%yz, 1),lbound(s%yz, 2)))
   case ('xu')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%xu(lbound(s%xu, 1),lbound(s%xu, 2)))
   case ('yu')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%yu(lbound(s%yu, 1),lbound(s%yu, 2)))
   case ('xv')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%xv(lbound(s%xv, 1),lbound(s%xv, 2)))
   case ('yv')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%yv(lbound(s%yv, 1),lbound(s%yv, 2)))
   case ('dsu')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dsu(lbound(s%dsu, 1),lbound(s%dsu, 2)))
   case ('dsv')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dsv(lbound(s%dsv, 1),lbound(s%dsv, 2)))
   case ('dsz')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dsz(lbound(s%dsz, 1),lbound(s%dsz, 2)))
   case ('dsc')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dsc(lbound(s%dsc, 1),lbound(s%dsc, 2)))
   case ('dnu')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dnu(lbound(s%dnu, 1),lbound(s%dnu, 2)))
   case ('dnv')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dnv(lbound(s%dnv, 1),lbound(s%dnv, 2)))
   case ('dnz')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dnz(lbound(s%dnz, 1),lbound(s%dnz, 2)))
   case ('dnc')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dnc(lbound(s%dnc, 1),lbound(s%dnc, 2)))
   case ('dsdnui')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dsdnui(lbound(s%dsdnui, 1),lbound(s%dsdnui, 2)))
   case ('dsdnvi')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dsdnvi(lbound(s%dsdnvi, 1),lbound(s%dsdnvi, 2)))
   case ('dsdnzi')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dsdnzi(lbound(s%dsdnzi, 1),lbound(s%dsdnzi, 2)))
   case ('alfaz')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%alfaz(lbound(s%alfaz, 1),lbound(s%alfaz, 2)))
   case ('alfau')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%alfau(lbound(s%alfau, 1),lbound(s%alfau, 2)))
   case ('alfav')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%alfav(lbound(s%alfav, 1),lbound(s%alfav, 2)))
   case ('sdist')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%sdist(lbound(s%sdist, 1),lbound(s%sdist, 2)))
   case ('ndist')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%ndist(lbound(s%ndist, 1),lbound(s%ndist, 2)))
   case ('dx')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dx)
   case ('dy')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dy)
   case ('xori')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%xori)
   case ('yori')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%yori)
   case ('alfa')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%alfa)
   case ('posdwn')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%posdwn)
   case ('thetamin')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%thetamin)
   case ('thetamax')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%thetamax)
   case ('nx')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%nx)
   case ('ny')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%ny)
   case ('zs01')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%zs01)
   case ('zs02')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%zs02)
   case ('zs03')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%zs03)
   case ('zs04')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%zs04)
   case ('xyzs01')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%xyzs01(lbound(s%xyzs01, 1)))
   case ('xyzs02')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%xyzs02(lbound(s%xyzs02, 1)))
   case ('xyzs03')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%xyzs03(lbound(s%xyzs03, 1)))
   case ('xyzs04')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%xyzs04(lbound(s%xyzs04, 1)))
   case ('tidelen')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%tidelen)
   case ('windlen')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%windlen)
   case ('zb')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%zb(lbound(s%zb, 1),lbound(s%zb, 2)))
   case ('zb0')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%zb0(lbound(s%zb0, 1),lbound(s%zb0, 2)))
   case ('theta')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%theta(lbound(s%theta, 1)))
   case ('theta_s')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%theta_s(lbound(s%theta_s, 1)))
   case ('ntheta')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%ntheta)
   case ('ntheta_s')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%ntheta_s)
   case ('dtheta')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dtheta)
   case ('dtheta_s')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dtheta_s)
   case ('theta0')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%theta0)
   case ('thetamean')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%thetamean(lbound(s%thetamean, 1),lbound(s%thetamean, 2)))
   case ('Fx')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%Fx(lbound(s%Fx, 1),lbound(s%Fx, 2)))
   case ('Fy')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%Fy(lbound(s%Fy, 1),lbound(s%Fy, 2)))
   case ('Sxy')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%Sxy(lbound(s%Sxy, 1),lbound(s%Sxy, 2)))
   case ('Syy')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%Syy(lbound(s%Syy, 1),lbound(s%Syy, 2)))
   case ('Sxx')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%Sxx(lbound(s%Sxx, 1),lbound(s%Sxx, 2)))
   case ('n')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%n(lbound(s%n, 1),lbound(s%n, 2)))
   case ('H')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%H(lbound(s%H, 1),lbound(s%H, 2)))
   case ('cgx')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%cgx(lbound(s%cgx, 1),lbound(s%cgx, 2),lbound(s%cgx, 3)))
   case ('cgy')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%cgy(lbound(s%cgy, 1),lbound(s%cgy, 2),lbound(s%cgy, 3)))
   case ('cx')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%cx(lbound(s%cx, 1),lbound(s%cx, 2),lbound(s%cx, 3)))
   case ('cy')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%cy(lbound(s%cy, 1),lbound(s%cy, 2),lbound(s%cy, 3)))
   case ('ctheta')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%ctheta(lbound(s%ctheta, 1),lbound(s%ctheta, 2),lbound(s%ctheta, 3)))
   case ('ee')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%ee(lbound(s%ee, 1),lbound(s%ee, 2),lbound(s%ee, 3)))
   case ('thet')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%thet(lbound(s%thet, 1),lbound(s%thet, 2),lbound(s%thet, 3)))
   case ('costh')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%costh(lbound(s%costh, 1),lbound(s%costh, 2),lbound(s%costh, 3)))
   case ('sinth')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%sinth(lbound(s%sinth, 1),lbound(s%sinth, 2),lbound(s%sinth, 3)))
   case ('sigt')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%sigt(lbound(s%sigt, 1),lbound(s%sigt, 2),lbound(s%sigt, 3)))
   case ('rr')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%rr(lbound(s%rr, 1),lbound(s%rr, 2),lbound(s%rr, 3)))
   case ('cgx_s')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%cgx_s(lbound(s%cgx_s, 1),lbound(s%cgx_s, 2),lbound(s%cgx_s, 3)))
   case ('cgy_s')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%cgy_s(lbound(s%cgy_s, 1),lbound(s%cgy_s, 2),lbound(s%cgy_s, 3)))
   case ('ctheta_s')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%ctheta_s(lbound(s%ctheta_s, 1),lbound(s%ctheta_s, 2),lbound(s%ctheta_s, 3)))
   case ('ee_s')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%ee_s(lbound(s%ee_s, 1),lbound(s%ee_s, 2),lbound(s%ee_s, 3)))
   case ('thet_s')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%thet_s(lbound(s%thet_s, 1),lbound(s%thet_s, 2),lbound(s%thet_s, 3)))
   case ('costh_s')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%costh_s(lbound(s%costh_s, 1),lbound(s%costh_s, 2),lbound(s%costh_s, 3)))
   case ('sinth_s')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%sinth_s(lbound(s%sinth_s, 1),lbound(s%sinth_s, 2),lbound(s%sinth_s, 3)))
   case ('k')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%k(lbound(s%k, 1),lbound(s%k, 2)))
   case ('c')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%c(lbound(s%c, 1),lbound(s%c, 2)))
   case ('cg')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%cg(lbound(s%cg, 1),lbound(s%cg, 2)))
   case ('sigm')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%sigm(lbound(s%sigm, 1),lbound(s%sigm, 2)))
   case ('wm')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%wm(lbound(s%wm, 1),lbound(s%wm, 2)))
   case ('hh')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%hh(lbound(s%hh, 1),lbound(s%hh, 2)))
   case ('zs')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%zs(lbound(s%zs, 1),lbound(s%zs, 2)))
   case ('zs0')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%zs0(lbound(s%zs0, 1),lbound(s%zs0, 2)))
   case ('tideinpt')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%tideinpt(lbound(s%tideinpt, 1)))
   case ('tideinpz')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%tideinpz(lbound(s%tideinpz, 1),lbound(s%tideinpz, 2)))
   case ('windinpt')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%windinpt(lbound(s%windinpt, 1)))
   case ('windvelts')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%windvelts(lbound(s%windvelts, 1)))
   case ('winddirts')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%winddirts(lbound(s%winddirts, 1)))
   case ('windxts')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%windxts(lbound(s%windxts, 1)))
   case ('windyts')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%windyts(lbound(s%windyts, 1)))
   case ('windsu')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%windsu(lbound(s%windsu, 1),lbound(s%windsu, 2)))
   case ('windnv')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%windnv(lbound(s%windnv, 1),lbound(s%windnv, 2)))
   case ('dzsdt')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dzsdt(lbound(s%dzsdt, 1),lbound(s%dzsdt, 2)))
   case ('dzsdx')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dzsdx(lbound(s%dzsdx, 1),lbound(s%dzsdx, 2)))
   case ('dzsdy')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dzsdy(lbound(s%dzsdy, 1),lbound(s%dzsdy, 2)))
   case ('dzbdx')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dzbdx(lbound(s%dzbdx, 1),lbound(s%dzbdx, 2)))
   case ('dzbdy')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dzbdy(lbound(s%dzbdy, 1),lbound(s%dzbdy, 2)))
   case ('dzbdt')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dzbdt(lbound(s%dzbdt, 1),lbound(s%dzbdt, 2)))
   case ('dzbnow')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dzbnow(lbound(s%dzbnow, 1),lbound(s%dzbnow, 2)))
   case ('uu')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%uu(lbound(s%uu, 1),lbound(s%uu, 2)))
   case ('vu')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%vu(lbound(s%vu, 1),lbound(s%vu, 2)))
   case ('uv')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%uv(lbound(s%uv, 1),lbound(s%uv, 2)))
   case ('vv')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%vv(lbound(s%vv, 1),lbound(s%vv, 2)))
   case ('qx')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%qx(lbound(s%qx, 1),lbound(s%qx, 2)))
   case ('qy')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%qy(lbound(s%qy, 1),lbound(s%qy, 2)))
   case ('sedero')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%sedero(lbound(s%sedero, 1),lbound(s%sedero, 2)))
   case ('dcbdx')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dcbdx(lbound(s%dcbdx, 1),lbound(s%dcbdx, 2)))
   case ('dcbdy')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dcbdy(lbound(s%dcbdy, 1),lbound(s%dcbdy, 2)))
   case ('dcsdx')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dcsdx(lbound(s%dcsdx, 1),lbound(s%dcsdx, 2)))
   case ('dcsdy')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dcsdy(lbound(s%dcsdy, 1),lbound(s%dcsdy, 2)))
   case ('ui')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%ui(lbound(s%ui, 1),lbound(s%ui, 2)))
   case ('vi')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%vi(lbound(s%vi, 1),lbound(s%vi, 2)))
   case ('E')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%E(lbound(s%E, 1),lbound(s%E, 2)))
   case ('R')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%R(lbound(s%R, 1),lbound(s%R, 2)))
   case ('urms')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%urms(lbound(s%urms, 1),lbound(s%urms, 2)))
   case ('D')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%D(lbound(s%D, 1),lbound(s%D, 2)))
   case ('Qb')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%Qb(lbound(s%Qb, 1),lbound(s%Qb, 2)))
   case ('ust')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%ust(lbound(s%ust, 1),lbound(s%ust, 2)))
   case ('ueu')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%ueu(lbound(s%ueu, 1),lbound(s%ueu, 2)))
   case ('vev')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%vev(lbound(s%vev, 1),lbound(s%vev, 2)))
   case ('vmagu')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%vmagu(lbound(s%vmagu, 1),lbound(s%vmagu, 2)))
   case ('vmageu')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%vmageu(lbound(s%vmageu, 1),lbound(s%vmageu, 2)))
   case ('vmagv')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%vmagv(lbound(s%vmagv, 1),lbound(s%vmagv, 2)))
   case ('vmagev')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%vmagev(lbound(s%vmagev, 1),lbound(s%vmagev, 2)))
   case ('u')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%u(lbound(s%u, 1),lbound(s%u, 2)))
   case ('v')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%v(lbound(s%v, 1),lbound(s%v, 2)))
   case ('ue')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%ue(lbound(s%ue, 1),lbound(s%ue, 2)))
   case ('ve')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%ve(lbound(s%ve, 1),lbound(s%ve, 2)))
   case ('hold')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%hold(lbound(s%hold, 1),lbound(s%hold, 2)))
   case ('wetu')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%wetu(lbound(s%wetu, 1),lbound(s%wetu, 2)))
   case ('wetv')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%wetv(lbound(s%wetv, 1),lbound(s%wetv, 2)))
   case ('wetz')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%wetz(lbound(s%wetz, 1),lbound(s%wetz, 2)))
   case ('wete')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%wete(lbound(s%wete, 1),lbound(s%wete, 2)))
   case ('hu')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%hu(lbound(s%hu, 1),lbound(s%hu, 2)))
   case ('hv')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%hv(lbound(s%hv, 1),lbound(s%hv, 2)))
   case ('hum')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%hum(lbound(s%hum, 1),lbound(s%hum, 2)))
   case ('hvm')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%hvm(lbound(s%hvm, 1),lbound(s%hvm, 2)))
   case ('vmag')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%vmag(lbound(s%vmag, 1),lbound(s%vmag, 2)))
   case ('ccg')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%ccg(lbound(s%ccg, 1),lbound(s%ccg, 2),lbound(s%ccg, 3)))
   case ('uwf')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%uwf(lbound(s%uwf, 1),lbound(s%uwf, 2)))
   case ('vwf')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%vwf(lbound(s%vwf, 1),lbound(s%vwf, 2)))
   case ('ustr')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%ustr(lbound(s%ustr, 1),lbound(s%ustr, 2)))
   case ('usd')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%usd(lbound(s%usd, 1),lbound(s%usd, 2)))
   case ('bi')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%bi(lbound(s%bi, 1)))
   case ('DR')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%DR(lbound(s%DR, 1),lbound(s%DR, 2)))
   case ('umean')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%umean(lbound(s%umean, 1),lbound(s%umean, 2)))
   case ('vmean')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%vmean(lbound(s%vmean, 1),lbound(s%vmean, 2)))
   case ('ur')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%ur(lbound(s%ur, 1),lbound(s%ur, 2)))
   case ('vardx')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%vardx)
   case ('D15')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%D15(lbound(s%D15, 1)))
   case ('D50')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%D50(lbound(s%D50, 1)))
   case ('D90')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%D90(lbound(s%D90, 1)))
   case ('sedcal')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%sedcal(lbound(s%sedcal, 1)))
   case ('ucrcal')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%ucrcal(lbound(s%ucrcal, 1)))
   case ('Tsg')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%Tsg(lbound(s%Tsg, 1),lbound(s%Tsg, 2),lbound(s%Tsg, 3)))
   case ('Susg')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%Susg(lbound(s%Susg, 1),lbound(s%Susg, 2),lbound(s%Susg, 3)))
   case ('Svsg')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%Svsg(lbound(s%Svsg, 1),lbound(s%Svsg, 2),lbound(s%Svsg, 3)))
   case ('Subg')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%Subg(lbound(s%Subg, 1),lbound(s%Subg, 2),lbound(s%Subg, 3)))
   case ('Svbg')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%Svbg(lbound(s%Svbg, 1),lbound(s%Svbg, 2),lbound(s%Svbg, 3)))
   case ('ceqbg')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%ceqbg(lbound(s%ceqbg, 1),lbound(s%ceqbg, 2),lbound(s%ceqbg, 3)))
   case ('ceqsg')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%ceqsg(lbound(s%ceqsg, 1),lbound(s%ceqsg, 2),lbound(s%ceqsg, 3)))
   case ('ua')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%ua(lbound(s%ua, 1),lbound(s%ua, 2)))
   case ('BR')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%BR(lbound(s%BR, 1),lbound(s%BR, 2)))
   case ('kb')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%kb(lbound(s%kb, 1),lbound(s%kb, 2)))
   case ('Tbore')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%Tbore(lbound(s%Tbore, 1),lbound(s%Tbore, 2)))
   case ('dzav')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dzav(lbound(s%dzav, 1),lbound(s%dzav, 2)))
   case ('maxzs')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%maxzs(lbound(s%maxzs, 1),lbound(s%maxzs, 2)))
   case ('minzs')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%minzs(lbound(s%minzs, 1),lbound(s%minzs, 2)))
   case ('L1')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%L1(lbound(s%L1, 1),lbound(s%L1, 2)))
   case ('Sk')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%Sk(lbound(s%Sk, 1),lbound(s%Sk, 2)))
   case ('As')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%As(lbound(s%As, 1),lbound(s%As, 2)))
   case ('gwhead')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%gwhead(lbound(s%gwhead, 1),lbound(s%gwhead, 2)))
   case ('gwheadb')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%gwheadb(lbound(s%gwheadb, 1),lbound(s%gwheadb, 2)))
   case ('gwlevel')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%gwlevel(lbound(s%gwlevel, 1),lbound(s%gwlevel, 2)))
   case ('gwheight')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%gwheight(lbound(s%gwheight, 1),lbound(s%gwheight, 2)))
   case ('gwbottom')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%gwbottom(lbound(s%gwbottom, 1),lbound(s%gwbottom, 2)))
   case ('gwu')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%gwu(lbound(s%gwu, 1),lbound(s%gwu, 2)))
   case ('gwv')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%gwv(lbound(s%gwv, 1),lbound(s%gwv, 2)))
   case ('gwqx')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%gwqx(lbound(s%gwqx, 1),lbound(s%gwqx, 2)))
   case ('gwqy')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%gwqy(lbound(s%gwqy, 1),lbound(s%gwqy, 2)))
   case ('gww')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%gww(lbound(s%gww, 1),lbound(s%gww, 2)))
   case ('gwcurv')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%gwcurv(lbound(s%gwcurv, 1),lbound(s%gwcurv, 2)))
   case ('dinfil')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dinfil(lbound(s%dinfil, 1),lbound(s%dinfil, 2)))
   case ('infil')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%infil(lbound(s%infil, 1),lbound(s%infil, 2)))
   case ('gw0back')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%gw0back(lbound(s%gw0back, 1),lbound(s%gw0back, 2)))
   case ('kturb')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%kturb(lbound(s%kturb, 1),lbound(s%kturb, 2)))
   case ('ero')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%ero(lbound(s%ero, 1),lbound(s%ero, 2),lbound(s%ero, 3)))
   case ('depo_im')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%depo_im(lbound(s%depo_im, 1),lbound(s%depo_im, 2),lbound(s%depo_im, 3)))
   case ('depo_ex')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%depo_ex(lbound(s%depo_ex, 1),lbound(s%depo_ex, 2),lbound(s%depo_ex, 3)))
   case ('nd')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%nd(lbound(s%nd, 1),lbound(s%nd, 2)))
   case ('pbbed')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%pbbed(lbound(s%pbbed, 1),lbound(s%pbbed, 2),lbound(s%pbbed, 3),lbound(s%pbbed, 4)))
   case ('dzbed')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dzbed(lbound(s%dzbed, 1),lbound(s%dzbed, 2),lbound(s%dzbed, 3)))
   case ('z0bed')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%z0bed(lbound(s%z0bed, 1),lbound(s%z0bed, 2)))
   case ('ureps')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%ureps(lbound(s%ureps, 1),lbound(s%ureps, 2)))
   case ('vreps')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%vreps(lbound(s%vreps, 1),lbound(s%vreps, 2)))
   case ('urepb')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%urepb(lbound(s%urepb, 1),lbound(s%urepb, 2)))
   case ('vrepb')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%vrepb(lbound(s%vrepb, 1),lbound(s%vrepb, 2)))
   case ('umwci')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%umwci(lbound(s%umwci, 1),lbound(s%umwci, 2)))
   case ('vmwci')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%vmwci(lbound(s%vmwci, 1),lbound(s%vmwci, 2)))
   case ('rolthick')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%rolthick(lbound(s%rolthick, 1),lbound(s%rolthick, 2)))
   case ('zswci')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%zswci(lbound(s%zswci, 1),lbound(s%zswci, 2)))
   case ('pres')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%pres(lbound(s%pres, 1),lbound(s%pres, 2)))
   case ('dU')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dU(lbound(s%dU, 1),lbound(s%dU, 2)))
   case ('dV')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dV(lbound(s%dV, 1),lbound(s%dV, 2)))
   case ('wb')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%wb(lbound(s%wb, 1),lbound(s%wb, 2)))
   case ('ws')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%ws(lbound(s%ws, 1),lbound(s%ws, 2)))
   case ('wscrit')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%wscrit(lbound(s%wscrit, 1),lbound(s%wscrit, 2)))
   case ('bedfriccoef')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%bedfriccoef(lbound(s%bedfriccoef, 1),lbound(s%bedfriccoef, 2)))
   case ('taubx')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%taubx(lbound(s%taubx, 1),lbound(s%taubx, 2)))
   case ('tauby')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%tauby(lbound(s%tauby, 1),lbound(s%tauby, 2)))
   case ('Df')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%Df(lbound(s%Df, 1),lbound(s%Df, 2)))
   case ('Dp')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%Dp(lbound(s%Dp, 1),lbound(s%Dp, 2)))
   case ('Sutot')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%Sutot(lbound(s%Sutot, 1),lbound(s%Sutot, 2)))
   case ('Svtot')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%Svtot(lbound(s%Svtot, 1),lbound(s%Svtot, 2)))
   case ('cctot')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%cctot(lbound(s%cctot, 1),lbound(s%cctot, 2)))
   case ('wi')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%wi(lbound(s%wi, 1),lbound(s%wi, 2)))
   case ('dUi')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dUi(lbound(s%dUi, 1),lbound(s%dUi, 2)))
   case ('zi')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%zi(lbound(s%zi, 1),lbound(s%zi, 2)))
   case ('nuh')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%nuh(lbound(s%nuh, 1),lbound(s%nuh, 2)))
   case ('cf')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%cf(lbound(s%cf, 1),lbound(s%cf, 2)))
   case ('cfu')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%cfu(lbound(s%cfu, 1),lbound(s%cfu, 2)))
   case ('cfv')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%cfv(lbound(s%cfv, 1),lbound(s%cfv, 2)))
   case ('D50top')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%D50top(lbound(s%D50top, 1),lbound(s%D50top, 2)))
   case ('D90top')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%D90top(lbound(s%D90top, 1),lbound(s%D90top, 2)))
   case ('structdepth')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%structdepth(lbound(s%structdepth, 1),lbound(s%structdepth, 2)))
   case ('zs0fac')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%zs0fac(lbound(s%zs0fac, 1),lbound(s%zs0fac, 2),lbound(s%zs0fac, 3)))
   case ('tdisch')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%tdisch(lbound(s%tdisch, 1)))
   case ('pdisch')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%pdisch(lbound(s%pdisch, 1),lbound(s%pdisch, 2)))
   case ('pntdisch')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%pntdisch(lbound(s%pntdisch, 1)))
   case ('qdisch')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%qdisch(lbound(s%qdisch, 1),lbound(s%qdisch, 2)))
   case ('idrift')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%idrift(lbound(s%idrift, 1)))
   case ('jdrift')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%jdrift(lbound(s%jdrift, 1)))
   case ('tdriftb')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%tdriftb(lbound(s%tdriftb, 1)))
   case ('tdrifte')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%tdrifte(lbound(s%tdrifte, 1)))
   case ('runup')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%runup(lbound(s%runup, 1)))
   case ('Hrunup')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%Hrunup(lbound(s%Hrunup, 1)))
   case ('xHrunup')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%xHrunup(lbound(s%xHrunup, 1)))
   case ('istruct')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%istruct(lbound(s%istruct, 1)))
   case ('iwl')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%iwl(lbound(s%iwl, 1)))
   case ('strucslope')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%strucslope(lbound(s%strucslope, 1)))
   case ('Dc')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%Dc(lbound(s%Dc, 1),lbound(s%Dc, 2)))
   case ('ph')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%ph(lbound(s%ph, 1),lbound(s%ph, 2)))
   case ('newstatbc')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%newstatbc)
   case ('dobs')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dobs(lbound(s%dobs, 1),lbound(s%dobs, 2)))
   case ('sig2prior')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%sig2prior(lbound(s%sig2prior, 1),lbound(s%sig2prior, 2)))
   case ('zbobs')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%zbobs(lbound(s%zbobs, 1),lbound(s%zbobs, 2)))
   case ('shobs')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%shobs(lbound(s%shobs, 1),lbound(s%shobs, 2)))
   case ('bwalpha')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%bwalpha(lbound(s%bwalpha, 1),lbound(s%bwalpha, 2)))
   case ('dcmdo')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dcmdo(lbound(s%dcmdo, 1),lbound(s%dcmdo, 2)))
   case ('dassim')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%dassim(lbound(s%dassim, 1),lbound(s%dassim, 2)))
   case ('cobs')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%cobs(lbound(s%cobs, 1),lbound(s%cobs, 2)))
   case ('shipxCG')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%shipxCG(lbound(s%shipxCG, 1)))
   case ('shipyCG')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%shipyCG(lbound(s%shipyCG, 1)))
   case ('shipzCG')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%shipzCG(lbound(s%shipzCG, 1)))
   case ('shipFx')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%shipFx(lbound(s%shipFx, 1)))
   case ('shipFy')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%shipFy(lbound(s%shipFy, 1)))
   case ('shipFz')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%shipFz(lbound(s%shipFz, 1)))
   case ('shipMx')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%shipMx(lbound(s%shipMx, 1)))
   case ('shipMy')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%shipMy(lbound(s%shipMy, 1)))
   case ('shipMz')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%shipMz(lbound(s%shipMz, 1)))
   case ('shipphi')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%shipphi(lbound(s%shipphi, 1)))
   case ('shipchi')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%shipchi(lbound(s%shipchi, 1)))
   case ('shippsi')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%shippsi(lbound(s%shippsi, 1)))
   case ('Dveg')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%Dveg(lbound(s%Dveg, 1),lbound(s%Dveg, 2)))
   case ('Fvegu')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%Fvegu(lbound(s%Fvegu, 1),lbound(s%Fvegu, 2)))
   case ('Fvegv')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%Fvegv(lbound(s%Fvegv, 1),lbound(s%Fvegv, 2)))
   case ('ududx')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%ududx(lbound(s%ududx, 1),lbound(s%ududx, 2)))
   case ('vdvdy')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%vdvdy(lbound(s%vdvdy, 1),lbound(s%vdvdy, 2)))
   case ('udvdx')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%udvdx(lbound(s%udvdx, 1),lbound(s%udvdx, 2)))
   case ('vdudy')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%vdudy(lbound(s%vdudy, 1),lbound(s%vdudy, 2)))
   case ('viscu')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%viscu(lbound(s%viscu, 1),lbound(s%viscu, 2)))
   case ('viscv')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%viscv(lbound(s%viscv, 1),lbound(s%viscv, 2)))
   case ('setbathy')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%setbathy(lbound(s%setbathy, 1),lbound(s%setbathy, 2),lbound(s%setbathy, 3)))
   case ('tsetbathy')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%tsetbathy(lbound(s%tsetbathy, 1)))
   case ('breaking')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%breaking(lbound(s%breaking, 1),lbound(s%breaking, 2)))
  end select

end subroutine get_var
