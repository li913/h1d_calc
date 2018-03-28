module FunMod_Van
implicit none
private
public ModelNum, xConv, tConv, FS, FH, FK, FQ, FC, FKcap, ReadUserDefinedPar
integer, parameter :: ModelNum = 0        ! 模型编号
real(4), parameter :: h0 = -6.3e6         ! 含水率为0时土壤吸力
real(8), parameter :: hL = h0*(1d0-1d-20) ! h取值的左边界
real(8), parameter :: hR = -1d-30         ! h取值的右边界
real(8):: xConv=100                       ! 米转为 Hydrus 单位的转换系数, 于input.for-Conversion中赋值
real(8):: tConv=1./(60.*60.*24.)          ! 秒转为 Hydrus 单位的转换系数, 于input.for-Conversion中赋值
real(8):: omega0
logical(4):: bSad
  contains
!计算θ(h)
pure function FQ( iModel, h, Par ) result(Res)
  implicit none
  integer(4), intent(in):: iModel
  real(4),    intent(in):: h, Par(:)
  real(4) Res
  real(8) Qr, Qs
  
  call SetPar( Par, Qr=Qr, Qs=Qs )
  if (h<hL) then
  	Res = 1e-37
  elseif(h<hR) then
    if(bSad) then
      Res = (Qs-Qr)*Scap( h, Par ) + Qr*Sad( h, Par )
    else
      !Res = Qs*Scap( h, Par )
      Res = (Qs-Qr)*Scap( h, Par ) + Qr
    end if
    Res = max( Res, 1e-37 )    
  else
    Res = Qs
  end if
end function FQ
!计算θ(h)的一阶导数
pure function FC( iModel, h, Par ) result(Res)
  implicit none
  integer(4), intent(in):: iModel
  real(4),    intent(in):: h, Par(:)
  real(4) Res
  real(8) Qr, Qs
  real(8) a, b
    
  if(h>=hR .OR. h<=hL) then
    Res = 0.0; return
  end if
  call SetPar( Par, Qr=Qr, Qs=Qs )
  a = Derivative_Scap( h, Par )  !Scap对h的导数  
  b = Derivative_Sad ( h, Par )  !Sad对h的导数
  if(bSad) then
    Res = (Qs-Qr)*a + Qr*b
  else
    !Res = Qs*a
    Res = (Qs-Qr)*a
  end if
end function FC

!计算Se(θ)
pure function FS( iModel, h, Par ) result(Res)
  implicit none
  integer(4), intent(in):: iModel
  real(4),    intent(in):: h, Par(:)
  real(4) Res
  real(8) Q, Qr, Qs
  
  call SetPar( Par, Qr=Qr, Qs=Qs )
  Q = FQ( iModel, h, Par )
  Res = (Q-Qr) / (Qs-Qr)
end function FS

!计算h(θ)
pure function FH( iModel, Qe, Par) result(Res)
  implicit none
  integer(4), intent(in):: iModel
  real(4),    intent(in):: Qe, Par(:)
  real(4) Res
  real(8) Q, Qr, Qs
   
  call SetPar( Par, Qr=Qr, Qs=Qs )
  Q = Qe * (Qs-Qr) + Qr
  Res = SolveFQ( Q, FQ, Par )
end function FH  

!计算K(h)
pure function FK( iModel, h, Par ) result(Res)
  implicit none
  integer(4), intent(in):: iModel
  real(4),    intent(in):: h, Par(:)
  real(4) Res
  real(8) Ks, omega
  real(8) s1, s2
  
  call SetPar( Par, Ks=Ks, omega=omega )
  if(h>hR) then
    Res = Ks; return
  end if  
  s1 = (1d0-omega)*Kcap( h, Par ) 
  s2 = omega*Kad( h, Par )
  Res = max( Ks*(s1+s2), 1d-38 ) 
end function FK
pure function FKcap( iModel, h, Par ) result(Res)
  implicit none
  integer(4), intent(in):: iModel
  real(4),    intent(in):: h, Par(:)
  real(4) Res
  real(8) Ks, omega
  real(8) s1, s2
  
  call SetPar( Par, Ks=Ks, omega=omega )
  if(h>hR) then
    Res = Ks*(1d0-omega); return
  end if
  s1 = (1d0-omega)*Kcap( h, Par ) 
  s2 = omega*Kad( h, Par )
  Res = max( Ks*s1, 1d-38 ) 
end function FKcap

!原毛细饱和度方程 Γ(h)
pure function gamma( h, Par ) result(Res)
  implicit none
  real(4), intent(in):: h, Par(:)
  real(8) Res
  real(8) alpha, rn, rm
  
  call SetPar( Par, alpha=alpha, rn=rn, rm=rm )
  Res = 1d0 + (-alpha*h)**rn
  Res = Res**(-rm)
end function gamma
! gamma函数Γ(h)对h的导数：dΓ/dh
pure function Derivative_gamma( h, Par ) result(Res)
  implicit none
  real(4), intent(in):: h, Par(:)
  real(8) Res
  real(8) alpha, rn, rm
  
  call SetPar( Par, alpha=alpha, rn=rn, rm=rm )
  Res = ( 1d0 + (-alpha*h)**rn )**(rm+1)
  Res = rm*rn*(alpha**rn) * (-h)**(rn-1d0) / Res
end function Derivative_gamma

!光滑因子b
pure function SmoothingPar( Par ) Result(Res)
  real(4), intent(in)   :: Par(:)
  real(8) Res
  real(8) Qr, Qs, rn
  
  call SetPar( Par, Qr=Qr, Qs=Qs, rn=rn )
  Res = 1d0 - exp( -(Qr/(Qs-Qr))**2 )
  Res = 0.1d0 + 0.2d0/(rn*rn)*Res
end function SmoothingPar

!新毛细饱和度方程
pure function Scap( h, Par ) result(Res)  
  implicit none
  real(4), intent(in):: h, Par(:)
  real(8) Res

  Res = gamma( h0, Par )
  Res = ( gamma( h , Par )-Res ) / ( 1d0-Res )
end function Scap
!Scap对h的导数
pure function Derivative_Scap( h, Par ) result(Res)
  implicit none
  real(4), intent(in):: h, Par(:)
  real(8) Res
  Res = Derivative_gamma( h, Par )     !gamma函数对h的导数
  Res = Res / ( 1d0-gamma( h0, Par ) )
end function Derivative_Scap

!薄膜水饱和度
pure function Sad( h, Par ) result(Res)
  implicit none
  real(4), intent(in):: h, Par(:) 
  real(8) Res
  real(8) alpha
  real(8), parameter :: x0 = log10(-h0)
  real(8) xa, x, b
  
  call SetPar( Par, alpha=alpha )
  x = log10( max(-h,1e-30) ); xa = log10( 1d0/max(alpha,1d-37) )
  b = SmoothingPar( Par ) !光滑因子b
  Res = 1d0 + exp((xa-x)/b)
  Res = x - xa + b*log(Res)
  Res = 1d0 + Res/(xa-x0)
end function Sad  
!薄膜水饱和度的导数：dSad/dh
pure function Derivative_Sad( h, Par ) result(Res)
  implicit none
  real(4), intent(in):: h, Par(:)
  real(8) Res
  real(8) alpha
  real(8), parameter :: x0 = log10(-h0)
  real(8) xa, x, b, e, dxdh
  
  call SetPar( Par, alpha=alpha )
  x = log10(-h); xa = log10( 1d0/max(alpha,1d-37) )
  b = SmoothingPar( par )    !光滑因子b
  e = exp( (xa-x)/b )        ! exp()
  dxdh = 1d0 / (h*log(10d0)) ! dx/dh, x=log10(-h)
  Res = dxdh / (xa-x0) / (1d0+e)
end function Derivative_Sad

!毛细相对导水率
pure function Kcap( h, Par ) result(Res)
  implicit none
  real(4), intent(in):: h, Par(:)
  real(8) Res
  real(8) rm, s, rL
  
  call SetPar( Par, rm=rm, rL=rL )
  s = Scap( h, Par )
  if(s<1d-37) then
    Res = 0d0
  else
    Res = 1d0 - ( 1d0 - s**(1d0/rm) )**rm
    Res = (s**rL) * Res*Res
  end if
end function Kcap

!薄膜相对导水率
pure function Kad( h, Par ) result(Res)
  implicit none
  real(4), intent(in):: h, Par(:)
  real(8) Res
  real(8) alpha, s
  
  call SetPar( Par, alpha=alpha )
  s = Sad( h, Par ); s = -1.5d0 * (1d0-s)
  Res = (-h0*alpha)**s
end function Kad

!气态导水率
elemental function Kvap( h, Q, Qs ) result(Res)
  implicit none
  real(8), intent(in):: h, Q, Qs
  real(8),  parameter:: temp = 20d0      !摄氏温度
  real(8),  parameter:: T0 = 273.15d0    !绝对零度
  real(8),  parameter:: M = 0.018015d0   !水的摩尔质量
  real(8),  parameter:: g = 9.81d0       !重力加速度
  real(8),  parameter:: R = 8.314d0      !通用气体常数
  real(8) h_m, T, Da, Qa, tau, D, row, rosv, Hr
  real(8) Res
  h_m = h / xConv !单位转换为米(m)
  T = Temp + T0   !转为开尔文温度  
  
  !计算水汽扩散系数 D
  Qa = max(Qs-Q,0d0) !空气体积含量
  tau= Qa**(7d0/3d0) / (Qs*Qs) !曲折因子
  Da = 2.12d-5 * (T/T0)**2 !T温度下水汽扩散系数
  D = tau * Qa * Da !水汽扩散系数
  !水和蒸汽的密度
  row = 1d3 * ( 1d0 - 7.37d-6*(Temp-4d0)**2 + 3.79d-8*(Temp-4d0)**3 )
  rosv= 1d-3 * exp( 31.3716d0 - 6014.79d0/T - 7.92495d-3*T ) / T
  !相对湿度
  Hr = exp( h_m*M*g/(R*T) )
  Res = D*rosv*M*g*Hr / (row*R*T)
  Res = Res*xConv/tConv !转为 Hydrus 的量纲         
end function Kvap

! 二分法求解非线性单调函数, 由θ反算h
! Q 函数值
pure function SolveFQ( Q, FQ, Par ) result(Res)
  implicit none
  real(8), intent(in):: Q
  real(4), intent(in):: Par(:)
  real(8) Res
  interface
    pure function FQ( iModel, h, Par ) result(Res)
      implicit none
      integer(4), intent(in):: iModel
      real(4),    intent(in):: h, Par(:)
      real(4) Res
    end function FQ
  end interface
  real(8),parameter:: eps_h = 1d-6, eps_Q = 1d-7 !误差限
  real(8) a, b, c, y1, y2, y3
  integer i
  
  a = hL; b = hR !解区间上下限
  if(Q<eps_Q) then !Q=0
    Res = a; return
  else if(Q>Par(2)-eps_Q) then !Q=Qs
    Res = b; return
  end if
  y1 = FQ( 0, real(a), Par ) - Q
  y2 = FQ( 0, real(b), Par ) - Q
  do i = 1, 10000 !最大迭代次数
    c = 0.5d0 * (a+b)
    y3 = FQ( 0, real(c), Par ) - Q
    if(abs(y3)<eps_h) exit !满足求解误差，退出循环
    if(y1*y3<0) then
      b = c
      y2 = FQ( 0, real(b), Par ) - Q
    else
      a = c
      y1 = FQ( 0, real(a), Par ) - Q
    end if
  end do
  Res = c
end function SolveFQ

!设置 Van 参数的值
pure subroutine SetPar( Par, Qr, Qs, alpha, rn, Ks, rL, rm, omega )
! SetPar( Par, Qr=Qr, Qs=Qs, alpha=alpha, rn=rn, Ks=Ks, rL=rL, rm=rm, omega=omega )
  implicit none
  real(4), intent(in):: Par(:)
  real(8),optional,intent(out):: Qr, Qs, alpha, rn, Ks, rL, rm, omega
  if(present(Qr))    Qr    = Par(1)
  if(present(Qs))    Qs    = Par(2)
  if(present(alpha)) alpha = Par(3)
  if(present(rn))    rn    = Par(4)
  if(present(Ks))    Ks    = Par(5)
  if(present(rL))    rL    = Par(6)
  if(present(rm))    rm    = (Par(4)-1d0) / Par(4)
  if(present(omega)) omega = omega0
end subroutine SetPar

!读取自定义参数
subroutine ReadUserDefinedPar(filename)
  implicit none
  character(*),intent(in):: filename
  namelist /UDP/ omega0, bSad
  integer i
  open(101,file=filename,status='old',iostat=i)
  if(i/=0) then
    omega0 = 2.5d-4
    bSad = .true.
    return
  end if
  read(101,*) omega0
  read(101,*) bSad
  close(101)
end subroutine ReadUserDefinedPar
end module FunMod_Van
  
  
  
!计算 hc 模块
module CalHc
implicit none
private
public solveHc
contains
  
!计算θ(h)
pure function FQ( h, Par ) result(Res)
  implicit none
  real(8), intent(in):: h
  real(4), intent(in):: Par(:)
  real(8) Res
  real(8) Qr, Qs, alpha, rn, rm
  real(8) a
  
  call SetPar( Par, Qr=Qr, Qs=Qs, alpha=alpha, rn=rn, rm=rm )
  if(h>=-1d-20) then 
    res = Qs; return
  end if
  a = 1d0 + (-alpha*h)**rn
  a = a**rm
  a = (Qs-Qr) / a
  res = Qr + a
end function FQ

!计算h(θ)
pure function FH( Q, Par) result(Res)
  implicit none
  real(8), intent(in):: Q
  real(4), intent(in):: Par(:)
  real(8) Res
  real(8) Qr, Qs, alpha, rn, rm
  real(8) a
  
  !直接求取
  call SetPar( Par, Qr=Qr, Qs=Qs, alpha=alpha, rn=rn, rm=rm )
  a = (Qs-Qr) / (Q-Qr)
  a = a**(1d0/rm) - 1d0
  a = a**(1d0/rn) / alpha
  Res = -a
end function FH 

!计算a
pure function Fa( Qc, Par ) result(Res)
  implicit none
  real(8), intent(in):: Qc
  real(4), intent(in):: Par(:)
  real(8) Qr, Qs, rn, rm
  real(8) Res, a
  
  call SetPar( Par, Qr=Qr, Qs=Qs, rn=rn, rm=rm )
  if(Qc<=Qr) then
    res = -1.0d37; return
  end if
  a = (Qc-Qr)/(Qs-Qr)
  a =  1d0  - a**(1d0/rm) 
  a = a * (Qc-Qr) * rm * rn
  res = -log10(exp(1d0)) / a
end function Fa

!计算残差
pure function FunHcErr( Q, xConv, Par ) result(Res)
  implicit none
  real(8), intent(in):: Q, xConv
  real(4), intent(in):: Par(:)
  real(8) Res
  real(8) Qs
  real(8) a, h
  
  call SetPar( Par, Qs=Qs )
  h = FH(Q,par)
  if(Q>Qs-1d-5) then
    res = 1d6; return
  end if
  a = Fa( Q, Par )
  res = a*Q - log10(-h) + 6.91d0-log10(xConv)
end function FunHcErr

!计算hc
 function SolveHc( Par, xConv, errVal ) result(Res)
  implicit none
  real(4), intent(in):: Par(:)
  real(8), intent(in):: xConv
  real(8), optional, intent(out):: errVal
  real(8) Res
  real(8) Qr, Qs
  integer,parameter::n=500
  real(8) err(1:n), d, Q
  integer i
  
  call SetPar( Par, Qr=Qr, Qs=Qs )
  d = 0.8*(Qs-Qr) / n
  Q = Qr
  do i = 1, n
    Q = Q + d
    err(i) = FunHcErr( Q, xConv, Par )
  end do
  i = sum( minloc(abs(err)) ) 
  
  !细化搜索
  if (i>=2) then
  	Res = Qr + (i-2)*d
  else
    Res = Qr + (i-1)*d
  end if
  Q = res
  d = 4*d / n
  do i = 1, n
    Q = Q + d
    err(i) = FunHcErr( Q, xConv, Par )
  end do
  i = sum( minloc(abs(err)) ) 
  if(present(errVal)) errVal = err(i)
  res = res + d*i
  Res = FH(res,par)
end function SolveHc

!设置 Van 参数的值
pure subroutine SetPar( Par, Qr, Qs, alpha, rn, Ks, rL, rm, omega )
  implicit none
  real(4), intent(in):: Par(:)
  real(8),optional,intent(out):: Qr, Qs, alpha, rn, Ks, rL, rm, omega
  if(present(Qr))    Qr    = Par(1)
  if(present(Qs))    Qs    = Par(2)
  if(present(alpha)) alpha = Par(3)
  if(present(rn))    rn    = Par(4)
  if(present(Ks))    Ks    = Par(5)
  if(present(rL))    rL    = Par(6)
  if(present(rm))    rm    = (Par(4)-1d0) / Par(4)
  if(present(omega)) omega = 2.5e-4
end subroutine SetPar
end module CalHc