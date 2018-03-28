module FunMod_Kosugi
implicit none
private
public ModelNum, xConv, tConv, FS, FH, FK, FQ, FC, FKcap
integer, parameter :: ModelNum = 0        ! ģ�ͱ��
real(8), parameter :: h0 = -6.3d6         ! ��ˮ��Ϊ0ʱ��������
real(8), parameter :: hL = h0*(1d0-1d-7)  ! hȡֵ����߽�
real(8), parameter :: hR = -1d-30         ! hȡֵ���ұ߽�
real(8):: xConv=100                       ! ��תΪ Hydrus ��λ��ת��ϵ��, ��input.for-Conversion�и�ֵ
real(8):: tConv=1./(60.*60.*24.)          ! ��תΪ Hydrus ��λ��ת��ϵ��, ��input.for-Conversion�и�ֵ
  contains
!�����(h)��һ�׵���
pure function FC( iModel, h, Par ) result(Res)
  implicit none
  integer(4), intent(in):: iModel
  real(4),    intent(in):: h, Par(:)
  real(8), parameter :: x0 = log10(-h0)
  real(8) Qr, Qs, hm, ha, Ks, sigma, tau, omega
  real(8) a, b
  real(4) Res  
  if(h>=hR .OR. h<=hL) then
    Res = 0.0; return
  end if
  call SetPar( Par, Qr, Qs, hm, ha, Ks, sigma, tau, omega )
  a = Derivative_gamma( real(h,8), hm, sigma ) !gamma������h�ĵ���
  a = a / ( 1d0-gamma(h0, hm, sigma) )         !Scap��h�ĵ���
  b = Derivative_Sad( real(h,8), Qr, Qs, ha, sigma )   !Sad��h�ĵ���
  Res = (Qs-Qr)*a + Qr*b
end function FC

!����Se(��)
pure function FS( iModel, h, Par ) result(Res)
  implicit none
  integer(4), intent(in):: iModel
  real(4),    intent(in):: h, Par(:)
  real(8) Q, Qr, Qs, hm, ha, Ks, sigma, tau, omega
  real(4) Res
  call SetPar( Par, Qr, Qs, hm, ha, Ks, sigma, tau, omega )
  Q = FQ( iModel, h, Par )
  Res = (Q-Qr) / (Qs-Qr)
end function FS

!����h(��)
pure function FH( iModel, Qe, Par) result(Res)
  implicit none
  integer(4), intent(in):: iModel
  real(4),    intent(in):: Qe, Par(:)
  real(8) Q, Qr, Qs, hm, ha, Ks, sigma, tau, omega
  real(4) Res 
  call SetPar( Par, Qr, Qs, hm, ha, Ks, sigma, tau, omega )
  Q = Qe * (Qs-Qr) + Qr
  Res = SolveFQ( iModel, Q, Par, FQ )
end function FH  

!����K(h)
pure function FK( iModel, h, Par ) result(Res)
  implicit none
  integer(4), intent(in):: iModel
  real(4),    intent(in):: h, Par(:)
  real(8) Qr, Qs, hm, ha, Ks, sigma, tau, omega
  real(8) Q, s1, s2
  real(4) Res
  call SetPar( Par, Qr, Qs, hm, ha, Ks, sigma, tau, omega )
  s1 = (1d0-omega)*Kcap_rel( real(h,8), hm, sigma, tau )
  s2 = omega*Kad_rel( real(h,8), Qr, Qs, ha, sigma )
  Q = FQ( iModel, h, Par )
  Res = Ks*(s1+s2) !+ Kvap( real(h,8), Q, Qs )
  Res = max(Res,1d-38)
end function FK
!����K(h)
pure function FKcap( iModel, h, Par ) result(Res)
  implicit none
  integer(4), intent(in):: iModel
  real(4),    intent(in):: h, Par(:)
  real(8) Qr, Qs, hm, ha, Ks, sigma, tau, omega
  real(8) Q, s1, s2
  real(4) Res
  call SetPar( Par, Qr, Qs, hm, ha, Ks, sigma, tau, omega )
  s1 = (1d0-omega)*Kcap_rel( real(h,8), hm, sigma, tau )
  s2 = omega*Kad_rel( real(h,8), Qr, Qs, ha, sigma )
  Q = FQ( iModel, h, Par )
  Res = Ks*s1 !+ Kvap( real(h,8), Q, Qs )
  Res = max(Res,1d-38)
end function FKcap

!�����(h)
pure function FQ( iModel, h, Par ) result(Res)
  implicit none
  integer(4), intent(in):: iModel
  real(4),    intent(in):: h, Par(10)
  real(8) Qr, Qs, hm, ha, Ks, sigma, tau, omega
  real(4) Res
  call SetPar( Par, Qr, Qs, hm, ha, Ks, sigma, tau, omega )
  if (h<hL) then
  	Res = 1d-37
  elseif(h<hR) then
    Res = (Qs-Qr)*Scap( real(h,8), hm, sigma ) + Qr*Sad( real(h,8), Qr, Qs, ha, sigma )
    Res = max( Res, 1d-37 )    
  else
    Res = Qs
  end if
end function FQ


! ԭëϸ���Ͷȷ��� ��(h)��(Kosugi, 1996)
elemental function gamma( h, hm, sigma ) result(Res)
  implicit none
  real(8), intent(in):: h, hm, sigma
  real(8) Res
  Res = log(-h/hm) / (sqrt(2d0)*sigma) !x
  Res = 0.5d0 * erfc(Res)
end function gamma
! gamma������(h)��h�ĵ�����d��/dh
elemental function Derivative_gamma( h, hm, sigma ) result(Res)
  implicit none
  real(8), intent(in):: h, hm, sigma
  real(8), parameter :: Pi = 3.1415926d0
  real(8) Res, x
  x = log(-h/hm) / (sqrt(2d0)*sigma)
  Res = -2d0 / sqrt(Pi) * exp(-x*x) !��˹��������
  Res = 0.5d0 * Res / (h*sqrt(2.0)*sigma)
end function Derivative_gamma

!�⻬����b
elemental function SmoothingPar( Qr, Qs, sigma ) Result(Res)
  real(8), intent(in)   :: Qr, Qs, sigma
  real(8) Res
  Res = 1d0 - exp( -(Qr/(Qs-Qr))**2 )
  Res = 0.1d0 + 0.07d0*sigma*Res
end function SmoothingPar

!��ëϸ���Ͷȷ���
elemental function Scap( h, hm, sigma ) result(Res)  
  implicit none
  real(8), intent(in):: h, hm, sigma
  real(8) Res
  Res = gamma( h0, hm, sigma ) 
  Res = ( gamma( h, hm, sigma )-Res ) / ( 1d0-Res )
end function Scap

!��Ĥˮ���Ͷ�
elemental function Sad( h, Qr, Qs, ha, sigma ) result(Res)
  implicit none
  real(8), intent(in):: h, Qr, Qs, ha, sigma
  real(8), parameter :: x0 = log10(-h0)
  real(8) Res, xa, x, b
  x = log10(-h); xa = log10(ha)
  b = SmoothingPar( Qr, Qs, sigma ) !�⻬����b
  Res = 1d0 + exp((xa-x)/b)
  Res = x - xa + b*log(Res)
  Res = 1d0 + Res/(xa-x0)
end function Sad 
!��Ĥˮ���Ͷȵĵ�����dSad/dh
elemental function Derivative_Sad( h, Qr, Qs, ha, sigma ) result(Res)
  implicit none
  real(8), intent(in):: h, Qr, Qs, ha, sigma
  real(8), parameter :: x0 = log10(-h0)
  real(8) Res, xa, x, b, c
  x = log10(-h); xa = log10(ha)
  b = SmoothingPar( Qr, Qs, sigma ) !�⻬����b
  b = exp( (xa-x)/b ) ! exp()
  c = 1d0 / (h*log(10d0)) ! dx/dh, x=log10(h)
  Res = c / (xa-x0) / (1d0+b)
end function Derivative_Sad

!ëϸ��Ե�ˮ��
elemental function Kcap_rel( h, hm, sigma, tau ) result(Res)
  implicit none
  real(8), intent(in):: h, hm, sigma, tau 
  real(8) Res, s, gh0, gh, a
  a = sigma / sqrt(2d0)
  gh0 = gamma( h0, hm, sigma )
  gh  = gamma( h , hm, sigma )
  s =   InverseErfc(2d0*gh0) + a; s = erf(s)
  Res = InverseErfc(2d0*gh)  + a; Res = erf(Res)
  Res = ( (s-Res) / (1+s) )**2 !����
  s = Scap( h, hm, sigma )
  if(s<1d-37) then
    Res = 0d0
  else
    Res = s**tau * Res
  end if
end function Kcap_rel

!��Ĥ��Ե�ˮ��
elemental function Kad_rel( h, Qr, Qs, ha, sigma ) result(Res)
  implicit none
  real(8), intent(in):: h, Qr, Qs, ha, sigma 
  real(8) Res, s
  s = Sad( h, Qr, Qs, ha, sigma ); s = -1.5d0 * (1d0-s)
  Res = (-h0/ha)**s
end function Kad_rel

!��̬��ˮ��
elemental function Kvap( h, Q, Qs ) result(Res)
  implicit none
  real(8), intent(in):: h, Q, Qs
  real(8),  parameter:: temp = 20d0      !�����¶�
  real(8),  parameter:: T0 = 273.15d0    !�������
  real(8),  parameter:: M = 0.018015d0   !ˮ��Ħ������
  real(8),  parameter:: g = 9.81d0       !�������ٶ�
  real(8),  parameter:: R = 8.314d0      !ͨ�����峣��
  real(8) h_m, T, Da, Qa, tau, D, row, rosv, Hr
  real(8) Res
  h_m = h / xConv !��λת��Ϊ��(m)
  T = Temp + T0   !תΪ�������¶�  
  
  !����ˮ����ɢϵ�� D
  Qa = max(Qs-Q,0d0) !�����������
  tau= Qa**(7d0/3d0) / (Qs*Qs) !��������
  Da = 2.12d-5 * (T/T0)**2 !T�¶���ˮ����ɢϵ��
  D = tau * Qa * Da !ˮ����ɢϵ��
  !ˮ���������ܶ�
  row = 1d3 * ( 1d0 - 7.37d-6*(Temp-4d0)**2 + 3.79d-8*(Temp-4d0)**3 )
  rosv= 1d-3 * exp( 31.3716d0 - 6014.79d0/T - 7.92495d-3*T ) / T
  !���ʪ��
  Hr = exp( h_m*M*g/(R*T) )
  Res = D*rosv*M*g*Hr / (row*R*T)
  Res = Res*xConv/tConv !תΪ Hydrus ������         
end function Kvap

! ���ַ��������Ե�������, �ɦȷ���h
! Q ����ֵ
pure function SolveFQ( iModel, Q, Par, FQ ) result(Res)
  implicit none
  integer(4), intent(in):: iModel
  real(8),    intent(in):: Q
  real(4),    intent(in):: Par(10)
  interface
    pure function FQ( iModel, h, Par ) result(Res)
      implicit none
      integer(4), intent(in):: iModel
      real(4),    intent(in):: h, Par(10)
      real(4) Res
    end function FQ
  end interface
  real(8),parameter:: eps_h = 1d-6, eps_Q = 1d-7 !�����
  real(8) a, b, c, y1, y2, y3
  real(4) Res
  integer i
  a = hL; b = hR !������������
  if(Q<eps_Q) then !Q=0
    Res = a; return
  else if(Q>Par(2)-eps_Q) then !Q=Qs
    Res = b; return
  end if
  y1 = FQ( iModel, real(a), Par ) - Q
  y2 = FQ( iModel, real(b), Par ) - Q
  do i = 1, 10000 !����������
    c = 0.5d0 * (a+b)
    y3 = FQ( iModel, real(c), Par ) - Q
    if(abs(y3)<eps_h) exit !����������˳�ѭ��
    if(y1*y3<0) then
      b = c
      y2 = FQ( iModel, real(b), Par ) - Q
    else
      a = c
      y1 = FQ( iModel, real(a), Par ) - Q
    end if
  end do
  Res = c
end function SolveFQ

! ���ַ����������ķ����� erfc-1(x)
! x ����ֵ
elemental function InverseErfc( x ) result(Res)
  implicit none
  real(8),    intent(in):: x
  real(8), parameter:: eps = 1d-7
  real(8) Res, a, b, c, y1, y2, y3
  integer i
  a = -1d37; b = -a !������������
  if(x>2d0-eps) then 
    Res = a; return
  else if(x<eps) then
    Res = b; return
  end if
  y1 = erfc(a) - x
  y2 = erfc(b) - x
  do i = 1, 10000 !����������
    c = 0.5d0 * (a+b)
    y3 = erfc(c) - x
    if(abs(y3)<eps) exit !����������˳�ѭ��
    if(y1*y3<0) then
      b = c
      y2 = erfc(b) - x
    else
      a = c
      y1 = erfc(a) - x
    end if
  end do
  Res = c
end function InverseErfc

!���� Kosugi ������ֵ
pure subroutine SetPar( Par, Qr, Qs, hm, ha, Ks, sigma, tau, omega )
  implicit none
  real(4), intent(in):: Par(:)
  real(8),intent(out):: Qr, Qs, hm, ha, Ks, sigma, tau, omega
  Qr = Par(1); Qs = Par(2);  hm = Par(3); ha = Par(4); Ks= Par(5); 
  sigma = Par(6); tau = Par(7); omega = Par(8); 
end subroutine SetPar
end module FunMod_Kosugi