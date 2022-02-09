function xdot = USV01( x,tao,current,d )
% USV01 xdot = USV( x,tao,taod ) returns the time derivative of 
% the state vector: x = [ u v r x y psi]'  for USV, where
% INPUT: 
% u=x(1): surge velocity (m/s)
% v=x(2): sway velocity (m/s)
% r=x(3): yaw velocity (rad/s)
% x=x(4): position in {n} x-direction (m)
% y=x(5): position in {n} y-direction (m)
% psai=x(6): yaw angle (rad)
% tao=[tx ty tn]':
% wind=[Vw betaw]'
% current=[Vc betac]'
% OUTPUT: 
% xdot=[udot vdot rdot xdot ydot psaidot Vcdot]':time derivative of state vector

% Author: Quyinsong
% Data: 13rd Jan 2022
% Reference: Improved line-of-sight trajectory tracking control of
% under-actuated AUV subjects to ocean currents and input saturation

% check input and state dimentions
if nargin ~=4,error('input number must be 4!');end
if length(x) ~=6,error('state x number must be 6!');end
if length(tao) ~=3,error('ctr input tao number must be 3!');end
if length(current) ~=2,error('current input tao number must be 2!');end
if length(d) ~=3,error('diturbance taod number must be 3!');end

%% USV state:
u=x(1);
v=x(2);
r=x(3);
V = [u v r]';
psai=x(6);
tu = tao(1);
tr = tao(3);
%% curremt
Vx=current(1);
Vy= current(2);
Vc = [Vx Vy 0]';
%% diturbance
fdu = d(1);
fdv = d(2);
fdr = d(3);
%% USV parameters 
m = 23.8; Xudot = -2; Nvdot = 0; Iz = 1.76; Yvdot = -10; 
Nrdot = -1; xg = 0.046; Yrdot = 0;
% ------------------------------------------------------
Xu=-0.72253;         Yv=-0.88965;          Nv=0.0313;
Xuu=-1.32742;        Yr=-7.25;             Nr=-1.900;
                     Yvv=-36.47287;        Nvv=3.95645;
                     Yrv=-0.805;           Nrv=0.130;
                     Yvr=-0.845;           Nvr=0.080;
                     Yrr=-3.45;            Nrr=-0.75;               
% ----------------------------------------------------
m11 = m-Xudot; 
m22 = m-Yvdot;
m23 = m*xg-Yrdot;
m32 = m*xg-Nvdot;
m33 = Iz-Nrdot;
% -----------------------------------------------------
c13 = -m*(xg*r+v); 
c23 = m*u;
c31 = -c13; c32 = -c23;
% -----------------------------------------------------
d11=-Xu-Xuu*abs(u);
d22=-Yv-Yvv*abs(v)-Yrv*abs(r);
d23=-Yr-Yvr*abs(v)-Yrr*abs(r);
d32=-Nv-Nvv*abs(v)-Nrv*abs(r);
d33=-Nr-Nvr*abs(v)-Nrr*abs(r);


%% matrix expression

Rbn=[ cos(psai) -sin(psai) 0;
      sin(psai) cos(psai)  0;
      0          0         1];
M=[m11  0    0;
   0   m22 m23;
   0   m32 m33];
Crb=[0     0     c13;
     0     0     c23;
    c31   c32     0 ];
Dvr=[d11  0     0;
     0    d22  d23;
     0    d32  d33];
% time derivatives

Vdot=M\(-Crb*V-Dvr*V+tao)+d ;
Xdot=Rbn*V+Vc;
xdot=[Vdot ;Xdot];

%% components expression

% detM2=m22*m33-m23*m32;
% m0=detM2;
% 
% fu = (-c13*r-d11*u)/m11;
% fv = (m23*c31*u+m23*c32*v-m33*c23*r-(m33*d22-m23*d32)*v-(m33*d23-m23*d33)*r-m23*tr)/m0;
% fr = (-m22*c31*u-m22*c32*v+m32*c23*r-(-m32*d22+m22*d32)*v-(-m32*d23+m22*d33)*r)/m0;
% 
% % time derivatives
% 
% xdot = [fu+tu/m11+fdu;
%         fv+fdv;
%         fr+tr*m22/m0+fdr;
%         u*cos(psai)-v*sin(psai)+Vx;
%         u*sin(psai)+v*cos(psai)+Vy;
%         r];

end



