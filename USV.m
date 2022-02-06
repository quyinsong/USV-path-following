function xdot = USV( x,tao,wind,current,d )
% USV xdot = USV( x,tao,taod ) returns the time derivative of 
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
% model equitions: kinetics: M*vdot+Crb*v+Nvr*vr=tao+wind+d  (vr=v-vc)
%                   xdot=Rbn*v

% check input and state dimentions
if nargin ~=5,error('input number must be 5!');end
if length(x) ~=6,error('state x number must be 6!');end
if length(tao) ~=3,error('ctr input tao number must be 3!');end
if length(wind) ~=2,error('wind input tao number must be 2!');end
if length(current) ~=2,error('current input tao number must be 2!');end
if length(d) ~=3,error('diturbance taod number must be 3!');end

% USV state:
u=x(1);
v=x(2);
r=x(3);
psai=x(6);
V=[u v r]';
% transfer matrix
Rbn=[ cos(psai) -sin(psai) 0;
      sin(psai) cos(psai)  0;
      0          0         1];
% wind
Vw=wind(1);
betaw=wind(2);
uwb=Vw*cos(betaw-psai);
vwb=Vw*sin(betaw-psai);
urw=u-uwb; vrw=v-vwb;
Vrw=sqrt(urw^2+vrw^2);
gamaw=-atan2(vrw,urw);
twind=1/2*0.01*Vrw^2*[-0.5*cos(gamaw)*1 
                      -0.7*sin(gamaw)*1 
                      -0.05*sin(2*gamaw)*1];
% curremt
Vc=current(1);
betac= current(2);
Vcn=Vc*[cos(betac) sin(betac) 0]';
Vcb=Rbn'*Vcn;
ucb=Vcb(1);
vcb=Vcb(2);
ur=u-ucb;
vr=v-vcb;
Vcr=[ur vr r]';
% USV parameters 
m11 = 25.8; m22 = 33.8; m33 = 2.76; m23 = 1.095; m32 =1.095;
Xu=0.72253;          Yv=-0.88965;          Nv=0.0313;
Xuu=-1.32742;        Yr=-7.25;             Nr=-1.9;
                     Yvv=-36.47297;        Nvv=3.95645;
                     Yrv=-0.805;           Nrv=0.13;
                     Yvr=-0.845;           Nvr=0.08;
                     Yrr=-3.45;            Nrr=-0.75;
n11=-Xu-Xuu*abs(ur);
n22=-Yv-Yvv*abs(vr)-Yrv*abs(r);
n23=-Yr-Yvr*abs(vr)-Yrr*abs(r);
n32=-Nv-Nvv*abs(vr)-Nrv*abs(r);
n33=-Nr-Nvr*abs(vr)-Nrr*abs(r);
% matrix
M=[m11  0    0;
   0   m22 m23;
   0   m32 m33];
Crb=[0             0       -m22*v-m23*r;
     0             0        m11*u;
     m22*v+m23*r  -m11*u       0        ];
Nvr=[n11  0     0;
     0    n22  n23;
     0    n32  n33];
% yaw force limit
% tu = tao(1);
% tv = tao(2);
% tr = tao(3);
% if tr >= 100000
%     tr_limit = 100000;
% else
%     tr_limit = tr;
% end
% tao = [tu tv tr_limit]';
% time derivatives

Vdot=M\(-Crb*V-Nvr*Vcr+tao+d) ;
Xdot=Rbn*V;
xdot=[Vdot ;Xdot];

% components expression
% detM2=m22*m33-m23*m32;
% xdot=[ ];

end

