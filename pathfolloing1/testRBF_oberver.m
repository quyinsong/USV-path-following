% RBFNN oberver test 
% Author: Quyinsong
% date: 10 Feb 2022
% Reference: Formation control with collision avoidance for underactuated
% surface vehicles , Guoqing Xia
clc;
clear;
close all;
% initial
Ustate = [0 0 0 0 0 0 ]';
Ustate_1 = [0 0 0 0 0 0]';
Ustate_2 = [0 0 0 0 0 0]';
tao = [30 0 5]';
current = [0 0]';
M = 7;
What = ones(M,3);
nhat = [0 0 0]';
vhat = [0 0 0]';
ts = 0.01;
tfinal = 50;
Ns = tfinal/ts;
% simulation
for k = 1:1:Ns
   tout(k,1) = (k-1)*ts;
   % USV
   d = [sin(k*ts) sin(k*ts) sin(k*ts)]';
   tao = [50 0 10*sin(k*ts)]';
   Ustatedot = USV01(Ustate,tao,current,d);
   Ustate = euler2(Ustatedot,Ustate,ts);
   Ustate_2= Ustate_1; Ustate_1 = Ustate;
   % USV parameters 
   u = Ustate(1); v = Ustate(2); r = Ustate(3);
   tr = tao(3);
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
    detM2=m22*m33-m23*m32;
    m0=detM2; 
    Mass=[m11  0    0;
        0   m22 m23;
        0   m32 m33];
    fu = (-c13*r-d11*u)/m11+d(1);
    fv = (m23*c31*u+m23*c32*v-m33*c23*r-(m33*d22-m23*d32)*v-(m33*d23-m23*d33)*r-m23*tr)/m0+d(2);
    fr = (-m22*c31*u-m22*c32*v+m32*c23*r-(-m32*d22+m22*d32)*v-(-m32*d23+m22*d33)*r)/m0+d(3);
    freal = [fu fv fr]';
   % RBF
   NNin = [Ustate(4:6)',Ustate_1(4:6)',Ustate_2(4:6)',tao']';
   N = length(NNin);
   b = 1.5*ones(M,1);
   c = 10*ones(N,M);
   for j = 1:1:M
       h(j,1) = exp(-norm(NNin-c(:,j))^2/2*b(j)^2);
   end
   NNout = What'*h;
   psai = Ustate(6);
   Rpsai = [cos(psai) -sin(psai) 0; sin(psai) cos(psai) 0; 0 0 1];
   gama = 1; rou = 5;
   What = What+ts*(gama*h*(nhat-Ustate(4:6))'*Rpsai-rou*What);
   % observer
   fhat = NNout;
   K01 = diag([1 1 1]);K02 = diag([100 100 100]);
   nhat = nhat+ts*(Rpsai*vhat-K01*(nhat-Ustate(4:6)));
   vhat = vhat+ts*(Mass\tao+fhat-K02*Rpsai'*(nhat-Ustate(4:6)));
   % output
   frealout(k,:) = freal';
   fhatout(k,:) = fhat';
   nout(k,:) = Ustate(4:6)';
   nhatout(k,:) = nhat';
   vout(k,:) = Ustate(1:3)';
   vhatout(k,:) = vhat';
   errorout(k,:) = [Ustate(1:3)'-vhat',Ustate(4:6)'-nhat'];
end
% plot
figure(1)
plot(nout(:,2),nout(:,1),'r-','linewidth',1.5);
title('船舶运动状态');
xlabel('E/m');ylabel('N/m');
figure(2)
plot(tout,nout(:,1),'r--',tout,nout(:,2),'b--',tout,nout(:,3),'k--',tout,nhatout(:,1),'r-',...
     tout,nhatout(:,2),'b-',tout,nhatout(:,3),'k-','linewidth',1.5); 
title('位置估计');
legend('x','y','psai','xhat','yhat','psaihat');
xlabel('time/s'); ylabel('real x y z/ estimate xhat yhat zhat');
figure(3)
plot(tout,vout(:,1),'r--',tout,vout(:,2),'b--',tout,vout(:,3),'k--',tout,vhatout(:,1),'r-',...
     tout,vhatout(:,2),'b-',tout,vhatout(:,3),'k-','linewidth',1.5);
title('速度估计');
legend('u','v','r','uhat','vhat','rhat');
xlabel('time/s'); ylabel('real u v r/ estimate uhat vhat rhat');
figure(4)
plot(tout,frealout(:,1),'b-',tout,fhatout(:,1),'r-','linewidth',1.5);
title('系统未建模动态估计');
legend('fu','fuhat');
xlabel('time/s');ylabel('fu/fuhat');
figure(5)
plot(tout,frealout(:,2),'b-', tout,fhatout(:,2),'r-','linewidth',1.5);
title('系统未建模动态估计');
legend('fv','fvhat');
xlabel('time/s');ylabel('fv/fvhat');
figure(6)
plot(tout,frealout(:,3),'b-',tout,fhatout(:,3),'r-','linewidth',1.5);
title('系统未建模动态估计');
legend('fr','frhat');
xlabel('time/s');ylabel('fr/frhat');
figure(7)
plot(tout,errorout(:,1),'r-',tout,errorout(:,2),'b-',tout,errorout(:,3),'k-',...
     tout,errorout(:,4),'c-',tout,errorout(:,5),'m-',tout,errorout(:,6),'y-');
title('估计误差');
legend('eu','ev','er','ex','ey','ez');
xlabel('time/s');ylabel('error');
