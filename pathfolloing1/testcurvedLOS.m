% curved path following test
% date: 11st Feb 2022
% Author: quyinsong
% Reference: Handbook of marine craft hydrodynamics and motion control
% Edition 1 : 12nd Feb 2022 
clc;
clear;
close all;
%% initial
ts = 0.01;
tfinal = 500;
Ns = tfinal/ts;
% USV 
Ustate = [0 0 0 25 1 pi/2]';
Ustate_1 = [0 0 0 22 0 0]';
Ustate_2 = [0 0 0 22 0 0]';
w = 0;
% kinematic
xxd_1 = 0;
beta_1 = 0;
psaif_1 = 0;
rd_1 = 0;
% RBF
M = 7;
What = ones(M,3);
nhat = [0 0 0]';
vhat = [0 0 0]';
fhat = [0 0 0]';
%% simulation
for k=1:1:Ns
   tout(k,1)=(k-1)*ts;
   % curved path
   K2 = 100;
   xd = 2*cos(w)+20;
   yd = 20*w;
   xd_dw = -2*sin(w);
   yd_dw = 20;
   xd_ddw = -2*cos(w);
   yd_ddw = 0;
   kc = abs(xd_dw*yd_ddw-yd_dw*xd_ddw)/sqrt(xd_dw^2+yd_dw^2)^3;  
   psaif = atan2(yd_dw,xd_dw); % psaif = w+pi/2; %  ???????????????? 
   s = cos(psaif)*(Ustate(4)-xd)+sin(psaif)*(Ustate(5)-yd);
   e = -sin(psaif)*(Ustate(4)-xd)+cos(psaif)*(Ustate(5)-yd);
   beta = atan2(vhat(2),vhat(1));
   xxsf = psaif-Ustate(6)-beta;
   Ud = sqrt(Ustate(1)^2+Ustate(2)^2)*cos(xxsf)+K2*s;
   w = ts*Ud/sqrt(xd_dw^2+yd_dw^2)+w;  
   % Kinematic
   deta = 1;
   K1 = 5;
   xxd = atan2(e,deta);
   xxddot = (xxd-xxd_1)/ts;
   xxd_1 = xxd;
   betadot = (beta-beta_1)/ts;
   beta_1 = beta;
   rd = kc*Ud-betadot+K1*(xxsf-xxd)-xxddot;
   rddot = (rd-rd_1)/ts;
   rd_1 = rd;
   % kinetic
    m = 23.8; Xudot = -2; Nvdot = 0; Iz = 1.76; Yvdot = -10; 
    Nrdot = -1; xg = 0.046; Yrdot = 0;
   % ----------------------------------------------------
    m11 = m-Xudot; 
    m22 = m-Yvdot;
    m23 = m*xg-Yrdot;
    m32 = m*xg-Nvdot;
    m33 = Iz-Nrdot;
   % -----------------------------------------------------
    u = vhat(1); v = vhat(2); r = vhat(3); frhat = fhat(3);
   % -----------------------------------------------------
    Mass=[m11  0    0;
          0   m22 m23;
          0   m32 m33];
    m0 = m22*m33-m23*m32;
   %------------------------------------------------------
   Kp = 1;
   tr =(-Kp*(r-rd)-frhat*m0/m22+rddot);
   % USV
   d = 2*randn(3,1);
   current = [0 0]';
   tao = [40 0 tr]';
   Ustatedot = USV01(Ustate,tao,current,d);
   Ustate = euler2(Ustatedot,Ustate,ts);
   Ustate_2= Ustate_1; Ustate_1 = Ustate;
   % RBF
   NNin = [Ustate(4:6)',Ustate_1(4:6)',Ustate_2(4:6)',tao']';
   N = length(NNin);
   b = 1.5*ones(M,1);
   c = 1*ones(N,M);
   for j = 1:1:M
       h(j,1) = exp(-norm(NNin-c(:,j))^2/2*b(j)^2);
   end
   NNout = What'*h;
   psai = Ustate(6);
   Rpsai = [cos(psai) -sin(psai) 0; sin(psai) cos(psai) 0; 0 0 1];
   gama = 0.01; rou = 0.005;
   What = What+ts*(gama*h*(nhat-Ustate(4:6))'*Rpsai-rou*What);
   % observer
   fhat = NNout;
   K01 = diag([1 1 1]);K02 = diag([10 10 10]);
   nhat = nhat+ts*(Rpsai*vhat-K01*(nhat-Ustate(4:6)));
   vhat = vhat+ts*(Mass\tao+fhat-K02*Rpsai'*(nhat-Ustate(4:6)));
   % out
   Ustateout(k,:) = Ustate';
   pathout(k,:) = [xd yd];
   seout(k,:) = [s e xxsf-xxd];
   taoout(k,:) = tao';
   frhatout(k,1) = frhat;
   nout(k,:) = Ustate(4:6)';
   nhatout(k,:) = nhat';
   vout(k,:) = Ustate(1:3)';
   vhatout(k,:) = vhat';
   errorout(k,:) = [Ustate(1:3)'-vhat',Ustate(4:6)'-nhat'];
end
%% plot
figure
plot(pathout(:,2),pathout(:,1),'b-',Ustateout(:,5),Ustateout(:,4),'r-','linewidth',2);
title('曲线路径跟踪效果图');
xlabel('E/m');ylabel('N/m');
figure
subplot(3,1,1); plot(tout,Ustateout(:,1),'r-','linewidth',2);title('u');xlabel('time/s');ylabel('u(m/s)');
subplot(3,1,2); plot(tout,Ustateout(:,2),'r-','linewidth',2);title('v');xlabel('time/s');ylabel('v(m/s)');
subplot(3,1,3); plot(tout,Ustateout(:,3),'r-','linewidth',2);title('r');xlabel('time/s');ylabel('r(rad/s)');
figure
subplot(3,1,1); plot(tout,Ustateout(:,4),'r-','linewidth',2);title('x');xlabel('time/s');ylabel('x(m)');
subplot(3,1,2); plot(tout,Ustateout(:,5),'r-','linewidth',2);title('y');xlabel('time/s');ylabel('y(m)');
subplot(3,1,3); plot(tout,Ustateout(:,6),'r-','linewidth',2);title('psai');xlabel('time/s');ylabel('psai(rad)');
figure
plot(tout,seout(:,1),'r-',tout,seout(:,2),'b-',tout,seout(:,3),'k-','linewidth',2);
title('跟踪误差');
xlabel('time/s');ylabel('误差/m');
legend('纵向误差s','横向误差e','角度误差xxsf~');
figure
plot(tout,taoout(:,1),'r-',tout,taoout(:,2),'g-',tout,taoout(:,3),'b-','linewidth',2);
title('控制力和力矩');
xlabel('time/s');ylabel('力/N');
legend('Tx','Ty','Tr');
figure
plot(tout,frhatout(:,1),'r-','linewidth',2);
title('frhat');
xlabel('time/s');ylabel('frhat');
figure
plot(tout,nout(:,1),'r--',tout,nout(:,2),'b--',tout,nout(:,3),'k--',tout,nhatout(:,1),'g-',...
     tout,nhatout(:,2),'m-',tout,nhatout(:,3),'y-','linewidth',2); 
title('位置估计');
legend('x','y','psai','xhat','yhat','psaihat');
xlabel('time/s'); ylabel('real x y z/ estimate xhat yhat zhat');
figure
plot(tout,vout(:,1),'r--',tout,vout(:,2),'b--',tout,vout(:,3),'k--',tout,vhatout(:,1),'g-',...
     tout,vhatout(:,2),'y-',tout,vhatout(:,3),'m-','linewidth',2);
title('速度估计');
legend('u','v','r','uhat','vhat','rhat');
xlabel('time/s'); ylabel('real u v r/ estimate uhat vhat rhat');
figure
plot(tout,errorout(:,1),'r-',tout,errorout(:,2),'b-',tout,errorout(:,3),'k-','linewidth',0.1);
title('估计误差uvr');
legend('eu','ev','er');
xlabel('time/s');ylabel('error');
figure
plot(tout,errorout(:,4),'c-',tout,errorout(:,5),'m-',tout,errorout(:,6),'y-','linewidth',0.1);
title('估计误差xypsai');
legend('ex','ey','epsai');
xlabel('time/s');ylabel('error');
