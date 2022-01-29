% line of sight guidence used in straight line path following

% control law use PID
clc
clear all
% USV parameters 
m11 = 25.8; m22 = 33.8; m33 = 2.76; m23 = 1.095; m32 =1.095;
Xu=0.72253;          Yv=-0.88965;          Nv=0.0313;
Xuu=-1.32742;        Yr=-7.25;             Nr=-1.9;
                     Yvv=-36.47297;        Nvv=3.95645;
                     Yrv=-0.805;           Nrv=0.13;
                     Yvr=-0.845;           Nvr=0.08;
                     Yrr=-3.45;            Nrr=-0.75;


% generate two points
xk =[5 5]';yk =[300 150]';
afak=atan2(yk(2)-xk(2),yk(1)-xk(1));
%initial
ts =0.01;
tfinal=50;
Ns=tfinal/ts;
x=[0.1 0.1 0 5 5 0]';x0=x;
ek_1=1; Ek(1)=ek_1;

psaid_1 = 0.1; psaid_2 = 0.05;
% simulation
disp('Simulation ...');
for k=1:1:Ns
    time(k)=(k-1)*ts;
    
    % LOS law
    Kp1=0.1;
    ye=-(x(4)-xk(1))*sin(afak)+(x(5)-xk(2))*cos(afak);
    YE(k)=ye;
    beta=atan2(x(2),x(1));
    psaid=afak+atan2(-Kp1*ye,1)-beta;
    % control law
    u = x(1);v=x(2);r= x(3);
    ek=x(6)-psaid;
    if k*ts<=0.5, Kd=0.8; else ,Kd=10;end 
    if ek<=0.05, Kp2=3;else,Kp2=1;end


    % matrix
    n11=-Xu-Xuu*abs(u);
    n22=-Yv-Yvv*abs(v)-Yrv*abs(r);
    n23=-Yr-Yvr*abs(v)-Yrr*abs(r);
    n32=-Nv-Nvv*abs(v)-Nrv*abs(r);
    n33=-Nr-Nvr*abs(v)-Nrr*abs(r);
    M=[m11  0    0;
       0   m22 m23;
      0   m32 m33];
    Crb=[0             0       -m22*v-m23*r;
         0             0        m11*u;
         m22*v+m23*r  -m11*u       0        ];
%     Nvr=[n11   0     0;
%           0   n22  n23;
%           0   n32  n33 ];
    c13=Crb(1,3);
    c23=Crb(2,3);
    c31=-c13;c32=-c23;
    m0 = m22*m33-m23*m32;
    fr = (m32*(c23*r+n22*v+n23*r)+m22*(c13*u+c23*v-n32*v-n33*r))/m0;
    psaidd = (psaid-2*psaid_1+psaid_2)/ts^2;
    psaid_2=psaid_1; psaid_1 = psaid;
    tpid=-Kp2*ek-Kd*(ek-ek_1)/ts-fr*m0/m22+psaidd; 
    ek_1=ek;
    tao=[20 0 tpid]';
    Ttao(k,:)=tao';
    Ek(k)=ek;
    % USV
    d = [1*randn(1,1) 2*randn(1,1) 2*randn(1,1)]';
    xdot=USV(x,tao,[0 0]',[0,0]',d);
    % state update
    x=euler2(xdot,x,ts);
    % store time series
    xout(1,:)=x0;
    xout(k,:)=x';
end
u=xout(:,1);
v=xout(:,2);
r=xout(:,3);
N=xout(:,4);
E=xout(:,5);
psai=xout(:,6);
% plot
disp('Plot ...');
for k=1:1:Ns
    pos =[N(k) E(k)]';
    if k==1
        modelplot(pos,psai(k));
    end
    if rem(k*ts,5)==0
        modelplot(pos,psai(k));
    end   
end
plot(E,N,'r','linewidth',2)
plot([xk(2) yk(2)],[xk(1) yk(1)],'b',E,N,'r--','linewidth',2)
hold off;
figure(2);
plot([xk(2) yk(2)],[xk(1) yk(1)],'b',E,N,'r--','linewidth',2)
xlabel('E');ylabel('N');
figure(3);
plot(time,psai*180/pi,'r','linewidth',2);
xlabel('time/s');ylabel('psai/deg');
figure(4);
plot(time,u,'r','linewidth',2)
xlabel('time/s');ylabel('u (m/s)');
figure(5);
plot(time,Ttao(:,1),'r',time,Ttao(:,3),'b','linewidth',2)
legend('surge force','yaw torch');
figure(6);
plot(time,YE,'r','linewidth',2)
xlabel('time/s');ylabel('YE (m)');

