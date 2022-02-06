% line of sight guidence used in straight line path following

% control law use PID
clc;
clear ;
close all;
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


%% initial
% generate two points
xk =[50 0]';yk =[50 80]';
afak=atan2(yk(2)-xk(2),yk(1)-xk(1));

ts =0.01;
tfinal=20;
Ns=tfinal/ts;
x=[0 0 0 2 5 pi/2]';
ek_1=2; 

psaid_1 = 0.1; psaid_2 = 0.05;

%% simulation
disp('Simulation ...');
for k=1:1:Ns
    time(k)=k*ts;
    
    % LOS law
    deta=2.5;
    ye=-(x(4)-xk(1))*sin(afak)+(x(5)-xk(2))*cos(afak);
    YE(k)=ye;
    beta=atan2(x(2),x(1));
    psaid=afak+atan2(-ye/deta,1)-beta;
    % control law
    u = x(1);v=x(2);r= x(3);
    ek=x(6)-psaid;
    Kd = 6; 
    Kp =4;
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
    m0 = m22*m33-m23*m32;
    %------------------------------------------------------
    fr = (-m22*c31*u-m22*c32*v+m32*c23*r-(-m32*d22+m22*d32)*v-(-m32*d23+m22*d33)*r)/m0;
    
    psaidd = (psaid-2*psaid_1+psaid_2)/ts^2;
    psaid_2=psaid_1; psaid_1 = psaid;
    tpid=(-Kp*ek-Kd*(ek-ek_1)/ts-fr*m0/m22+psaidd);     
    ek_1=ek;
    tao=[20 0 tpid]';
    Ttao(k,:)=tao';
    Ek(k)=ek;
    % USV
    d = [0 0 0]';
    xdot=USV01(x,tao,[0,0]',d);
    % state update
    x=euler2(xdot,x,ts);
    % store time series
    
    xout(k,:)=x';
end

%% plot

u=xout(:,1);
v=xout(:,2);
r=xout(:,3);
N=xout(:,4);
E=xout(:,5);
psai=xout(:,6);

disp('Plot ...');
for k=1:1:Ns
    pos =[N(k) E(k)]';
    if k==1
        modelplot(pos,psai(k));
    end
    if rem(k*ts,1)==0
        modelplot(pos,psai(k));
    end   
end
plot(E,N,'r','linewidth',1)
plot([xk(2) yk(2)],[xk(1) yk(1)],'b',E,N,'r--','linewidth',1)
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

