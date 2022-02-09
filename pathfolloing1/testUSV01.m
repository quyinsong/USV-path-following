% Author: Quyinsong
% Data: 14th Jan 2022
% test the USV function
clc 
clear all
close all
%% initial
ts=0.1;                 % sample time
tfinal =100;             % simulation final time
Ns =tfinal/ts;          % step number of simulation

current=[0 0]';    % current

tao=[60 0 0]'; 
tao0=tao;
d=[0 0 0]';
x=[0 0 0 40 0 pi/2]';
x0=x;
%% simulation start
disp('Simulation ... ');
for k=1:1:Ns
    time(1)=0;
    time(k+1)=k*ts;
    
    if x(6)*180/pi>=360
        x(6)=x(6)-2*pi;
    end
    if x(6)*180/pi<=-360
        x(6)=x(6)+2*pi;
    end
    
    if k*ts>=10 && k*ts<30
       tao = [60 0 5]';
    end
    if k*ts>=30 && k*ts<40
       tao = [60 0 0]';
    end
    if k*ts>=40
        tao = [60 0 8]';
    end

    Ttao(1,:)=tao0';
    Ttao(k+1,:)=tao';
    % time derivatives
    xdot=USV01(x,tao,current,d);
    % update states
    x=euler2(xdot,x,ts);
    % store time series
    xout(1,:)=x0;
    xout(k+1,:)=x';
    
end

%% plot
u=xout(:,1);
v=xout(:,2);
r=xout(:,3);
N=xout(:,4);
E=xout(:,5);
psai=xout(:,6);

disp('plot ...');
for k=1:1:Ns
    pos =[N(k) E(k)]';
    if k==1
        modelplot(pos,psai(k));
    end
    if rem(k,10)==0
        modelplot(pos,psai(k));
    end   
end
plot(E,N,'r','linewidth',2)
hold off;
figure(2);
plot(time,psai*180/pi,'r','linewidth',2);
xlabel('time/s');ylabel('psai/deg');
figure(3);
plot(E,N,'r','linewidth',2)
xlabel('E');ylabel('N');
figure(4);
plot(time,u,'r','linewidth',2)
xlabel('time/s');ylabel('u (m/s)');
figure(5);
plot(time,Ttao(:,1),'r',time,Ttao(:,2),'k',time,Ttao(:,3),'b','linewidth',2)
xlabel('time/s');ylabel('T (m/s)');


    