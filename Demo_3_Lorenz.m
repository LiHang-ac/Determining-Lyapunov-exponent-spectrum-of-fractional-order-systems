%% README
%
% Program to compute Lyapunov exponents of Fractional-Order Systems (FOS)
% defined by Gr¨¹nwald¨CLetnikov (G-L) derivative.
%
% Example: fractional-order Lorenz system.
% 
% Output:
% x,y,z - system responses;
% t - Time series of system responses;
% LE - Lyapunov Exponents;
% T - Time series of Lyapunov exponents;
%
% Author: Hang Li
% Email: Lihang.ac@hotmail.com
% 
% Cite as:
% Hang Li, Yongjun Shen, et al. Determining Lyapunov exponents of 
% fractional-order systems: A general method based on memory principle, 
% Chaos, Solitons & Fractals, 2023, 168: 113167. 
% Doi:10.1016/j.chaos.2023.113167

%% system parameters and simulation conditions
clear;clc;close all;
SIGMA = 16;
RHO = 45.92;
BETA = 4;


h=5e-4;
h_norm=10*h;
N=h_norm/h;

tn=100-h;
t=0:h:tn;
n=length(t);
T=0:h_norm:tn;

%% define the order
q1=0.985;q2=0.99;q3=0.98;
% q1=0.985;q2=0.97;q3=0.98;
% q1=1;q2=1;q3=1;
%% fractional-order binomial coefficient
cp1=1; cp2=1; cp3=1;
for j=1:n
    c1(j)=(1-(1+q1)/j)*cp1;
    c2(j)=(1-(1+q2)/j)*cp2;
    c3(j)=(1-(1+q3)/j)*cp3;
    cp1=c1(j); cp2=c2(j); cp3=c3(j);
end

%% initialization
x(1) = 20;
y(1) = 10;
z(1) = 50;

k=1;
f11(k)=1;         f12(k)=0;        f13(k)=0;
f21(k)=0;         f22(k)=1;        f23(k)=0;
f31(k)=0;         f32(k)=0;        f33(k)=1;

J=[ f11(k),f12(k),f13(k);...
    f21(k),f22(k),f23(k);...
    f31(k),f32(k),f33(k)];

SUM=[0;0;0];
LE=[];

%% main loop
for k=2:n
    x(k)=(SIGMA*(y(k-1) - x(k-1)))*h^q1-calmem(x,c1,k);
    y(k)=(RHO*x(k) - y(k-1) -x(k)*z(k-1))*h^q2-calmem(y,c2,k);
    z(k)=(x(k)*y(k) - BETA*z(k-1))*h^q3-calmem(z,c3,k);
    
    f11(k)=(SIGMA*(f21(k-1) - f11(k-1)))*h^q1-calmem(f11,c1,k);
    f12(k)=(SIGMA*(f22(k-1) - f12(k-1)))*h^q1-calmem(f12,c1,k);
    f13(k)=(SIGMA*(f23(k-1) - f13(k-1)))*h^q1-calmem(f13,c1,k);
    
    f21(k)=(RHO*f11(k-1) - f21(k-1) -f11(k-1)*z(k-1)-x(k-1)*f31(k-1))*h^q2-calmem(f21,c2,k);
    f22(k)=(RHO*f12(k-1) - f22(k-1) -f12(k-1)*z(k-1)-x(k-1)*f32(k-1))*h^q2-calmem(f22,c2,k);
    f23(k)=(RHO*f13(k-1) - f23(k-1) -f13(k-1)*z(k-1)-x(k-1)*f33(k-1))*h^q2-calmem(f23,c2,k);
    
    f31(k)=(f11(k-1)*y(k-1)+x(k-1)*f21(k-1) - BETA*f31(k-1))*h^q3-calmem(f31,c3,k);
    f32(k)=(f12(k-1)*y(k-1)+x(k-1)*f22(k-1) - BETA*f32(k-1))*h^q3-calmem(f32,c3,k);
    f33(k)=(f13(k-1)*y(k-1)+x(k-1)*f23(k-1) - BETA*f33(k-1))*h^q3-calmem(f33,c3,k);
    
    J=[ f11(k),f12(k),f13(k);...
        f21(k),f22(k),f23(k);...
        f31(k),f32(k),f33(k)];
    
    if mod(k,N)==0                  %% G-S orthonormalization procedure
        [J,E]=GSR(J);
        SUM=SUM+log(E);
        f11(k)=J(1,1);         f12(k)=J(1,2);        f13(k)=J(1,3);
        f21(k)=J(2,1);         f22(k)=J(2,2);        f23(k)=J(2,3);
        f31(k)=J(3,1);         f32(k)=J(3,2);        f33(k)=J(3,3);
       
        LE=[LE,SUM/(k*h)];
        
        fprintf('completed: %3.2f%%\n',k/n*100)
        fprintf('LE=%2.3f  ,  %2.3f  ,  %2.3f\n',LE(:,end));
    end
    
end
LE(:,end)

%% FIGURE
figure(1)
plot3(x,y,z,'k','LineWidth',0.5);

figure(2)
plot(T,LE,'LineWidth',1);

figure(3)
plot(t,x,'k-','LineWidth',0.1);


