%% README
%
% Program to compute Lyapunov exponents of Fractional-Order Systems (FOS)
% defined by Gr¨¹nwald¨CLetnikov (G-L) derivative.
%
% Example: 4D fractional-order Chen system with low effective order.
% 
% Output:
% x1,x2,x3,x4 - system responses;
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
a=36;
b=3;
c=28;
d=-16;
ka=0.5;

h=1e-3;
h_norm=10*h;
N=h_norm/h;

tn=40-h;
t=0:h:tn;
n=length(t);
T=0:h_norm:tn;
%% define the order
Q=0.95;

% Q=0.88;

% Q=0.78;

q1=Q;q2=Q;q3=Q;q4=Q;

%% Fractional-order binomial coefficient
cp1=1; cp2=1; cp3=1; cp4=1;
for j=1:n
    c1(j)=(1-(1+q1)/j)*cp1;
    c2(j)=(1-(1+q2)/j)*cp2;
    c3(j)=(1-(1+q3)/j)*cp3;
    c4(j)=(1-(1+q4)/j)*cp4;
    cp1=c1(j); cp2=c2(j); cp3=c3(j); cp4=c4(j);
end

%% initialization
x1(1) = 0;
x2(1) = 0;
x3(1) = 8;
x4(1) = 6;

k=1;
f11(k)=1;         f12(k)=0;        f13(k)=0;        f14(k)=0;
f21(k)=0;         f22(k)=1;        f23(k)=0;        f24(k)=0;
f31(k)=0;         f32(k)=0;        f33(k)=1;        f34(k)=0;
f41(k)=0;         f42(k)=0;        f43(k)=0;        f44(k)=1;


J=[ f11(k),f12(k),f13(k),f14(k);...
    f21(k),f22(k),f23(k),f24(k);...
    f31(k),f32(k),f33(k),f34(k);...
    f41(k),f42(k),f43(k),f44(k)];

SUM=[0;0;0;0;];
LE=[];

%% main loop
for k=2:n
    x1(k)=(a*(x2(k-1)-x1(k-1)))*h^q1-calmem(x1,c1,k);
    x2(k)=(d*x1(k)-x1(k)*x3(k-1)+c*x2(k-1)-x4(k-1))*h^q2-calmem(x2,c2,k);
    x3(k)=(x1(k)*x2(k)-b*x3(k-1))*h^q3-calmem(x3,c3,k);
    x4(k)=(x1(k)+ka)*h^q4-calmem(x4,c4,k);
    
    f11(k)=(a*(f21(k-1)-f11(k-1)))*h^q1-calmem(f11,c1,k);
    f12(k)=(a*(f22(k-1)-f12(k-1)))*h^q1-calmem(f12,c1,k);
    f13(k)=(a*(f23(k-1)-f13(k-1)))*h^q1-calmem(f13,c1,k);
    f14(k)=(a*(f24(k-1)-f14(k-1)))*h^q1-calmem(f14,c1,k);
    
    f21(k)=(d*f11(k)-f11(k)*x3(k-1)-x1(k)*f31(k-1)+c*f21(k-1)-f41(k-1))*h^q2-calmem(f21,c2,k);
    f22(k)=(d*f12(k)-f12(k)*x3(k-1)-x1(k)*f32(k-1)+c*f22(k-1)-f42(k-1))*h^q2-calmem(f22,c2,k);
    f23(k)=(d*f13(k)-f13(k)*x3(k-1)-x1(k)*f33(k-1)+c*f23(k-1)-f43(k-1))*h^q2-calmem(f23,c2,k);
    f24(k)=(d*f14(k)-f14(k)*x3(k-1)-x1(k)*f34(k-1)+c*f24(k-1)-f44(k-1))*h^q2-calmem(f24,c2,k);
    
    f31(k)=(f11(k)*x2(k)+x1(k)*f21(k)-b*f31(k-1))*h^q3-calmem(f31,c3,k);
    f32(k)=(f12(k)*x2(k)+x1(k)*f22(k)-b*f32(k-1))*h^q3-calmem(f32,c3,k);
    f33(k)=(f13(k)*x2(k)+x1(k)*f23(k)-b*f33(k-1))*h^q3-calmem(f33,c3,k);
    f34(k)=(f14(k)*x2(k)+x1(k)*f24(k)-b*f34(k-1))*h^q3-calmem(f34,c3,k);
    
    f41(k)=(f11(k))*h^q4-calmem(f41,c4,k);
    f42(k)=(f12(k))*h^q4-calmem(f42,c4,k);
    f43(k)=(f13(k))*h^q4-calmem(f43,c4,k);
    f44(k)=(f14(k))*h^q4-calmem(f44,c4,k);
    
    
    J=[ f11(k),f12(k),f13(k),f14(k);...
        f21(k),f22(k),f23(k),f24(k);...
        f31(k),f32(k),f33(k),f34(k);...
        f41(k),f42(k),f43(k),f44(k)];
    
    if mod(k,N)==0
        [J,E]=GSR(J);
        SUM=SUM+log(E);
        f11(k)=J(1,1);         f12(k)=J(1,2);        f13(k)=J(1,3);        f14(k)=J(1,4);
        f21(k)=J(2,1);         f22(k)=J(2,2);        f23(k)=J(2,3);        f24(k)=J(2,4);
        f31(k)=J(3,1);         f32(k)=J(3,2);        f33(k)=J(3,3);        f34(k)=J(3,4);
        f41(k)=J(4,1);         f42(k)=J(4,2);        f43(k)=J(4,3);        f44(k)=J(4,4);
        

        LE=[LE,SUM/(k*h)];
        
        fprintf('completed: %3.2f%%\n',k/n*100)
        fprintf('LLE=%2.3f  ,  %2.3f  ,  %2.3f  ,  %2.3f\n',LE(:,end));
    end
    
end
LE(:,end)

%% FIGURE
figure(1)
plot3(x1,x2,x3,'LineWidth',1);

figure(2)
plot(T,LE,'LineWidth',1);

figure(3)
plot(t,x2,'LineWidth',1);

