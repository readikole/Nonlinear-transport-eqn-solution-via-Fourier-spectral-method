%% solving the Transport equation using the Chebyshev spectral 
% differentiation matrix 
clear all;close all;clc
% defining the computational parameters
N = 64; L =10; dt = 0.0001;
% computing the first order differentiation matrix via Fourier Method and
% the x value

[x,D]=fourier_matrix(N);x=(2*x-2*pi)';
Tmax=4; tplot =0.075;

plotgap= round(tplot/dt); dt = tplot/plotgap;
nplots = round(Tmax/tplot);

% defining the initial condition
u = sech(x); uold = sech(x-dt);

% Main calculation
udata = [u; zeros(nplots,N)]; t=0; tdata=0;

%% main parameter of interest
k =0;
% picking the values of k, we have deduced that the region where the
% propagating wave is when -7<k<0. Meaning the interval where the wave
% doesnot break is outside that interval


%% loop
coeff_fourier = 2*pi/(2*L);
% timing the loop
t0 = cputime;
for i =1:nplots
    for n = 1:plotgap
        t = t+dt;
        v_tilda = D*u';
        unew = uold -2*dt*sin(t).*(1./(1+k*u.^2)).*coeff_fourier.*v_tilda';
        uold =u; u=unew;
    end
    udata(i+1,:)=u;tdata =[tdata;t];
end
runningtime = cputime-t0 

%% plotting the numerical results
figure(1)
%surf(x,tdata,udata)
surf(x, tdata,udata); view(10,70)
axis([-2*pi 2*pi 0 Tmax min(min(udata)) max(max(udata))])
title('The numerical solution of the nonlinear transport equation by the Fourier spectral method')
xlabel('x');ylabel('t'); zlabel('u(x,t)');
%figure(2)
%imagesc(udata)