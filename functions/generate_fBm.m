clear all
clc

N = 2^22;
dt = 1;

H = 1.2;


[W] = fbm1d(H,N);
W = [0;diff(W)];

VH = gamma(1-2*H)*cos(pi*H)/(pi*H);

%W = VH*W;
t = [0:N-1]*dt;

figure(1);clf;hold all;
plot(t,W,'b-')

f = [-N/2 : N/2-1]*1/(N*dt);
PSD = @(x) dt/N * abs(fftshift(fft(x))).^2;

figure(2);clf;hold all;
plot(f,PSD(W),'b-')
plot(f,0.001*f.^(-2*H+1),'r-')
set(gca,'XScale','log')
set(gca,'YScale','log')