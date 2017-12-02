close all
clear all
%% define variables


timeStep=1.1/100;
%% Load Theory
% load('TheoryEta.mat')
etaResult=load('Elevation.mat','-ascii');
%% Load Data I
%dirName = strcat('.','/surfaceElevation/0.0119048/');
Depth=mean(etaResult(:,2));
WP1=etaResult(:,2);
WP2=etaResult(:,3);
WP3=etaResult(:,4);
WP4=etaResult(:,5);
WP1=interp1(etaResult(:,1),WP1,(timeStep:timeStep:timeStep*(length(etaResult(:,1)))));
WP2=interp1(etaResult(:,1),WP2,(timeStep:timeStep:timeStep*(length(etaResult(:,1)))));
WP3=interp1(etaResult(:,1),WP3,(timeStep:timeStep:timeStep*(length(etaResult(:,1)))));
WP4=interp1(etaResult(:,1),WP4,(timeStep:timeStep:timeStep*(length(etaResult(:,1)))));

%Replace NaN at start with mean
WP1(isnan(WP1))=Depth;
WP2(isnan(WP2))=Depth;
WP3(isnan(WP3))=Depth;
WP4(isnan(WP4))=Depth;
etaResult=etaResult-Depth;
%%
[frequency, xi, xr] = threeprobereflectionanalysis(WP3,WP2,WP4,0.15,0.225,Depth,timeStep);

figure
plot(frequency,xr)
title('Reflected Wave')

figure
plot(frequency,xi)
title('Incident Wave')
[C,I] = max(xi);
peak_freq = frequency(I)/(2*pi);
amp = C;
reflection_coeff = xr(I)/xi(I)
%% 
