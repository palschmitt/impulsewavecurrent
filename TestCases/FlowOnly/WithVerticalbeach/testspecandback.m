clear all
close all

addpath('~/octavescripts/matlab2tikz/src/')
addpath('/home/pal/OpenFOAM/pal-4.x/FrequMethod/DCSeaCalibration/functions/')
timeStep=0.1;
Tol=1E-5; %RMS Error of surface elevation

%Bound crazy predictions
bound=100;
probepos=2.75;
Nattempts=4;
%Load targetdata
%load('smalltestwave.mat')
load('smallregular.mat')
Time1=Time;
Elevation1=Elevation-mean(Elevation);
%load('Targetdata.mat');

Time=(min(Time1):timeStep:max(Time1))-min(Time1);
targetEta=(interp1(Time1-min(Time1),Elevation1,Time));

[targetetafreq, targetetaamp, targetetavar, targetetaphase] = spectrum_from_trace(targetEta, timeStep);

Ampmat=zeros(length(targetetaamp),length(targetetaamp));
Ampmat(:,1)=targetetaamp;
phimat=zeros(length(targetetaphase),length(targetetaphase));
phimat(:,1)=targetetaphase;


L = 2*length(targetetaamp);
%df=diff(targetetafreq(1:2))
df = 1/(2*timeStep*(L));
[time_trace, dt] = trace_from_spectrum(Ampmat, phimat, df);
etanew=time_trace(:,1);
etanew(isnan(etanew))=0.;

 figure
 subplot (3, 1, 1)
plot(Time,targetEta,'r')
hold on
%newtime=(0:dt:dt*(length(etanew)-1));
newtime=(0:dt:dt*(length(etanew)-1));
plot(newtime,etanew,'b')
legend('Orig','FromSpec')
subplot (3, 1, 2)
plot(targetetafreq,targetetaamp)
subplot (3, 1, 3)
plot(targetetafreq,targetetaphase)
