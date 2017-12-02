%Create random trace to initialise AI learning
close all
clear all
endTime=16
timeStep=0.01
Period=0.8;
Ampmax=0.01;
delay=3;
Time=(0.0:timeStep:endTime)';
pkg load signal
%Wavegroup
Amp=[0.007 0.007 0.007];
Period=[0.7 1.15 1.3];
%phase=Period/(0.25*pi)-2*pi./Period*8
phase=Period/(0.25*pi)
Elevation=zeros(length(Time),1);
for i=1:length(Amp)
    Elevation=Elevation+Amp(i)*sin(2*pi/Period(i)*Time+phase(i));
end
%Elevation=2*Ampmax/6*sin(2*pi/Period/1.2*Time)+Ampmax/8*sin(2*pi/Period/1.5*Time)+Ampmax/3*sin(2*pi/Period/2*Time);
%Elevation=Ampmax*sin(2*pi/Period*Time)

%Ramp up over ramp
%ramptime=4;
%L=length(Time);
%t = tukeywin(L,endTime-2*ramptime);
%Elevation=t.*Elevation;
%Elevation(1:floor(ramptime/timeStep))=Elevation(1:floor(ramptime/timeStep)).*t(1:floor(ramptime/timeStep));

plot(Time,Elevation);

save('-V7','Targeteta.mat','Time','Elevation')