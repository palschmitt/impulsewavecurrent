%Create random trace to initialise AI learning
close all
clear all
endTime=30
timeStep=1.1/100;%0.01
delay=0 % Add 0 to end for creating delay data
Umax=1
Period=2;
pkg load signal
t=(0.0:timeStep:endTime-delay)';
%addedtime=(endTime-delay+timeStep:timeStep:endTime)';
%U=Umax.*(max(t)/2-abs((max(t)/2-t)))./max(t).*sin(2*pi/Period*t);

%U=Umax.*sin(2*pi/Period*t);
% irregular waves
%U=-2*Umax/6*sin(2*pi/Period/1.2*t)+Umax/8*sin(2*pi/Period/0.7*t)+Umax/3*sin(2*pi/Period/2*t)+Umax/6*sin(2*pi/Period/1.*t);
U=ones(ceil(endTime/timeStep),1)*2;
for i=1:200
  U(i)=U(i)*t(i)/t(200);
end
%U=[U; zeros(size(addedtime))];
%t=[t; addedtime];

plot(U);
U=num2cell(U);
t=num2cell(t);
Output=[cell2mat(t),cell2mat(U),zeros(size(cell2mat(U))),zeros(size(cell2mat(U))) ];
plot(cell2mat(t),cell2mat(U))
fid = fopen('Wavemaker.dat', 'w');
fprintf(fid,'%i\n', size(t,1));
fprintf(fid,'(\n')
fprintf(fid,'(%f ( %f %f %f ))\n',Output' )
fprintf(fid,')\n')
fclose(fid)
save('-V7','Initialwavemakerinput.mat','t','U')