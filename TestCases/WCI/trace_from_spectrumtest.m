%function [time_trace, dt] = trace_from_spectrum(amp, phase, df)

%          *********************Description********************* 
% Calculates the time trace and sample interval from the input spectral 
% amplitude and phase components
% Note: 
% - [df] = Hz
% Author: DC
pkg load signal

 amp = zeros(64);
 phase = zeros(64);
 df = 0.2;
 freq = [0:1:length(amp)].*df;
 freq(5)
 amp(5) = 12;

 
 
 
L = 2*length(amp);
freq = [0:1:(L/2) - 1]*df;

Z = (L/2).*amp.*(cos(phase) + sqrt(-1).*(sin(phase)));
for i = L/2 + 1:L 
    j = L - (i - 1);
    Z(i) = Z(j)';    
end

dt = 1/(2*df*(L));
time_trace = real(ifft(Z));
plot(time_trace(:,1))