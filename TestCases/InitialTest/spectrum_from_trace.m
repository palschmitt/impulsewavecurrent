function [freq, amp, var, phase] = spectrum_from_trace(time_trace, dt)

%          *********************Description********************* 
% Calculates the amplitude and energy density spectra (amp & var respectively)
% for the input time trace signal sampled at dt seconds. 
% Note: 
% - time_trace should have a mean value of 0 to avoid severe low frequency
%   contamination
% - [amp] = [time_trace]
% - [var] = [time_trace]^2/[freq]
% Author: DC

% dt = 0.1;
% time = [0:1:99].*0.1;
% time_trace = 14*sin(2*pi*0.2.*time);

L = length(time_trace);
Y=fft(time_trace);
[Ywelch, f ]= pwelch(time_trace,[],[],[],1./dt);
amp = abs(Y(1:(L/2)))./(L/2);
freq = [0:1:L/2 - 1]./(L*dt);
phase = angle(Y(1:(L/2)));
df = freq(2) - freq(1);
var = (amp.^2)./(2*df);

