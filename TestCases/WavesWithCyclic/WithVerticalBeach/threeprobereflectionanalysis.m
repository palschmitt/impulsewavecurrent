function [frequency, xi, xr] = threeprobereflectionanalysis(y1, y2, y3, d12,d23, h, dt);
% usage  [frequency, xi, xr] = threeprobereflectionanalysis(y1, y2, y3, d12,d23, h, dt)
% yi datasets in columns		  		Probe1			Probe2 				Probe3
% 										|-------d12-----|-------------d23---------|
% h waterdepth			---> incoming wave direction
% dt timestep
 
% constants
g = 9.81;

%number of data points
n = length(y1);

%calculate nyquist frequency
df = 1 / (n*dt);

% separation ratio
d13 = d12+d23;
mu = d12/d13;

% limit to 0.05 < kx < 0.45
kmax = 0.45 * 2 * pi / min(d12,d23);
fmax = sqrt( g * kmax * tanh(kmax * h)) / (2.0 * pi);
kmin = 0.05 * 2 * pi / d13;
fmin = sqrt( g * kmin * tanh(kmin * h)) / (2.0 * pi);

if (fmin<=0)
        error ('(fmin<0)');
    end

nfreq = fix(fmax/df);
mfreq = fix(fmin/df);
if (nfreq<mfreq)
        error ('(nfreq<mfreq)');
    end
if (mfreq<=0)
        error ('(mfreq<0)');
    end



%apply fourier transform to signals
y1_fft = fft(y1);
y2_fft = fft(y2);
y3_fft = fft(y3);

y1_fft(1) = []; %remove offset
y2_fft(1) = []; %remove offset
y3_fft(1) = []; %remove offset

options = optimset;

% calculate wave number and shoaling coefficients
for ifreq 			= 1:nfreq-mfreq+1
   %Replaced 
   %frequency(ifreq) = 2 * pi * (mfreq + ifreq - 1) * df;
   %sig = frequency(ifreq);
   %k(ifreq) = fzero(@rootk, sig^2/g, options, sig, h);
   %with this because rootk did not work
   frequency(ifreq) = 2 * pi * (mfreq + ifreq - 1) * df;
   sig 				= frequency(ifreq);
   num 				= sig^2/g;
   fun 				= @(var)abs(g*var*tanh(h*var)) - sig^2;
   k(ifreq) 		= fzero(fun, num);
 end;

kx(1,:) = k .* 0;
kx(2,:) = k .* d12;
kx(3,:) = k .* d13;

% wave component amplitudes at each wave probe position
A(:,1) = abs(y1_fft(mfreq:nfreq))/(0.5*n);
A(:,2) = abs(y2_fft(mfreq:nfreq))/(0.5*n);
A(:,3) = abs(y3_fft(mfreq:nfreq))/(0.5*n);

% wave component phases at each wave probe position
phase(:,1) = angle(y1_fft(mfreq:nfreq));
phase(:,2) = angle(y2_fft(mfreq:nfreq));
phase(:,3) = angle(y3_fft(mfreq:nfreq));

% number of frequencies
nf = nfreq-mfreq+1;

s1 = sum(exp(complex(0,2*kx)),1);
s2 = sum(exp(complex(0,-2*kx)),1);
s3 = sum(A'.*exp(complex(0,phase'+kx)),1);
s4 = sum(A'.*exp(complex(0,phase'-kx)),1);
s5 = s1 .* s2 - 9;
%Adjusted xi and xr to be in accordance with wave direction as in documentation
xr = abs((s2.*s3 - 3*s4) ./ s5);
xi = abs((s1.*s4 - 3*s3) ./ s5);
end