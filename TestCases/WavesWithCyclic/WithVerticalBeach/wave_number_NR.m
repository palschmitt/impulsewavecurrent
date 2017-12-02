function [k, L] = wave_number_NR(w, h)

g = 9.81;
k_1 = (w^2)/g; % use deep water wavenumber as first estimate

k = fzero(@(x) (w^2) - g*x*tanh(x*h), k_1); % Newton-Raphson
k = abs(k);

L = 2*pi/k; % wavelength
