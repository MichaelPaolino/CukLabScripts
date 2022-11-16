function S_w = dft1D(S_t,t,w)
%Discrete fourier transform of a time domain data set into frequency domain
%S_t = data matrix of size [nt nx] where nx is an extra dim
%t = delays in S_t of size [nt 1], unit of ps
%w = freq to tranform to of size [nw 1], unit of THz
%
%S_w = output data matrix of size [nw nx]

dt = t(2)-t(1);

S_w = dt*((S_t.')*exp(-1i*2*pi*t*(w.'))).';