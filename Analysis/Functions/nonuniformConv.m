function vOut = nonuniformConv(vIn,tIn,filterIn,dt)
% NONUNIFORMCONV convolves input data on a non-uniform grid with an input
% filter. This function works by interpolating input data on a uniform
% grid, performing the convolution, and interpolating back onto the
% non-uniform grid.
%
% vOut = nonuniformConv(vIn,tIn,filterIn,dt)
%   vOut [nt,1] is convolution of vIn [nt,1] on nonuniform grid tIn [nt,1] 
%   with uniformly sampled filterIn [t',1] with step dt (scalar).
%
%Note: 
%-This function requires vectors as inputs to work correctly.
%-filterIn [nt',1] MUST be centered around 0 and evaluated on a uniform 
%   grid (-t'):dt:(+t').
%-This function performs slowly when tIn spans too many orders of magnitude
%   of dt, i.e. when you have dt in the 10s of fs range but data out to us.
%   In this case, split your data into an early and late portion and use 
%   this function on the early portion.
%
% See Also: interp1, conv

nt = length(filterIn); %length of the input filter
n = sum(filterIn); %normalization for convolution

%uniformly spaced grid with spacing dt. First value of tConv and tIn are the same
tConv = min(tIn):dt:(max(tIn)+nt*dt); 

%This does the convolution. Index managment is implicit and not very obvious... 
%Index managment is done via interpolating over the time grid. A few key points:
%1) tConv grid is generated to contain the delay range tIn plus the filter length
%2) 'extrap' allow for vIn to be linearly extrapolated beyond tIn
%3) 'spline' removes any sharp edges if vIn is undersampled by tIn
%4) 'same' ensures that we do not keep data before the first value of tIn
%   and last value of tConv. This only works if filterIn is centered!
%   This also allows to use tConv as the x-values for back-interpolation to the
%   non-uniform grid.
vOut = conv(interp1(tIn,vIn,tConv,'spline','extrap'),filterIn,'same')/n;
vOut = interp1(tConv,vOut,tIn,'linear','extrap'); %back-interpolate to non-uniform grid