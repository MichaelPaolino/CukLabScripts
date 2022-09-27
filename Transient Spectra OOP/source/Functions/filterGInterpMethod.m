function [vOut, mInterp, mExtrap] = filterGInterpMethod(vIn,dInterp,dExtrap)
% isGInterpMethod is a helper function for griddedInterpolant wrappers.
% This function takes a set of input varargin and psitionally filters them 
% for valid interpolation and extrapolation method keywords.
%
% Valid interpolation methods are:
%   'linear', 'nearest', 'next', 'previous', 'pchip', 'cubic', 'spline', 
%   and'makima'
%
% Valid extrapolation methods are:
%   'linear', 'nearest', 'next', 'previous', 'pchip', 'cubic', 'spline', 
%   'makima', and 'none'
%
% [vOut, mInterp, mExtrap] = filterGInterpMethod(vIn,dInterp,dExtrap)
%   Filters cell array of inputs vIn into non-associated inputs vOut and
%   interpolation method mInterp and extrapolation method mExtrap.
%   vIn is filtered positionally by looking for the methods at the end of
%   the cell array. This function sets the interpolation method first, and
%   the extrapolation method second. If interpolation methods or 
%   extrapolation methods are not found, the default dInterp and dExtrap 
%   methods are returned for mInterp and mExtrap instead.
%   
%   Example:
%       vIn = {'wavelengths',[400:5:500],'spline'}
%       dInterp = 'linear'
%       dExtrap = 'linear'
%   Returns:
%       vOut = {'wavelengths',[400:5:500]}
%       mInterp = 'spline'
%       mExtrap = 'linear'
% 
% See Also: GRIDDEDINTERPOLANT

% Define filter keywords
mInterpFilt = {'linear', 'nearest', 'next', 'previous', 'pchip', 'cubic', 'spline', 'makima'};
mExtrapFilt = [mInterpFilt, {'none'}];

% set default values for outputs
vOut = vIn;
mInterp = dInterp;
mExtrap = dExtrap;

% Filter for first keyword
for vInd = 1:numel(vOut)
    % Check if input is a char
    if ischar(vOut{vInd})
        % Check if input is an interpolation method
        if any(strcmp(vOut{vInd}, mInterpFilt))
            % put found index into mInterp and remove from vIn
            mInterp = vOut{vInd};
            vOut(vInd) = [];
            break;
        end
    end
end

% If the first keyword was found, look for the second
if numel(vIn) ~= numel(vOut)
    for vInd = 1:numel(vOut)
        % Check if input is a char
        if ischar(vOut{vInd})
            % Check if input is an interpolation method
            if any(strcmp(vOut{vInd}, mExtrapFilt))
                % put found index into mInterp and remove from vIn
                mExtrap = vOut{vInd};
                vOut(vInd) = [];
                break;
            end
        end
    end
end
