function [sameP, relP1, relP2] = comparePaths(P1, P2, varargin)
% COMAPREPATHS compares two paths and returns their common path and
% relative paths away from the common path. This function works with paths
% generated with any platform-dependent path separator.
%
% [sameP, relP1, relP2] = comparePaths(P1, P2)
%   Returns the common path, sameP char array, and relative difference
%   paths, relP2 and relP2 char arrays, between two input paths, P1 and P2
%   char arrays. Outputs are given with a platform-specific path separator.
%
% [sameP, relP1, relP2] = comparePaths(P1, P2, fs)
%   Same call as above except the outputs are given with fs as the path
%   separator.
%
% See Also: 

% Validate inputs are char arrays
assert(ischar(P1) && ischar(P2), 'Input paths must be char arrays.');

% Set the file separator
switch nargin
    case 2
        fs = filesep();
    case 3
        assert(ischar(varargin{1}),'File separator expected to be of type char');
        fs = varargin{1};
    otherwise
        error('comparePaths expects 2 or 3 inputs')
end

% Split paths into cell arrays (makes platform-independent)
p1Split = splitPath(P1);
p2Split = splitPath(P2);

% Find common path inds between two split paths by using logical inds
numelShortest = min([numel(p1Split),numel(p2Split)]);
isSame1 = false(numel(p1Split),1);
isSame2 = false(numel(p2Split),1);

for ii = 1:numelShortest
    isSame1(ii) = strcmp(p1Split{ii},p2Split{ii});
    isSame2(ii) = isSame1(ii);
end

% Use strjoin and fs defined above to rebuild common and relative paths
sameP = strjoin(p1Split(isSame1),fs);
relP1 = strjoin(p1Split(~isSame1),fs);
relP2 = strjoin(p2Split(~isSame2),fs);