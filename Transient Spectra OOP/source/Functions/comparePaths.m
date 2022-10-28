function [sameP, relP1, relP2] = comparePaths(P1, P2, varargin)
% COMPAREPATHS compares two paths and returns their common path and
% relative paths away from the common path. This function works with paths
% generated with any platform-dependent path separator. Paths that start
% with '/' are treated as absolute paths where '/' is the root symbol on
% MacOS and Linux platforms, equivalent to e.g. 'C:' on Windows platforms. 
% Both root, '/', and the drive letter, e.g. 'C:', are treated as path
% elements by this function.
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
% See Also: fileparts, filesep, fullfile, splitPath, split, strjoin

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
    if ~isSame1(ii)
        break;
    end
end

% Use strjoin and fs defined above to rebuild common and relative paths
% Take subset of the split paths
sameP = p1Split(isSame1);
relP1 = p1Split(~isSame1);
relP2 = p2Split(~isSame2);

% use strjoin to rebuild paths, but include special case to handle root in (MacOS and Linux)
if ~isempty(sameP) && strcmp(sameP(1),'/')
    sameP(1) = [];
    sameP = ['/' strjoin(sameP,fs)];
else
    sameP = strjoin(sameP,fs);
end

if ~isempty(relP1) && strcmp(relP1(1),'/')
    relP1(1) = [];
    relP1 = ['/' strjoin(relP1,fs)];
else
    relP1 = strjoin(relP1,fs);
end

if ~isempty(relP2) && strcmp(relP2(1),'/')
    relP2(1) = [];
    relP2 = ['/' strjoin(relP2,fs)];
else
    relP2 = strjoin(relP2,fs);
end