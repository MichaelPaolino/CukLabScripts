function sz = sizePadded(s,dims)
% SIZEPADDED is a wrapper for size that returns size of s for specific 
% dims. Use this function in older versions of MATLAB where size does not 
% support sizes of multiple dims.
%
% sz = SIZEPADDED(s, dims)
%   Returns size vector sz which contains the size of s for specified
%   vector dims. sz is padded with 1 for any value of dims which is larger
%   than the number of dims in s.
%
% See Also: SIZE

dims = dims(:);
szP = ones(1,max(dims));
sz = size(s);
szP(1:numel(sz)) = sz;
sz = szP(dims);