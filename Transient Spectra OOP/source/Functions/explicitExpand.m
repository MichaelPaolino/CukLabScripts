function arrayOut = explicitExpand(arrayIn, arraySize)
% EXPLICITEXPAND explicitly expands singleton dims of an input array to a
% target array size. This function checks that the size of the input array 
% matches the target size, except for singleton dims, beofre attempting the
% expansion. If the size does not match, this function throws an error.
%
% Use this function for vectorization when MATLAB's implicit expansion is
% not supported for your specific use case.
%
% arrayOut = EXPLICITEXPAND(arrayIn, arraySize)
%   For example:
%   If arrayIn is [5] and arraySize is [5,5], this function returns a 5x5 
%       matrix whose elements are 5.
%   If arrayIn is 1:5 and arraySize is [5,5], this function returns a 5x5
%       matrix whose rows are 1:5.
%   If arrayIn is 1:3 and arraySize is [3,6], this function returns an
%       error because size(arrayIn) is 1x3 and arraySize corresponds to 3x6
%
% See Also: SIZE, REPMAT, RESHAPE, PERMUTE

% Get the size vector of arrayIn, padded with singleton dims to match the length of arraySize
tmpSize = size(arrayIn);
inSize = ones(numel(arraySize),1);
inSize(1:numel(tmpSize)) = tmpSize;

% ensure that inSize matches arraySize, except for singleton dims
assert(all(or(inSize(:) == arraySize(:),inSize(:) == 1)),'Non-singleton array dims of arrayIn must match arraySize.');

% Use repmat to explicitly expand arrayIn to be of size arraySize
arrayOut = repmat(arrayIn,(arraySize(:)./inSize(:))');

end