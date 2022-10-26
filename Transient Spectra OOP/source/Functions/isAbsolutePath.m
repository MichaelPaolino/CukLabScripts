function tf = isAbsolutePath(pathIn)
% ISABSOLUTEPATH returns a logical array the same size as the input path
% array that indicates whether the input paths are absolute or relative
% paths. This function works by checking for the drive letter, 'C:', or 
% root, '/', at the start of the path.

% Validate inputs and convert to cell array
if ischar(pathIn)
    pathIn = {pathIn};
elseif iscell(pathIn)
    assert(all(cellfun('isclass',pathIn,'char'),'all'),'Cell input for pathIn expected to contain char arrays.');
else
    error('Invalid input type for pathIn. Expected char or cell array of char');
end

% Get size of pathIn for easy vectorization
sz = size(pathIn);
pNumel = numel(pathIn);
pathIn = pathIn(:);

tf = false(pNumel,1);

% Determine if path is an absolute path
for pInd = 1:pNumel
    if numel(pathIn{pInd})>=1 && strcmp(pathIn{pInd}(1),'/')
        tf(pInd) = true;
    elseif numel(pathIn{pInd})>=2 && strcmp(pathIn{pInd}(2),':')
        tf(pInd) = true;
    end
end

% reshape logical array to be the same size as input
tf = reshape(tf,sz);