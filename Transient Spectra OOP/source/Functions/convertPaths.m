function pathsOut = convertPaths(pathsIn)
% CONVERTPATHS converts input paths file separators into platform-dependent
% file separators. This function accepts both char and cell array inputs.
% For cell array inputs, the output is the same size as the input.
%
% pathsOut = convertPaths(pathsIn)
%   Converts file separators in char or cell array pathsIn into
%   platform-dependent file separators, pathsOut.
%
% See Also: splitPath, filesep, fullfile, splitPath

    % Convert all inputs to cell
    isPathCell = iscell(pathsIn);
    if ~isPathCell
        pathsIn = {pathsIn};
    end
    
    % Get input size
    sz = size(pathsIn);
    pathsNumel = numel(pathsIn);
    pathsIn = pathsIn(:);
    pathsOut = cell(pathsNumel,1);
        
    % Convert all paths to platform dependent file sep
    for i = 1:pathsNumel
        tmp = splitPath(pathsIn{i});
        pathsOut{i} = fullfile(tmp{:});
    end
    
    % Convert pathsOut to input size
    pathsOut = reshape(pathsOut,sz);
    
    % Convert to non-cell if not cell
    if ~isPathCell
        pathsOut = pathsOut{1};
    end
end