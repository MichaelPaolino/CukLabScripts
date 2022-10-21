function fileParts = splitPath(pathIn)
% SPLITPATH splits a file path into component parts at platform-dependent
% file sperators. This method searches for all platform-dependent
% separators.
%
% fileParts = splitPath(pathIn)
%   Returns a cell array of char, fileParts, which are the path components
%   char array pathIn. 
%
% See Also: FILEPARTS, FILESEP, FULLFILE, SPLIT

fileSepChars = {'/','\'};

% Validate inputs and convert to cell array
if ischar(pathIn)
    fileParts = {pathIn};
elseif iscell(pathIn)
    fileParts = pathIn(:);
    assert(all(cellfun('isclass',fileParts,'char')),'Cell input for pathIn expected to contain char arrays.');
else
    error('Invalid input type for pathIn. Expected char or cell array of char');
end

% Loop over file seperator characters
for ii = 1:numel(fileSepChars)
    
    % Loop over elements of the fileParts cell array
    for jj = 1:numel(fileParts)
        % Use split to split elements of the fileParts cell array by the
        % file seperator. This replaces the current element of fileParts
        % with a cell array of variable length.
        fileParts{jj} = split(fileParts{jj},fileSepChars{ii});
    end
    
    % Flatten nested cell arrays
    fileParts = vertcat(fileParts{:});
end


    