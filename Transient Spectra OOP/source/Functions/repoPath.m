function pathOut = repoPath(varargin)
% REPOPATH returns the CukLab google drive data repository path.
%
% pathOut = repoPath();
%   Returns char array pathOut for the google drive using a
%   platform-specific separator
%
% pathOut = repoPath(fs);
%   Returns char array pathOut for the google drive using char fs as the
%   file separator

    switch nargin
        case 0
            fs = filesep();
        case 1
            assert(ischar(varargin{1}),'File separator expected to be of type char');
            fs = varargin{1};
        otherwise
        error('repoPath() expects 0 or 1 inputs')
    end
    
    % Check if user has a custom repository path
    if isfile(fullfile(userpath(),'getCustomRepoPath.m'))
        tmp = splitPath(getCustomRepoPath());
        if ~isempty(tmp) && strcmp(tmp{1},'/')
            tmp(1) = [];
            pathOut = ['/',strjoin(tmp,fs)];
        else
            pathOut = strjoin(tmp,fs);
        end
    else
        pathOut = strjoin({'G:','.shortcut-targets-by-id','1O21NYESEEKxBx8k6l21vr3Feo26hcErZ'},fs);
    end
end