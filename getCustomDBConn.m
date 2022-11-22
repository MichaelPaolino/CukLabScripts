function cRepoPath = getCustomRepoPath()
% GETCUSTOMREPOPATH causes repoPath() to return a custom repository path. 
% Modify this file if do not wish to use the default repository path that 
% ships with CukLabScripts. This may be useful if you have a local copy of 
% the repository that is set up to periodically synchronize with the full 
% data repository.
%
% To use getCustomRepoPath(), place this function in your MATLAB user
% folder, returned by userpath() and typically located under 
% %USERPROFILE%\Documents\MATLAB. Modify the cRepoPath to correspond to
% your custom repository path
%
% cRepoPath = getCustomRepoPath()
%   Returns a custom repository path.
%
% See Also: repoPath, userpath

% Modify path below for your custom repository path
cRepoPath = 'G:\.shortcut-targets-by-id\1O21NYESEEKxBx8k6l21vr3Feo26hcErZ';