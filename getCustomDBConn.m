function conn = getCustomDBConn()
% GETCUSTOMDBCONN causes dbConn() to return a custom database connection. 
% Modify this file if do not wish to use the default database connection
% that ships with CukLabScripts. This may be useful if you do not have MS
% Access ODBC connection capabilities, such as on a macOS.
%
% To use getCustomDBConn(), place this function in your MATLAB user
% folder, returned by userpath() and typically located under 
% %USERPROFILE%\Documents\MATLAB. Modify the conn variable to return a
% custom database connection.
%
% conn = getCustomDBConn()
%   Returns a custom database connection.
%
% See Also: database, dbConn, userpath

% Modify path below for your custom repository path
db_path = fullfile(filesep,'Users','suryansh','Work','CUB','TC research group','Data','New Data','Raman.accdb');
url = ['jdbc:ucanaccess://' db_path];
conn = database('','','','net.ucanaccess.jdbc.UcanaccessDriver',url);
