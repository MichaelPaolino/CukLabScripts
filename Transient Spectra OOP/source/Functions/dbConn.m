function conn = dbConn()
% DBCONN returns a database connection to the Cuk Lab database.
% The default below can be overriden in getCustomDBConn.m
%
% conn = dbConn()
%   Returns a database connection object with optional overwrite.
%
    
    % Check if user has a custom repository path
    if isfile(fullfile(userpath(),'getCustomDBConn.m'))
        conn = getCustomDBConn();
    else
        conn = database('TRRaman DB','','');
    end
end