function fkStruct = genFKStruct(fkCol, pkTable, pkCol, varargin)
% GENFKSTRUCT generates the foreign key input structure for
% transientSpectra dbCommit. Use this structure to set foreign key values
% by the display column of the primary key table.
%
% Note: This function does not check for valid input tables or columns 
% against the database
%
% Definitions:
% Foreign Key: column that points to a primary key column of another table
% Primary Key: column that contains a unique indentifier for a row
%
% Foreign Key Table: Table name that contains the foreign key column.
% Primary Key Table: Table name that the foreign key points to.
% Primary Key Display Column: A column that contains user-identifiable
%   information in the primary key table
%
% fkStruct = genFKStruct(fkCol, pkTable, pkCol, pkDispCol, pkDispVal);
%   Generates a struct with fields fkTable, fkCol, pkTable, pkCol,
%   pkDispCol, and pkDispVal. Use the struct as an input for dbCommit.
%   Inputs can be cell arrays of char or char array. When inputs are cell
%   arrays, this function returns a struct array the same size as the input
%   cell array. Any singleton inputs are expanded to match the size of the
%   largest cell array. Mixed cell/char array inputs are allowed.
%
% todo: add WHERE clause option
%
% See Also: TRANSIENTSPECTRA, DBCOMMIT
    
    % parse varargin for user input
    switch nargin
        case 4  % User supplied a WHERE statement
            whrStmt = varargin{1};
            % Ensure that the WHERE statement contains 'WHERE'
            if iscell(whrStmt)
               for ii = 1:numel(whrStmt)
                   assert(contains(whrStmt{ii},'WHERE'),'WHERE statement must contain WHERE clause.');
               end
            else
               assert(contains(whrStmt,'WHERE'),'WHERE statement must contain WHERE clause.');
            end
        case 5  % User supplied a value column and its value
            % Convert value element into cell if not already a cell array
            if ~iscell(varargin{2})
                varargin{2} = {varargin{2}};
            end
            
            % Add '...' to strings and convert double to strings
            for ii = 1:numel(varargin{2})
               if isnumeric(varargin{2}{ii})
                   varargin{2}{ii} = num2str(varargin{2}{ii});
               else
                   varargin{2}{ii} = ['''' char(varargin{2}{ii}) ''''];
               end
            end
            
            % Generate WHERE statement
            whrStmt = strcat({'WHERE '}, varargin{1}, {' = '}, varargin{2});
        otherwise
            error('genFKStruct requires 4 or 5 arguments.');
    end

    % convert input args into a single cell array for looping
    argsIn = {fkCol, pkTable, pkCol, whrStmt};
    argsNumel = zeros(size(argsIn));   %this will hold the size of each argument
    
    % loop over all arguments, convert to cell array, get and store size
    argsIn = struct2cell(ensureCellVals(cell2struct(argsIn,{'v1','v2','v3','v4'},2)));    %ensure cell array
        
    for ii = 1:numel(argsIn)
        argsIn{ii} = argsIn{ii}(:); %ensure column cell array
        argsNumel(ii) = numel(argsIn{ii}); %number of elements of each cell array
    end
    
    % max number of rows in the struct to output
    maxNumel = max(argsNumel);
    
    % check to make sure sizes are valid
    assert(all(any([argsNumel(:)==maxNumel, argsNumel(:)==1],2),1),...
           ['The dims of the input arguments cannot be expanded to a matching size. ',...
           'Make sure the input argument cell arrays have either 1 element or the same maximum number of elements']);
       
    % make sure all inputs are the same size as the max
    for ii = 1:numel(argsIn)
        argsIn{ii} = explicitExpand(argsIn{ii},[maxNumel,1]);
    end
    
    % build output structure
    fkStruct = struct('fkCol',argsIn{1},...
                      'pkTable',argsIn{2},...
                      'pkCol',argsIn{3},...
                      'WHERE',argsIn{4});
