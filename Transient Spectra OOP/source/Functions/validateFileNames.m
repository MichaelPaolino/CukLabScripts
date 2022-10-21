function namesOut = validateFileNames(namesIn)
    % VALIDATEFILENAMES ensures that the input names are all unique and do not
    % contain invalid file name characters. The following conversions are made 
    % in addition to MATLABs name validation:
    % ' ' -> '_'
    % '.' -> 'dot'
    % '%' -> 'per'
    % '\' or '/' -> s
    % ';' -> col
    % ''' -> ''
    % ',' -> ''
    %
    % Note: names are trimmed at 62 chars due to matlab name validation
    %
    % namesOut = VALIDATEFILENAMES(namesIn)
    %   Converts elements of cell array namesIn into unique and valid names and
    %   outputs them as namesOut.
    
    replaceRules = {'.','dot';...
                    ' ','_';...
                    '%','per';...
                    '/','fs';...
                    '\','bs';...
                    '''','';...
                    ',','';...
                    newline,'_';...
                    char(13),'_'};
    
    % function specific repalce rules
    for ii = 1:size(replaceRules,1)
        namesIn = strrep(namesIn,replaceRules{ii,1},replaceRules{ii,2});
    end
    
    % todo: splice file names into 62 char subsets for matlab name validation
    
    % Matlab name validation replace rules (max name length is 42 chars)
    namesIn = matlab.lang.makeValidName(strcat('x', namesIn));   %prefix char to perserve 1st numeric char
    namesIn = cellfun(@(c) c(2:end),namesIn,'UniformOutput',false); %remove prefix
    
    % Generate unique names
    namesOut = matlab.lang.makeUniqueStrings(namesIn);

end