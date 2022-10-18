function namesOut = validateFileNames(namesIn)
    % VALIDATEFILENAMES ensures that the input names are all unique and do not
    % contain invalid file name characters. The following conversions are made 
    % in addition to MATLABs name validation:
    % ' ' -> '_'
    % '.' -> 'p'
    % '%' -> 'per'
    % '\' or '/' -> s
    % ';' -> col
    % ''' -> ''
    %
    % namesOut = VALIDATEFILENAMES(namesIn)
    %   Converts elements of cell array namesIn into unique and valid names and
    %   outputs them as namesOut.
    
    replaceRules = {'.','p';...
                    ' ','_';...
                    '%','per';...
                    '/','s';...
                    '\','s';...
                    '''','';...
                    newline,'_';...
                    char(13),'_'};
    
    % function specific repalce rules
    for ii = 1:size(replaceRules,1)
        namesIn = strrep(namesIn,replaceRules{ii,1},replaceRules{ii,2});
    end

    % Matlab name validation replace rules
    % namesIn = matlab.lang.makeValidName(namesIn);

    % Generate unique names
    namesOut = matlab.lang.makeUniqueStrings(namesIn);

end