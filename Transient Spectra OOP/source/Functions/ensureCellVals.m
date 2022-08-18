function structOut = ensureCellVals(structIn)
% ENSURECELLVALS converts all non-cell elements of the input struct to 
% cells. This means if the element is already a cell, a nested cell will
% not be created. Meanwhile, if the element is not a cell, its contents
% will be wrapped in a 1x1 cell array.
%
% structOut = ensureCellVals(structIn)
%   Converts all non-cell elements of structIn into cells

% get field names
fn = fieldnames(structIn);

% if input is a struct array, convert into vector
structSize = size(structIn);
structNumel = numel(structIn);
structIn = structIn(:);

% loop over field names and array elements
for fInd = 1:numel(fn)
    for aInd = 1:structNumel
        % if the element is not a cell convert it to a cell
        if ~iscell(structIn(aInd).(fn{fInd}))
            structIn(aInd).(fn{fInd}) = {structIn(aInd).(fn{fInd})};
        end
    end
end 

% restore original struct size
structOut = reshape(structIn,structSize);

end