function objOut = validateObjData(dataIn, objClass)
% VALIDATEOBJDATA validates and returns and object from a .mat load struct.
% This function returns an error if the loaded data does not contain the
% target object type, returns a warning if the data contains more than one
% variable, and returns a warning if the data contains an object array. In
% these cases, the user should load data directly using load() instead of
% the object constructor import method.
% 
% objOut = validateObjData(dataIn, objClass)
%   Returns the 1st object element of class objClass from .mat load struct
%   or cell dataIn.

% Convert dataIn into cell
if ~iscell(dataIn)
    if isstruct(dataIn)
       dataIn = struct2cell(dataIn); 
    else
       error('Expected dataIn to be a cell array or struct.'); 
    end
end

% Verify class of dataIn elements
isObjClass = false(1,numel(dataIn));
for ii = 1:numel(dataIn)
    isObjClass = isa(dataIn{ii},objClass);
end

% Make sure obj of objClass exists in .mat file
assert(any(isObjClass),'Attempted to import a .mat file that did not contain an object of the target class.');

% Keep only the 1st element that is of type objClass
dataIn = dataIn(isObjClass); %subset that is type objClass
dataIn = dataIn{1}; %1st element that is objClass

% Warn user there was more than one variable
if numel(isObjClass) > 1
   warning('Imported a .mat file with more than one variable. If you need all variables call load(...) directly.'); 
end

% Warn user there was more than one objClass
if sum(isObjClass) > 1
    warning('Imported a .mat file that contains more than one object of the target class. Keeping only the first variable encountered. If you need all variables call load(...) directly.');
end

% warn user that there is an object array
if numel(dataIn) > 1
    warning('Imported a .mat file that contains an object array. Keeping only the 1st element. If you need all elements call load(...) directly.');
end

% keep 1st element of object array and return
objOut = dataIn(1); 