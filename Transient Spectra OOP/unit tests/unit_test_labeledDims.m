% Unit tests for nameRule class

%% test constructors
% default constructor
myDims = labeledDims();

% constructor with double data
myDims = labeledDims(1:5);

% constructor with double sublcass
myDims = labeledDims(doubleWithUnits(1:5));

%% All tests pass
disp('All tests passed!');

%% functions specific to test script