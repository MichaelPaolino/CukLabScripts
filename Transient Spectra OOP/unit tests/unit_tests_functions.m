%Unit tests for functions
%% findGridIntercepts.m
x = [0:0.01:2.99 3.000001];
y = sin(2*pi*x);
xInt = findGridIntercepts(x,y,0);
soln = [0, 0.5, 1, 1.5, 2, 2.5, 3];
assert(all((xInt(:)-soln(:)) < 1e-6), 'All intercepts were not found for sin');

y = [zeros(1,length(x)-1) 1];
xInt = findGridIntercepts(x,y,0);

%% nearestVal.m
x = 1:5;
[vals,ind] = nearestVal(x,1.1);
assert(vals==1 && ind==1,'failed to find nearest val for a vector with 1 target val');

[vals,ind] = nearestVal(x,[1.1, 2.1]);
assert(all(vals(:)==[1;2]) && all((ind)==[1;2]),...
    'failed to find nearest val for a vector with 2 target vals');

[vals,ind] = nearestVal(repmat(x(:),1,2),[1.1, 2.1]);
assert(all(vals==[1,1;2,2],'all') && all(ind==[1,1;2,2],'all'),...
    'failed to find nearest val for a matrix with 2 target vals');

[vals,ind] = nearestVal(repmat(x(:),1,2,2),1.2);
assert(all(size(vals)==[1,2,2]) && all(size(ind)==[1,2,2]),...
    'failed to return correct vals and ind size for a 3d array with 1 target val');
assert(all(vals(:)==ones(4,1)) && all(ind(:)==ones(4,1)),...
    'failed to find nearest val for a 3d array with 1 target val');

[vals,ind] = nearestVal(repmat(x(:),1,2,2),[1.2, 2.2]);
assert(all(size(vals)==[2,2,2]) && all(size(ind)==[2,2,2]),...
    'failed to return correct vals and ind size for a 3d array with 2 target val');
assert(all(vals(:)==[1,2,1,2,1,2,1,2]') && all(ind(:)==[1,2,1,2,1,2,1,2]'),...
    'failed to find nearest val for a 3d array with 2 target val');

[vals,ind] = nearestVal(repmat(x(:),1,2,2),[1.2, 1.3],'unique',false,'trim',false);
assert(all(size(vals)==[2,2,2]) && all(size(ind)==[2,2,2]),...
    'failed to return correct vals and ind size for a 3d array with 2 target val with unique set to false');
assert(all(vals(:)==[1,1,1,1,1,1,1,1]') && all(ind(:)==[1,1,1,1,1,1,1,1]'),...
    'failed to find nearest val for a 3d array with 2 target val with unique set to false');

[vals,ind] = nearestVal(repmat(x(:),1,2,2),[1.2, 1.3]);
assert(all(size(vals)==[1,2,2]) && all(size(ind)==[1,2,2]),...
    'failed to return correct vals and ind size for a 3d array with 2 non-unique target val');
assert(all(squeeze(vals(1,:,:))==[1,1;1,1],'all') && all(squeeze(vals(1,:,:))==[1,1;1,1],'all'),...
    'failed to find nearest val for a 3d array with 2 non-unique target val');

[vals,ind] = nearestVal(zeros(1,1,3),0);
assert(all(size(vals)==[1,1,3]),'failed to preserve size of array')

[vals,ind] = nearestVal(zeros(1,1,3),[0,0],'trim',false);
assert(all(size(vals)==[2,1,3]),'failed to preserve size of array with two targetVals')

[vals,ind] = nearestVal([0,1],0);
assert(all(size(vals)==[1,1]),'failed to do a row vector search');

[vals,ind] = nearestVal([0,1],0,'forceMat',true);
assert(all(size(vals)==[1,2]),'failed to do a column search on a row vector');

[vals,ind] = nearestVal([1:10; 6:15]',1:0.1:10,'threshold',0.1);
assert(all(vals(:,1)==(1:10)') && all(isnan(vals(1:5,2))) && all(vals(6:10,2)==(6:10)'),...
       'failed to return correct vals with dense targetVals');
assert(all(ind(:,1)==(1:10)') && all(isnan(ind(1:5,2))) && all(ind(6:10,2)==(1:5)'),...
       'failed to return correct ind with dense targetVals');
   
%% genFKStruct.m
fStruct = genFKStruct('UserID','Users','ID','FirstName',1);
tCell = {'UserID','Users','ID','WHERE FirstName = 1'};
assert(isequal(struct2cell(fStruct),tCell'), 'Failed to generate correct fk struct');

fStruct = genFKStruct('UserID','Users','ID','FirstName','Ilya');
tCell = {'UserID','Users','ID','WHERE FirstName = ''Ilya'''};
assert(isequal(struct2cell(fStruct),tCell'), 'Failed to generate correct fk struct');

fStruct = genFKStruct('UserID','Users','ID','FirstName',{'Ilya','Michael'});
tCell = {'UserID','Users','ID','WHERE FirstName = ''Ilya''';...
         'UserID','Users','ID','WHERE FirstName = ''Michael'''};
assert(isequal(struct2cell(fStruct),tCell'), 'Failed to generate correct fk struct');

fStruct = genFKStruct('UserID','Users','ID','FirstName',{'Ilya',1});
tCell = {'UserID','Users','ID','WHERE FirstName = ''Ilya''';...
         'UserID','Users','ID','WHERE FirstName = 1'};
assert(isequal(struct2cell(fStruct),tCell'), 'Failed to generate correct fk struct');

fStruct = genFKStruct({'UserID','SampleID'},...
                      {'Users','Samples'},...
                       'ID',...
                      {'FirstName','ShortName'},...
                      {'Ilya','Sample S'});
tCell = {'UserID','Users','ID','WHERE FirstName = ''Ilya''';...
         'SampleID','Samples','ID','WHERE ShortName = ''Sample S'''};
assert(isequal(struct2cell(fStruct),tCell'), 'Failed to generate correct fk struct');

%% splitPath.m
cellOut = {'C:','Users','User','Documents','MATLAB','CukLabScripts','Transient Spectra OOP'};

pathOut = splitPath('C:\Users\User\Documents\MATLAB\CukLabScripts\Transient Spectra OOP');
assert(isequal(cellOut(:),pathOut(:)),'Failed to split path');

pathOut = splitPath('C:/Users/User/Documents/MATLAB/CukLabScripts/Transient Spectra OOP');
assert(isequal(cellOut(:),pathOut(:)),'Failed to split path');

pathOut = splitPath({'C:/Users/User/Documents/MATLAB/CukLabScripts/Transient Spectra OOP'});
assert(isequal(cellOut(:),pathOut(:)),'Failed to split path');

% not a realistic use case, just to see if it works...
pathOut = splitPath({'C:/Users\User/Documents\MATLAB/CukLabScripts\Transient Spectra OOP'});
assert(isequal(cellOut(:),pathOut(:)),'Failed to split path');

% check relative vs absolute paths
% Absolute path on MacOS
cellOut = {'/','Users','User','Documents','MATLAB','CukLabScripts','Transient Spectra OOP'};
pathOut = splitPath('/Users/User/Documents/MATLAB/CukLabScripts/Transient Spectra OOP');
assert(isequal(cellOut(:),pathOut(:)),'Failed to split path');

% Relative path on Windows
cellOut = {'Users','User','Documents','MATLAB','CukLabScripts','Transient Spectra OOP'};
pathOut = splitPath('Users\User\Documents\MATLAB\CukLabScripts\Transient Spectra OOP');
assert(isequal(cellOut(:),pathOut(:)),'Failed to split path');

%% comparePaths.m
path1 = fullfile(repoPath(),'test1','test2','test3.mat');
path2 = fullfile(repoPath(),'test4','test5');

[sameP, relP1, relP2] = comparePaths(path1, path2);
assert(strcmp(sameP,repoPath()) && strcmp(relP1, fullfile('test1','test2','test3.mat')) && strcmp(relP2, fullfile('test4','test5')),...
    'Failed to compare paths');

[sameP, relP1, relP2] = comparePaths(repoPath(), path2);
assert(strcmp(sameP,repoPath) && isempty(relP1) && strcmp(relP2, fullfile('test4','test5')),...
    'Failed to compare paths');

[sameP, relP1, relP2] = comparePaths(path1, repoPath());
assert(strcmp(sameP,repoPath()) && strcmp(relP1, fullfile('test1','test2','test3.mat')) && isempty(relP2),...
    'Failed to compare paths');

[sameP, relP1, relP2] = comparePaths(path1, path2,'/');
assert(strcmp(sameP,repoPath('/')) && strcmp(relP1, strjoin({'test1','test2','test3.mat'},'/')) && strcmp(relP2, strjoin({'test4','test5'},'/')),...
    'Failed to compare paths');

% test with linux root '/'
[sameP, relP1, relP2] = comparePaths('/asdf/asdf/p', '/asdf/bsdf/p','/');
assert(strcmp(sameP,'/asdf') & strcmp(relP1,'asdf/p') & strcmp(relP2,'bsdf/p'),'Failed to compare paths');

%% isAbsolutePath.m
tf = isAbsolutePath({'/asdf\asdf\asdf','asdf\asdf\asdf','C:\asdf\asdf\asdf'});
assert(all(tf==[1,0,1]),'Absolute Path Failed');

tf = isAbsolutePath({'/asdf\asdf\asdf','asdf\asdf\asdf','C:\asdf\asdf\asdf'}');
assert(all(tf==[1,0,1]'),'Absolute Path Failed');

tf = isAbsolutePath('/asdf\asdf\asdf');
assert(tf,'Absolute Path Failed');

%% convertPaths.m
% todo: write assert statements...
pathOut = convertPaths('/Users/Raman/Documents/LabVIEW Code/RamanProgram/Raman-OOP/Scripts/MATLAB/MATLAB Data API');
pathOut = convertPaths('/Users/Raman/Documents/LabVIEW Code/RamanProgram/Raman-OOP/Scripts/MATLAB/MATLAB Data API');
pathOut = convertPaths(repmat({'/Users/Raman/Documents/LabVIEW Code/RamanProgram/Raman-OOP/Scripts/MATLAB/MATLAB Data API'},3,2));