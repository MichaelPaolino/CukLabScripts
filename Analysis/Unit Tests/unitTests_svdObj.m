%% Load phosphate metadata
conn = database('Raman DB','','');
dataList = conn.fetch(strjoin(...
    ["SELECT ID, FileName, FilePath, ShortName",...
     "FROM TRPhononRemoved",...
     "WHERE ID IN (63, 7) OR ID BETWEEN 66 AND 74",...
     "ORDER BY ID;"],' '));
 
conn.close();

%% Load phosphate data into object
% Extract pH values and sort in ascending order
c = regexp(dataList.ShortName,'(?<=pH )\d+\.?\d*','match');
[pH, idx] = sort(cellfun(@(c) str2double(c), vertcat(c{:})));
dataList = dataList(idx,:);

% Load data
myTR = wlTR(fullfile(repoPath(), dataList.FilePath, dataList.FileName),'loadType','cListFile','shortName',dataList.ShortName);

%% Clean data and calculate SVD in SVD object
% should not return any errors
[s,l,t] = getNumeric(myTR.trim('wavelengths',[375,700]));  %return "cleaned" data from myTR object array
trSVD = svdObj(s);   % perform the SVD analysis
trSVD = svdObj(s,3);   % perform the SVD analysis with 3 components
trSVD = trSVD.setC(2); % set componment number back to 2

%% test flip basis components
% should not return any errors
trSVD = trSVD.resetC(); %resets flips
trSVD = trSVD.permuteC(1:trSVD.nSets,[2,1]); %swaps absorptive with emissive
trSVD = trSVD.flipC(1,[1,7,9,11],true);
trSVD = trSVD.flipC(2,[1,2,3,4,6,7,8,9,10,11],true);

%% Test plotting basis components
figure; trSVD.previewr(2);
figure; trSVD.previewr(4,'xVals',l,'yVals',t,'xLabel','nm','yLabel','ps');
figure; trSVD.preview('xVals',l,'yVals',t,'xLabel','nm','yLabel','ps');
figure; trSVD.previewc('xVals',l,'yVals',t,'xLabel','nm','yLabel','ps');

%% Test plotting constrained basis components after chaning nC
trSVD = trSVD.setC(3);
figure; trSVD.previewc('xVals',l,'yVals',t,'xLabel','nm','yLabel','ps');
figure; trSVD.preview('xVals',l,'yVals',t,'xLabel','nm','yLabel','ps');
trSVD = trSVD.setC(2);
figure; trSVD.previewc('xVals',l,'yVals',t,'xLabel','nm','yLabel','ps');

%% Test error and warning generation related to component number
% trSVD = trSVD.setC(2);
% trSVD = trSVD.permuteC(1,3:-1:1);   %throws error

trSVD = trSVD.setC(3);
trSVD = trSVD.permuteC(1,3:-1:1);   
trSVD = trSVD.setC(2);              %throws warning

%% Test constrain and rotation analysis
% shoud not return errors
trSVD = trSVD.setC(3);
trSVD = trSVD.constrain(trSVD.U(:,:,1));
trSVD.rotAnalysis();

%% Extract SVD calculation results into workspace
i = 1;
U = trSVD.U(:,:,i);
V = trSVD.V(:,:,i);
S = trSVD.S(:,:,i);
Ur = trSVD.Ur(:,1:2,i);
Vr = trSVD.Vr(:,1:2,i);
Sr = trSVD.Sr(1:2,1:2,i);
T = trSVD.T(:,:,i);
Xr = trSVD.Xr(:,:,i);
X = trSVD.X(:,:,i);

%% check that display transformations work on constrained spectra
Mr = Ur*Sr*Vr';
M = U*S*V';
assert(all((M-Mr)==0,'all'),'Failed to reconstruct M.');

Mr = Ur*Xr*inv(Xr)*Sr*Vr';
M = U*X*inv(X)*S*V';
assert(all((M-Mr)==0,'all'),'Failed to reconstruct M.');

%% Use this to display failing M and Mr
figure; 
subplot(1,3,1);
    imagesc(s(:,:,1));
subplot(1,3,2);
    imagesc(Mr);
subplot(1,3,3);
    imagesc(M);