e a%% Load phosphate metadata
conn = database('Raman DB','','');

% 0.1% phosphate data in Nat. Mat. paper
dataList = conn.fetch(strjoin(...
    ["SELECT ID, FileName, FilePath, ShortName",...
     "FROM TRPhononRemoved",...
     "WHERE ID IN (63, 24) OR ID BETWEEN 66 AND 74",...
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
[s,l,t] = getNumeric(myTR.trim('wavelengths',[375,700]));  %return "cleaned" data from myTR object array
trSVD = svdObj(s);   % perform the SVD analysis with default 2 components

%% Plot results (raw)
% Unflipped (raw) SVD preview
figure; previewr(trSVD,2,'xVals',l,'yVals',t,'xLabel','nm','yLabel','ps');

% Preview raw singular values by returning the diagonal of obj.Sr by calling obj.Srd
figure; plot(1:trSVD.nS,trSVD.Srd);
xlabel('Component #'); xlim([1, 5]);
ylabel('Singular Value');

%% Flip basis components
% Note: in general the raw SVD results will be data set specific. Adding or
% removing a spectrum can change other spectra svd flips. I recommend
% trying to come up with a rule to establish which spectra need to be
% flipped or swapped. Keep in mind this rule may be different for different
% types of data sets.
trSVD = trSVD.resetC(); %resets flips and permutes
trSVD = trSVD.permuteC(1:trSVD.nSets,[2,1]); %swaps absorptive with emissive so that absorptive is 1st and emissive is 2nd
trSVD = trSVD.flipC(1,squeeze(trSVD.Vr(end,1,:)<0),true);   % flip 1st component based on raw kinetic traces--must end as a number greater than zero (emissive growth)
trSVD = trSVD.flipC(2,squeeze(trSVD.Vr(end,2,:)>0),true); % flip 2st component based on raw kinetic traces--must end as a number less than zero (absorptive decay)

% Plot 'display' results after flipping and swapping
figure; preview(trSVD,'xVals',l,'yVals',t,'xLabel','nm','yLabel','ps');

%% Rotation Analysis
% Do the rotation analysis
[th, R] = trSVD.rotAnalysis();

% Plot the constrained spectra and kinetics to the 10th data (pH 13)
figure;
previewc(trSVD.constrain(trSVD.U(:,:,10)),'xVals',l,'yVals',t,'xLabel','nm','yLabel','ps',...
         'legend',cellfun(@(c) strcat({'pH '}, num2str(c)), num2cell(pH)));

% Plot the rotation angle w.r.t. pH 7 and ph 13
figure;
    plot(pH, 180/pi*th(:,[1,10]));
xlabel('pH');
ylabel('angle (deg)');
legend('pH 7','pH 13');