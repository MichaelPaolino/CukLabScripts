%% Define data and load into fsrs object array
% Google drive path
gDrive = 'G:\.shortcut-targets-by-id\1O21NYESEEKxBx8k6l21vr3Feo26hcErZ\';
dataFolder = 'CukResearchGroup\FSRS\Raw\22_02_08Feb\';

% Use a cell array to load multiple data sets at once.
dataList = {'8.0 uJ','22-02-08_14h02m19s_Sample_BSP_Air.mat'
            '6.0 uJ','22-02-08_14h43m22s_Sample_BSP_Air.mat';...
            '4.0 uJ','22-02-08_14h48m43s_Sample_BSP_Air.mat';...
            '2.0 uJ','22-02-08_14h59m03s_Sample_BSP_Air.mat';...
            '1.0 uJ','22-02-08_15h10m21s_Sample_BSP_Air.mat';...
            '0.51 uJ','22-02-08_15h15m39s_Sample_BSP_Air.mat';...
            '0.21 uJ','22-02-08_15h20m54s_Sample_BSP_Air.mat';...
            '0.11 uJ','22-02-08_15h24m00s_Sample_BSP_Air.mat'
            };

% Use fullfile to append file name cell array to to gDrive string (char array)
% Use array () indexing instead of cell {} indexing to return a cell array
% any transientSpectra and descendent class accepts cell arrays as inputs
% and returns an object array that is the same size as the input cell array
myFSRS = fsrs(fullfile(gDrive,dataList(:,2)),'shortName',dataList(:,1));

%% Alternatively, one could load by query
gDrive = 'G:\.shortcut-targets-by-id\1O21NYESEEKxBx8k6l21vr3Feo26hcErZ\';

% You may need to change this depending on how you set up the 64-bit ODBC DSN for the MS access database
dbDSN = 'Raman Access';  

% This query selects ID, Filename, Filepath and Description columns but
% parses description for the raman pump fluence and calls the resutls
% ShortName. The Query filters results so that they are:
% 1. 400 nm fluence is in the description
% 2. the last unique entry is returned
myQuery = ['SELECT ID, FileName, FilePath, MID(Description, 21, 4) AS ShortName ',...
           'FROM FSRSAcquisition ',...
           'WHERE Description LIKE ''400 nm fluence%''' ,...
           'AND ID IN (SELECT MAX(ID) FROM FSRSAcquisition GROUP BY Description);'];

% Open the database connection, execute the query above, close the database connection
conn = database(dbDSN, '','');
dataList = conn.fetch(myQuery);
conn.close();

% 1. Load using table. table.[fieldName] returns a cell array of the [fieldName]
% 2. fullfile concatonates the different path parts into one path, and is smart enough to handle cell arrays
% 3. strtrim trims white space in string cell array
% 4. strcat concatonates string cell array with another string
myFSRS = fsrs(fullfile(gDrive, dataList.FilePath, dataList.FileName),...
              'shortName', strcat(strtrim(dataList.ShortName), ' uJ'));

%% Costmetic processing
myFSRS = myFSRS.stitch();
myFSRS = myFSRS.average();
myFSRS = myFSRS.findRamanPumpNm();
myFSRS = myFSRS.trim('wavelengths',[-50,3800]);

%% Plot results after cosmetic processing
figure;
myFSRS.plotSpectra()

%% Correct the baseline--maybe this should be a method in the fsrs class?
fsrsBLind = myFSRS(:); %This will hold unnomalized baseline-corrected data
fsrsBLindNmd = fsrsBLind; %This will hold normalized baseline-corrected data

for objInd = 1:numel(fsrsBLind)
   % manually extract data from fsrs object for custom processing
   y = fsrsBLind(objInd).spectra.data(:); % sample spectra
   k = ones(size(y)); % 'background' spectra--ones because we don't know what the true background is
   
   % calculate the subtractor y' = y - x*k where x is a multiplicative error in baseline
   x = AslsSolventSubt(y,k,'maxRec',25);
   
   % update the fsrs data and normalize by fluence
   fsrsBLind(objInd).spectra.data = (y - x.*k);
   fsrsBLindNmd(objInd).spectra.data = fsrsBLind(objInd).spectra.data/str2double(dataList.ShortName{objInd});
end

% Plot unnormalized
figure;
fsrsBLind.plotSpectra()

% Plot normalized
figure;
fsrsBLindNmd(1:6).plotSpectra()