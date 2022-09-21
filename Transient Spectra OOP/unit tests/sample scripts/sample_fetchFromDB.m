%% Define data repo path
% google drive location for CukResearchGroup folder
dataPath = 'G:\.shortcut-targets-by-id\1O21NYESEEKxBx8k6l21vr3Feo26hcErZ\';

% ODBC data source name for raman.accdb
dbDSN = 'Raman Access';

%% Connect to database and query data
% Establish ODBC connection to data source defined above
conn = database(dbDSN, '','');

% fetch data using an SQL query
dataList = conn.fetch(['SELECT FilePath, FileName, ShortName '...
                       'FROM TRPhononRemoved '...
                       'WHERE ShortName LIKE ''%pH 13 %'' ',...
                       'AND ShortName NOT LIKE ''%OC%'' ',...
                       'AND ShortName NOT LIKE ''%D2O%'' ',...
                       'AND (ShortName NOT LIKE ''%mJ/cm2%'' OR ShortName LIKE ''%0.04 mJ/cm2%'') ',...
                       'AND (ShortName NOT LIKE ''%[%]%'' OR ShortName LIKE ''%0.1[%]%'');']);

% Close the ODBC connection
conn.close;

%% Load data into wlTR object
pH13Data = wlTR(fullfile(dataPath,dataList.FilePath,dataList.FileName),...
                'loadType','cListFile','shortName',dataList.ShortName);

%% Plot TR spectra at 10 ps
figure;
pH13Data.plotSpectra('delays',10);

%% Plot TR traces at 400 nm
figure;
subplot(1,2,1);
    pH13Data.plotKinetics('wavelengths',400,'legend',false);
    xlim([-5,50]);
    ylim([-3,0.5]);
subplot(1,2,2);
    pH13Data.plotKinetics('wavelengths',400);
    set(gca,'XScale','log');
    set(gca,'YTick',[]);
    ylabel('');
    xlim([50,4000]);
    ylim([-3,0.5]);