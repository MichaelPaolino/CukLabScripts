%Unit tests for transientSpectra class
%% constructors
myTS = transientSpectra();

myTS = transientSpectra('Sample_BSP_Water_rpts_schemes_gpos_delays.mat');

myTS = transientSpectra('Sample_BSP_Water_rpts_schemes_gpos_delays.mat','shortName','asdf');

myTS = transientSpectra({'Sample_BSP_Water_rpts_schemes_gpos_delays.mat','Sample_BSP_Water_rpts_schemes_gpos_delays.mat','Sample_BSP_Water_rpts_schemes_gpos_delays.mat'});

% Load by array test
inputTable = {'spectra 1', 'Sample_BSP_Water_rpts_schemes_gpos_delays.mat';...
              'spectra 2', 'Sample_BSP_Water_rpts_schemes_gpos_delays.mat';...
              'spectra 3', 'Sample_BSP_Water_rpts_schemes_gpos_delays.mat'};    

myTS = transientSpectra(inputTable(:,2),'shortName',inputTable(:,1));

%% Units
myTS = testObjArray();

%test setting units
myTS = myTS.setUnits('nm','ps','OD');

%test getting units
[wlUnit,delayUnit,sigUnit] = myTS.getUnits;
assert(ischar(wlUnit),'failed to return units as char array');
assert(strcmp(wlUnit,'nm') && strcmp(delayUnit,'ps') && strcmp(sigUnit,'OD'),'Failed to return correct units');

%% Scheme operations
myTS = testObjArray();

%test get scheme
myTS2 = myTS.getScheme('ES Raman');
assert(all(strcmp([myTS2.schemes],'ES Raman')),'Failed to return correct scheme with direct name call');
myTS2 = myTS.getScheme(1);
assert(all(strcmp([myTS2.schemes],myTS2(1).schemes)),'Failed to return correct scheme with direct name call');

%test common scheme with the same objects
schemeList = myTS.getCommonLabels('schemes');
assert(all(strcmp(schemeList,{'ES - GS Raman';'ES Raman';'Electronic Background';'GS Raman';'Transient Reflectance'})),...
    'Failed to return correct common schemes with all objects the same');

%test common scheme with mixed schemes
myTS2 = transientSpectra('Sample_None_Methanol.mat');
myTS3 = [myTS; myTS2];
schemeList = myTS3.getCommonLabels('schemes');
assert(all(strcmp(schemeList,'GS Raman')),...
    'Failed to return correct common schemes with mixed scheme object array');

%test split schemes with various calls
[myTS2,schemeList] = myTS.splitSchemes();
assert(all(size(myTS2)==[4,5]),'splitScheme failed to return correct object size');
assert(all(strcmp(schemeList,{'ES - GS Raman';'ES Raman';'Electronic Background';'GS Raman';'Transient Reflectance'})),...
    'Failed to return correct schemes for splitSchemes()');

[myTS2,schemeList] = myTS.splitSchemes('-search','raman','Transient Reflectance');
assert(all(size(myTS2)==[4,4]),'splitScheme failed to return correct object size');
assert(all(strcmp(schemeList,{'ES - GS Raman';'ES Raman';'GS Raman';'Transient Reflectance'})),...
    'Failed to return correct schemes for splitSchemes -search and transient reflectance');

[myTS2,schemeList] = myTS.splitSchemes('-search','raman','-drop');
assert(all(size(myTS2)==[4,2]),'splitScheme failed to return correct object size');
assert(all(strcmp(schemeList,{'Electronic Background';'Transient Reflectance'})),...
    'Failed to return correct schemes for splitSchemes -search -drop');

%% Test getLabel()
myTS = testObjArray();
myTS = myTS.getLabel('gPos',1);
assert(all(size(myTS(1).spectra.data)==[1340,19,2,1,5]),'Failied to return correct object size');
assert(myTS(1).gPos == 400,'Failed to return correct grating position');
assert(myTS(2).gPos == 400,'Failed to return correct grating position');



%% Class conversion tests
% Class conversion tests
myFSRS = fsrs('22-10-21_10h07m42s_MeOH_Fluence_8p75_uJ,_initial_methanol_FSRS_spectrum','loadType','tsObj');
myTS = transientSpectra(repmat(myFSRS,5,1));

figure; myFSRS.plotSpectra();
figure; myTS.plotSpectra();

%% All tests pass
disp('All tests passed!');

%% support functions
function obj = testObjArray()
    inputTable = {'spectra 1', 'Sample_BSP_Water_rpts_schemes_gpos_delays.mat';...
                  'spectra 2', 'Sample_BSP_Water_rpts_schemes_gpos_delays.mat';...
                  'spectra 3', 'Sample_BSP_Water_rpts_schemes_gpos_delays.mat';...
                  'spectra 4', 'Sample_BSP_Water_rpts_schemes_gpos_delays.mat'};

    obj = transientSpectra(inputTable(:,2),'shortName',inputTable(:,1));
end