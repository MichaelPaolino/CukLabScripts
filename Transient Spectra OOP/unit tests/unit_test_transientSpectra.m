%Unit tests for transientSpectra class
%% constructors
myTS = transientSpectra();

myTS = transientSpectra('Sample_BSP_Water_rpts_schemes_gpos_delays.mat');

myTS = transientSpectra('Sample_BSP_Water_rpts_schemes_gpos_delays.mat','shortName','asdf');

myTS = transientSpectra({'Sample_BSP_Water_rpts_schemes_gpos_delays.mat','Sample_BSP_Water_rpts_schemes_gpos_delays.mat','Sample_BSP_Water_rpts_schemes_gpos_delays.mat'});

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
schemeList = myTS.getCommonSchemes;
assert(all(strcmp(schemeList,{'ES - GS Raman';'ES Raman';'Electronic Background';'GS Raman';'Transient Reflectance'})),...
    'Failed to return correct common schemes with all objects the same');

%test common scheme with mixed schemes
myTS2 = transientSpectra('Sample_None_Methanol.mat');
myTS3 = [myTS; myTS2];
schemeList = myTS3.getCommonSchemes;
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

%% support functions
function obj = testObjArray()
    inputTable = {'spectra 1', 'Sample_BSP_Water_rpts_schemes_gpos_delays.mat';...
                  'spectra 2', 'Sample_BSP_Water_rpts_schemes_gpos_delays.mat';...
                  'spectra 3', 'Sample_BSP_Water_rpts_schemes_gpos_delays.mat';...
                  'spectra 4', 'Sample_BSP_Water_rpts_schemes_gpos_delays.mat'};

    obj = transientSpectra(inputTable(:,2),'shortName',inputTable(:,1));
end