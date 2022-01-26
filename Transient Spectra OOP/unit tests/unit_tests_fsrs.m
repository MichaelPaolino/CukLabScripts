%% This is the unit test for the FSRS class

%% constructor tests
% construct default object
myFSRS = fsrs();

% load FSRS object from live tweaking data
myGsFsrs = myFSRS.loadPath('22-01-10_14h11m45stransmission_bsp_air_GS_raman_max_pump_power.mat');

% construct FSRS object from acquisition data using an alternative call
[~,~,myTR] = loadPath(fsrs(), '22-01-07_15h07m32s_Sample_BSP_Air.mat');

%% plotSpectra tests on acquisition data
%load a test data set (acquisition)
[~,~,myTR] = loadPath(fsrs(),'22-01-07_15h07m32s_Sample_BSP_Air.mat');

% plot all FSRS spectra in acquisition data
figure;
myTR.plotSpectra();

% plot a subset of FSRS spectra in the delay vector without the legend
figure;
myTR.plotSpectra('no legend');

% plot a subset of FSRS spectra in the delay vector
figure;
myTR.plotSpectra('delays',[0,0.5,1]+145.6);

% plot a subset of FSRS spectra in the delay vector with non-unique values
figure;
myTR.plotSpectra('delays',[0,0.001,0.5,1]+145.6);

% plot a subset of FSRS spectra in the delay vector while passing parameters to plot
figure;
myTR.plotSpectra('delays',[0,0.5,1]+145.6,'-k');

%% plotSpectra tests on live tweaking data
%load a test data set (live tweaking)
myFSRS = loadPath(fsrs(),'22-01-10_14h11m45stransmission_bsp_air_GS_raman_max_pump_power.mat');

% plot all FSRS spectra in live tweaking data
figure;
myFSRS.plotSpectra('no legend');

%% plotTrace tests
%load a test data set (acquisition)
[~,~,myTR] = loadPath(fsrs(),'22-01-07_15h07m32s_Sample_BSP_Air.mat');

figure;
myTR.plotTrace('wavelengths',[380,400,420,440]);

figure;
myTR.plotTrace('wavelengths',[350:5:500],'no legend');