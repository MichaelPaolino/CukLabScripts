%% Load non-phonon removed raw data
conn = database('Raman DB','','');

% pH 13 and 14 phosphate data acquired in fall 2022
dataList = conn.fetch(strjoin(...
    ["SELECT ID, FileName, FilePath, ShortName, Description",...
     "FROM TRAcquisition",...
     "WHERE ID IN (787, 799)",...
     "ORDER BY ID;"],' '));
 
conn.close();

%% Load phosphate data into object
% Extract pH values and sort in ascending order
c = regexp(dataList.ShortName,'(?<=pH )\d+\.?\d*','match');
[pH, idx] = sort(cellfun(@(c) str2double(c), vertcat(c{:})));
dataList = dataList(idx,:);

% Load data
trRaw = wlTR(fullfile(repoPath(), dataList.FilePath, dataList.FileName),'shortName',strcat({'pH '}, vertcat(c{:})));

% Load chirp parameters
load(fullfile(repoPath(),'CukResearchGroup','TR-VIS','chirpFitParam22Jul28.mat'));

%% Check for pruning 
figure; trRaw.plotKinetics('wavelengths',560,'average',false,'index',true); %nothing to prune here!

%% Verify chirp parameters...
% This is optional. I reccommend using OC spectra to determine chirp
% parameters. Below is just a crazy check to make sure chirp is not
% completely off.
[trRaw, ~] = trRaw.fitChirp('wavelengths',[375:5:700],...
    'fitFun',@(x,y) sum(lsqFitExpGrow(x, y).*[0,0,0,0,0,1]));

figure;
l = 375:700;
plot(l,polyval(chirpFit,l)-polyval(chirpFit,560),...
     l,polyval(trRaw(1).chirpParams,l)-polyval(trRaw(1).chirpParams,560),...
     l,polyval(trRaw(2).chirpParams,l)-polyval(trRaw(2).chirpParams,560));
 xlabel('Wavelength (nm)');
 ylabel('Group Delay (ps)');
 legend('OC','pH 13','pH 14');

%% Remove phonons
% Use OC chirp parameters collected earlier from file. Using fit above
% doesn't make a big difference.
trRaw = trRaw.setChirp(chirpFit);

% Calculate and remove phonons. See documentation (f1 on removePhonons) for
% options. Default works well for 0V CC data SrTiO3.
[trPR, phononFit] = trRaw.removePhonons();

%% Cosmetic processing
% Append raw data and phonon removed data and suffix phonon removed data
% with phonon removed
tr = [trRaw, trPR];
tr(:,2) = tr(:,2).setShortName(strcat(tr(:,2).getShortName, {' phonon removed'}));

% Average, stitch, correct chirp
tr = tr.average();
tr = tr.stitch();
tr = tr.correctChirp();

% Slightly more accurate version of finding t0--does not improve things
tr = tr.findT0('wlAr',400:50:700,'fitFun',@(x,y) sum(lsqFitExpGrow(x, y).*[0,0,0,0,0,1]));
tr = tr.correctT0();

% Interpolate all object data on the same delay so that plots below show same delays
tr = tr.interp('delays',tr(1,1).delays.data);

%% Plot results
% Kinetics to check quality of phonon removal and early-time CAWs subtraction
figure;
tr.plotKinetics('wavelengths',[400,600]);
xlim([-5,100]);

% Spectra to check quality of phonon removal--note phase offset between pH
% 13 and pH 14 in spectra. This is not caused by t0, maybe spectrometer
% calibration or interpolation?
figure;
tr.plotSpectra('delays',[100,500]);

% contours
figure;
cvals = -5:0.25:2; % contour values
for objInd = 1:numel(tr)
    subplot(2,2,objInd);
        contourf(tr(objInd).wavelengths.data,tr(objInd).delays.data,tr(objInd).spectra.data',cvals);
        title(tr(objInd).shortName);
        caxis([min(cvals), max(cvals)]);
        xlim([375,700]); xlabel('wavelengths (nm)');
        ylim([0.1,4000]);  ylabel('delay (ps)');
        set(gca,'YScale','log');
end

%% Export data to Igor
savePath = fullfile(repoPath(),'CukResearchGroup','Software','Acquisition and Processing Software','Examples','TR-VIS','Procd','Files and Scripts','22_11_16Nov pH 13 14 example','Igor Files');
tr.export(fullfile(savePath,'example_pH13_pH14_phononRmvd_Igor.mat'),'kinetics',[400,600],'spectra',[10,100,500,1000]);

%% Export data to database--Run first and double check stage tables!!!
% Establish ODBC connection to data source defined above
conn = database('Raman DB', '','');
savePath = fullfile(repoPath(),'CukResearchGroup','Software','Acquisition and Processing Software','Examples','TR-VIS','Procd','Files and Scripts','22_11_16Nov pH 13 14 example','Procd Data');

% assign foreign key values
fkCols = genFKStruct('UserID','Users','ID','FirstName','Ilya');

% non-phonon removed data
procdTable = tr(:,1).dbStageSave(conn,'TRProcd',...
                                'path',savePath,...
                                'prefix','Example',...
                                'fkStruct',fkCols);
                                %'update',[168 169]);    %Use this line when updating the database instead of inserting

% phonon removed data
procdTablePR = tr(:,2).dbStageSave(conn,'TRPhononRemoved',...
                                'path',savePath,...
                                'prefix','Example',...
                                'fkStruct',fkCols);
                                %'update',[120 121]);    %Use this line when updating the database instead of inserting

%% Save and commit data--Run only after double checking stage tables!!
[commitTable, IDStr] = tr(:,1).dbCommitSave(conn,procdTable);
[commitTablePR, IDStrPR] = tr(:,2).dbCommitSave(conn,procdTablePR);

% display the IDs that were generated so that the user can copy/paste them into the dbStageSave 'update' arguments 
disp(IDStr);
disp(IDStrPR);

%% Export commited data to associative tables--Run first and double check stage tables!!!
[assocStage, assocView] = tr(:,1).dbStageAssoc(conn,'assocTRAcquisitionProcd','TRAcquisition','TRProcd',commitTable);
[assocStagePR, assocViewPR] = tr(:,2).dbStageAssoc(conn,'assocTRAcquisitionPhononRemoved','TRAcquisition','TRPhononRemoved',commitTablePR);

%%  Commit associative data--Run only after double checking stage tables!!
[assocCommit, IDStrAssoc] = tr(:,1).dbCommitAssoc(conn, assocStage);
[assocCommitPR, IDStrAssocPR] = tr(:,2).dbCommitAssoc(conn, assocStagePR);

% display the IDs that were generated so that the user can copy/paste them into the dbStageSave 'update' arguments 
disp(IDStrAssoc);
disp(IDStrAssocPR);

%% Close the ODBC connection
conn.close