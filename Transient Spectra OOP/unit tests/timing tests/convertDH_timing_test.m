%old load method
%fsrs('Sample_S_pH_13100mMNaOH_incomplete');

%new load method
%sizes
nArray = length(dh_array);
nSchemes = length(dh_array(1).proc_data);

%get indicies of field names for dh_array
dhFields = fieldnames(dh_array);
procDataInd = strcmp(dhFields,'proc_data');
wlcInd = strcmp(dhFields,'wlc');

%get indicies of field names for proc_data inside dh_array
procDataFields = fieldnames(dh_array(1).proc_data);
spectraField = strcmp(procDataFields,'data');
spectraStdField = strcmp(procDataFields,'std');

%get indicies of field names related to iteration info
rptInd = strcmp(dhFields,'Repeat');
gposInd = strcmp(dhFields,'Grating_Position');
delayInd = strcmp(dhFields,'Delay');

%convert nested double data inside struct to cell array and then to matrix
%this is ~2 order of magnitude faster than doing it in a for loop...
cellDHArray = struct2cell(dh_array); 
cellProcData = cellfun(@(s) struct2cell(s),cellDHArray(procDataInd,1,:),'UniformOutput',false); 
cellProcData = reshape(vertcat(cellProcData{:}),[],nArray,nSchemes);

wlc_tmp = permute(cell2mat(cellDHArray(wlcInd,1,:)),[2,3,1]); %[pixel, array ind]
spectra_tmp = reshape(cell2mat(cellProcData(spectraField,:,:)),[],nArray,nSchemes); %[pixel, array ind, data scheme]
spectrastd_tmp = reshape(cell2mat(cellProcData(spectraStdField,:,:)),[],nArray,nSchemes); %[pixel, array ind, data scheme]

rpts_tmp = permute(cell2mat(cellDHArray(rptInd,1,:)),[2,3,1]); %[1, array ind]
gpos_tmp = permute(cell2mat(cellDHArray(gposInd,1,:)),[2,3,1]); %[1, array ind]
delays_tmp = permute(cell2mat(cellDHArray(delayInd,1,:)),[2,3,1]); %[1, array ind]

%%