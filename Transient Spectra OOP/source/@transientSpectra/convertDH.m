function obj = convertDH(obj, dh_static, dh_array)
%converts a data_holder object to a FSRS object. Need to come up with some
%strategy for managing all the different versions until a standardized file
%format is created for live tweaking and acquisition. May use the availalbe
%fields as a version signature for now with specific implementations for
%different versions.

    %first determine whether this is a live tweaking or acquisition
    %data_holder
    if ~isfield(dh_array,'Repeat')
        isAcquisition = false;
    else
        isAcquisition = true;
        if ~isfield(dh_static,'SpectCalib')
            labelSpectCalib = 'x';
        else
            labelSpectCalib = 'SpectCalib';
        end
    end
    
    %get dh_array information
    nArray = length(dh_array);
    obj.sizes.nPixels = length(dh_array(1).wlc);
    obj.sizes.nSchemes = length(dh_array(1).proc_data);

    %initialize tmp variables for storing spectral information
    rpts_tmp = zeros(nArray, 1);
    if isAcquisition
    	gpos_tmp = zeros(nArray, 1);
        delays_tmp = zeros(nArray, 1);
    end
    
    wlc_tmp = zeros(obj.sizes.nPixels, nArray);
    spectra_tmp = zeros(obj.sizes.nPixels, nArray, obj.sizes.nSchemes);
    spectrastd_tmp = zeros(obj.sizes.nPixels, nArray, obj.sizes.nSchemes);
    
    %get indicies of field names for dh_array proc_data inside dh_array
    dhFields = fieldnames(dh_array);
    procDataFields = fieldnames(dh_array(1).proc_data);
    
    %convert nested double data inside struct to cell array and then to matrix
    %this is 2 order of magnitude faster than doing it in a for loop...
    cellDHArray = struct2cell(dh_array); 
    cellProcData = cellfun(@(s) struct2cell(s),cellDHArray(strcmp(dhFields,'proc_data'),1,:),'UniformOutput',false); 
    cellProcData = reshape(vertcat(cellProcData{:}),[],nArray,obj.sizes.nSchemes);
    
    %extract cell data into matrix
    %wlc_tmp = permute(cell2mat(cellDHArray(strcmp(dhFields,'wlc'),1,:)),[2,3,1]); %[pixel, array ind]
    spectra_tmp = reshape(cell2mat(cellProcData(strcmp(procDataFields,'data'),:,:)),[],nArray,obj.sizes.nSchemes); %[pixel, array ind, data scheme]
    spectrastd_tmp = reshape(cell2mat(cellProcData(strcmp(procDataFields,'std'),:,:)),[],nArray,obj.sizes.nSchemes); %[pixel, array ind, data scheme]
    
    if isAcquisition
        rpts_tmp = permute(cell2mat(cellDHArray(strcmp(dhFields,'Repeat'),1,:)),[2,3,1]); %[1, array ind]
        gpos_tmp = permute(cell2mat(cellDHArray(strcmp(dhFields,'Grating_Position'),1,:)),[2,3,1]); %[1, array ind]
        delays_tmp = permute(cell2mat(cellDHArray(strcmp(dhFields,'Delay'),1,:)),[2,3,1]); %[1, array ind]
    else
        rpts_tmp = 1:nArray;
    end
    
%     %Loop through data holder array and convert nested struct
%     %arrays to temporary simple arrays
%     for ii = 1:nArray
%         
%         if isAcquisition
%             rpts_tmp(ii) = dh_array(ii).Repeat;
%             gpos_tmp(ii) = dh_array(ii).Grating_Position;
%             delays_tmp(ii) = dh_array(ii).Delay;
%         else
%             rpts_tmp(ii) = ii;
%         end
%         
%         wlc_tmp(:,ii) = dh_array(ii).wlc;
%         
%         for jj = 1:obj.sizes.nSchemes
%             spectra_tmp(:,ii,jj) = dh_array(ii).proc_data(jj).data; %[pixel, array ind, data scheme]
%             spectrastd_tmp(:,ii,jj) = dh_array(ii).proc_data(jj).data; %[pixel, array ind, data scheme]
%         end
%     end
    
    %define unit rules for spectra, delays, and wavelengths.
    %spectra
    spectraRules = doubleWithUnits([],'OD','\DeltaAbs. (OD)');
    spectraRules = spectraRules.addRule('mOD','\DeltaAbs. (mOD)',@(f) 1e3*f, @(f) 1e-3*f);
    
    %delays
    delayRules = doubleWithUnits([],'ps','Delay (ps)');
    delayRules = delayRules.addRule('fs','Delay (fs)',@(f) 1e3*f, @(f) 1e-3*f);
    delayRules = delayRules.addRule('ns','Delay (ns)',@(f) 1e-3*f, @(f) 1e3*f);
    delayRules = delayRules.addRule('us','Delay (\ms)',@(f) 1e-6*f, @(f) 1e6*f);
    
    %wavelengths
    wlRules = doubleWithUnits([],'nm','Wavelength (nm)');            
    wlRules = wlRules.addRule('um','Wavelength (\mm)',@(f) 1e3*f, @(f) 1e-3*f);
    wlRules = wlRules.addRule('eV','Energy (eV)',@(f) 1239.8./f, @(f) 1239.8./f);
    wlRules = wlRules.addRule('ecm-1','Wavenumber (cm^{-1})',@(f) 1e7./f, @(f) 1e7./f);
    
    % Assign unit rules to object data, starting with an empty doubleWithUnits
    obj.spectra = doubleWithUnits([],spectraRules);
    obj.spectra_std = doubleWithUnits([],spectraRules);
    obj.delays = doubleWithUnits([],delayRules);
    obj.t0 = doubleWithUnits(0,delayRules);
    obj.wavelengths = doubleWithUnits([],wlRules);
    
    %Format repeat and grating position temporary arrays into final formats
    obj.sizes.nRpts = max(rpts_tmp);
    if isAcquisition
        obj.gPos = unique(gpos_tmp);
        obj.sizes.nGPos = length(obj.gPos);
        
        %Add delays to object inclusive of incomplete data acquisition
        %todo: explicity assign number of delays in acquisition program
        obj.sizes.nDelays = sum(and((gpos_tmp==gpos_tmp(1)),rpts_tmp==rpts_tmp(1)));
        obj.delays.data = nan(obj.sizes.nDelays*obj.sizes.nRpts*obj.sizes.nGPos,1);
        obj.delays.data(1:length(delays_tmp)) = delays_tmp;
        obj.delays.data = reshape(obj.delays.data, obj.sizes.nDelays, obj.sizes.nRpts, obj.sizes.nGPos); %[delays, repeats, grating positions]
        
        %loop through wavelength calibrations for multiple grating positions
        obj.wavelengths.data = zeros(obj.sizes.nPixels, obj.sizes.nGPos);
        for ii = 1:obj.sizes.nGPos
            %assumes that the grating positions are in the same order as they were acquired
            obj.wavelengths(:,ii) = dh_static.(labelSpectCalib)(ii).Wavelengths;   
        end
    else
        obj.gPos = dh_static.Grating_Position;
        obj.sizes.nGPos = length(obj.gPos);
        obj.delays.data = dh_static.Delay;
        obj.sizes.nDelays = 1;
        obj.wavelengths.data = dh_static.Wavelengths(:);
    end

    %get list of acquired schemes
    obj.schemes = {dh_array(1).proc_data.dataScheme};
    obj.schemes = obj.schemes(:);
    
    %finally format raw data into final format
    obj.spectra.data = nan(obj.sizes.nPixels,obj.sizes.nDelays*obj.sizes.nRpts*obj.sizes.nGPos,obj.sizes.nSchemes);
    obj.spectra.data(:,1:size(spectra_tmp,2),:) = spectra_tmp;
    obj.spectra.data = reshape(obj.spectra.data, obj.sizes.nPixels, obj.sizes.nDelays, obj.sizes.nRpts, obj.sizes.nGPos, obj.sizes.nSchemes); %[pixels, delays, rpts, grating pos, schemes]
    
    obj.spectra_std.data = nan(obj.sizes.nPixels,obj.sizes.nDelays*obj.sizes.nRpts*obj.sizes.nGPos,obj.sizes.nSchemes);
    obj.spectra_std.data(:,1:size(spectrastd_tmp,2),:) = spectrastd_tmp;
    obj.spectra_std.data = reshape(obj.spectra_std.data, obj.sizes.nPixels, obj.sizes.nDelays, obj.sizes.nRpts, obj.sizes.nGPos, obj.sizes.nSchemes); %[pixels, delays, rpts, grating pos, schemes]

    %Set a default short name, long name, and description
    obj.description = dh_static.Description;
    
    %Set default units
    obj = obj.setUnits('nm','ps','mOD');
end