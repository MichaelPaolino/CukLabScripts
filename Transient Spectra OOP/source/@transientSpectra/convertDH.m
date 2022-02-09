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

    %Loop through data holder array and convert nested struct
    %arrays to temporary simple arrays
    for ii = 1:nArray
        
        if isAcquisition
            rpts_tmp(ii) = dh_array(ii).Repeat;
            gpos_tmp(ii) = dh_array(ii).Grating_Position;
            delays_tmp(ii) = dh_array(ii).Delay;
        else
            rpts_tmp(ii) = ii;
        end
        
        wlc_tmp(:,ii) = dh_array(ii).wlc;
        
        for jj = 1:obj.sizes.nSchemes
            spectra_tmp(:,ii,jj) = dh_array(ii).proc_data(jj).data; %[pixel, array ind, data scheme]
            spectrastd_tmp(:,ii,jj) = dh_array(ii).proc_data(jj).data; %[pixel, array ind, data scheme]
        end
    end
    
    %Set units for delays and wavelengths
    obj.delays = obj.delays.changeBaseName('ps','Delay (ps)');
    obj.wavelengths = obj.wavelengths.changeBaseName('nm','Wavelength (nm)');
    
    %Format repeat and grating position temporary arrays into final formats
    obj.sizes.nRpts = max(rpts_tmp);
    if isAcquisition
        obj.gPos = unique(gpos_tmp);
        obj.sizes.nGPos = length(obj.gPos);
        
        %this assumes the program finished collecting data. todo: add acquisition settings to dh_static
        obj.sizes.nDelays = nArray/obj.sizes.nRpts/obj.sizes.nGPos; 
        obj.delays.data = reshape(delays_tmp, obj.sizes.nDelays, obj.sizes.nRpts, obj.sizes.nGPos); %[delays, repeats, grating positions]
        
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
    obj.schemes = cell(obj.sizes.nSchemes,1);
    for ii = 1:obj.sizes.nSchemes
       obj.schemes(ii) = {dh_array(1).proc_data.dataScheme}; 
    end
    
    %set units for spectra
    obj.spectra = obj.spectra.changeBaseName('OD','\DeltaAbs. (OD)');
    obj.spectra_std = obj.spectra_std.changeBaseName('OD','\DeltaAbs. (OD)');
    
    %finally format raw data into final format
    obj.spectra.data = reshape(spectra_tmp, obj.sizes.nPixels, obj.sizes.nDelays, obj.sizes.nRpts, obj.sizes.nGPos, obj.sizes.nSchemes); %[pixels, delays, rpts, grating pos, schemes]
    obj.spectra_std.data = reshape(spectrastd_tmp, obj.sizes.nPixels, obj.sizes.nDelays, obj.sizes.nRpts, obj.sizes.nGPos, obj.sizes.nSchemes); %[pixels, delays, rpts, grating pos, schemes]

    %Set a default short name, long name, and description
    obj.desc.description = dh_static.Description;

end