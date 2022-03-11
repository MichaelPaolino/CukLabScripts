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
    
    %set units for spectra
    obj.spectra = obj.spectra.changeBaseName('OD','\DeltaAbs. (OD)');
    obj.spectra_std = obj.spectra_std.changeBaseName('OD','\DeltaAbs. (OD)');
    
    %finally format raw data into final format
    obj.spectra.data = nan(obj.sizes.nPixels,obj.sizes.nDelays*obj.sizes.nRpts*obj.sizes.nGPos,obj.sizes.nSchemes);
    obj.spectra.data(:,1:size(spectra_tmp,2),:) = spectra_tmp;
    obj.spectra.data = reshape(obj.spectra.data, obj.sizes.nPixels, obj.sizes.nDelays, obj.sizes.nRpts, obj.sizes.nGPos, obj.sizes.nSchemes); %[pixels, delays, rpts, grating pos, schemes]
    
    obj.spectra_std.data = nan(obj.sizes.nPixels,obj.sizes.nDelays*obj.sizes.nRpts*obj.sizes.nGPos,obj.sizes.nSchemes);
    obj.spectra_std.data(:,1:size(spectrastd_tmp,2),:) = spectrastd_tmp;
    obj.spectra_std.data = reshape(obj.spectra_std.data, obj.sizes.nPixels, obj.sizes.nDelays, obj.sizes.nRpts, obj.sizes.nGPos, obj.sizes.nSchemes); %[pixels, delays, rpts, grating pos, schemes]

    %Set a default short name, long name, and description
    obj.description = dh_static.Description;

end