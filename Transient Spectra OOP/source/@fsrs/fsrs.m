%% thoughts
%data holder strategies:
%1. use child classes to specify types of data that will contain experiment
%specific data access methods. Question to be addressed is how data will be stored in
%the intermediate when there is data processing?
%   a. FSRS_raw -> used to extract raw data from dh class
%   experiment
%   b. FSRS_proc? -> used to store and extract processed data
%   format
%   c. TR_raw -> used to extract raw data from dh class
%   d. TR_proc -> used to store and extract processed data
%2. **try this one first** create a new class that takes data_holder as an input and has its own
%member data
%   a. FSRS -> constructor converts dh data into FSRS spectra. Class contains
%   override methods for common functions, such as plot, contour, etc
%   b. TR -> similar to FSRS...
%
%3. create a general-purpose SQL base extraction method

classdef fsrs < transientSpectra
    properties

        ramanPumpNm = 400;  %raman pump wavelength. todo: add conversion method to wavenumbers
        %implement conversions to plot methods
        

    end
   
    methods
        %%LOAD METHODS%%

        %use transientSpectra constructor. Matlab calls this implicitly
        
        %override superclass convertDH method to convert multi-scheme data
        %into specific scheme objects
        function [objGsFsrs, objEsFsrs, objTR] = convertDH(obj, dh_static, dh_array)
            
            %call parent method first for generic conversion of dh_array to transientSpectra
            obj = convertDH@transientSpectra(obj, dh_static, dh_array);
            
            %group various schemes into new objects
            schemeList = {'GS Raman','ES Raman','Transient Reflectance'};
            nSchemes = length(schemeList);
            tsArray = repmat(fsrs(),nSchemes,1);
            
            %todo: redesign how this works...
            for ii = 1:nSchemes
                loc = strcmp(obj.schemes,schemeList(ii));
                if any(loc)
                    tsArray(ii) = obj;
                    tsArray(ii).spectra = obj.spectra(:,:,:,:,loc);
                    tsArray(ii).spectra_std = obj.spectra(:,:,:,:,loc);
                end
            end
            
            objGsFsrs = tsArray(1);
            objEsFsrs = tsArray(2);
            objTR = tsArray(3);
        end
        
    end
end

