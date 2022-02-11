classdef fsrs < transientSpectra
    properties
        ramanPumpNm = 400;  %raman pump wavelength in nm
        %implement conversions methods
    end

    % Constructor, load, get and set methods with custom implementation
    methods
        %%CONSTRUCTOR/LOAD METHODS%%

        %use transientSpectra constructor. Matlab calls this implicitly
        
        %override superclass convertDH method to convert multi-scheme data
        %into specific scheme objects. Call this to load multiple objects
        function [objGsFsrs, objEsFsrs, objTR] = convertDH(obj, dh_static, dh_array)
            
            %call parent method first for generic conversion of dh_array to transientSpectra
            obj = convertDH@transientSpectra(obj, dh_static, dh_array);
            
            %add unit rules to spectra, delays, etc. todo: add to method
            obj.spectra = obj.spectra.addRule('mOD','\DeltaAbs. (mOD)',...
                                                @(f) 1e3*f, @(f) 1e-3*f);
            obj.spectra = obj.spectra.addRule('%Gain','Raman Gain (%)',...
                                                @(f) 1e2*(10.^f-1), @(f) log10(1+1e-2*f));
            obj.spectra = obj.spectra.addRule('ppmGain','Raman Gain (ppm)',...
                                                @(f) 1e6*(10.^f-1), @(f) log10(1+1e-6*f));
            
            obj.delays = obj.delays.addRule('fs','Delay (fs)',@(f) 1e3*f, @(f) 1e-3*f);
            obj.delays = obj.delays.addRule('ns','Delay (ns)',@(f) 1e-3*f, @(f) 1e3*f);
            obj.delays = obj.delays.addRule('us','Delay (\ms)',@(f) 1e-6*f, @(f) 1e6*f);
            
            obj.wavelengths = obj.wavelengths.addRule('um','Wavelength (\mm)',....
                                                        @(f) 1e3*f, @(f) 1e-3*f);
            obj.wavelengths = obj.wavelengths.addRule('eV','Energy (eV)',...
                                                        @(f) 1239.8./f, @(f) 1239.8./f);
            obj.wavelengths = obj.wavelengths.addRule('ecm-1','Wavenumber (cm^{-1})',...
                                                        @(f) 1e7./f, @(f) 1e7./f);
            obj.wavelengths = obj.wavelengths.addRule('rcm-1','Raman Shift (cm^{-1})',...
                                                        @(f) 1e7*(1/obj.ramanPumpNm-1./f),...
                                                        @(f) 1./(1/obj.ramanPumpNm-1e-7*f));
            
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
            
            %set output objects
            objGsFsrs = tsArray(1);
            objEsFsrs = tsArray(2);
            objTR = tsArray(3);
            
            %set output object units
            objGsFsrs = objGsFsrs.setUnits('rcm-1','ps','%Gain'); %change to raman shift, ps, and raman gain
            objEsFsrs = objGsFsrs.setUnits('rcm-1','ps','%Gain'); %change to raman shift, ps, and raman gain
            objTR = objTR.setUnits('eV','ps','mOD'); %change to raman shift, ps, and raman gain
            
        end
        
        %%GET/SET METHODS%% 
        
        %Set the raman pump nm. This updates the raman shift unit
        %definition, which is why it requires an explicit set method.
        function obj = set.ramanPumpNm(obj,newNm)
            obj.ramanPumpNm = newNm;
            %todo: figure out how to do this for all properties that might
            %have unit changes. Maybe cosmetic struct has .wavelengths,
            %.delays and .spectra properties
            obj.wavelengths = obj.wavelengths.updateRule('rcm-1','Raman Shift (cm^{-1})',...
                @(f) 1e7*(1/obj.ramanPumpNm-1./f),...
                @(f) 1./(1/obj.ramanPumpNm-1e-7*f));
        end
        
    end
    
    % Methods specific to the fsrs class that cannot be implemented in the transientSpectra class
    methods
        
        function [obj, pumpNm] = findRamanPumpNm(obj, varargin)
        % FINDRAMANPUMPNM automatically finds the raman pump wavelength in nm.
        % The algorithm works by searching a wavelength region defined by guess +/- 
        % threshold and searching for the maximum point. Because the fundamental
        % peak is an image of the pump scatter on the spectrometer slit, it does 
        % not have a clean shape. Instead of fitting the peak, the algroithm 
        % attempts to determine the peak center by its FWHM. If the FWHM cannot be 
        % resolved, the algorithm takes the maximum point as the center. The final 
        % pump wavelength is the average over all valid grating positions.
        %
        % By default, the guess wavelength is obj.ramanPumpNm and the threshold is
        % 10 nm. These defaults can be overriden. Note: the search is performed in
        % nm units, but the output object will have the same units as the input.
        %
        % obj = obj.FINDRAMANPUMPNM()
        %   Search and update the raman pump center wavelength (nm) using the
        %   default search region of obj.ramanPumpNm +/- 10 nm.
        %
        % obj = obj.FINDRAMANPUMPNM('name1', value1, 'name2', value2,...)
        %   Search and update the raman pump center wavelength (nm) with options
        %   set by name-value pair. Name-value paris are:
        %   'pump guess', newGuess: double newGuess updates the search region
        %       center to newGuess +/- threshold, where threshold defaults to 10 nm
        %   'threshold', newThresh: double newThresh updates the search region
        %       range to currentGuess +/- newThresh, where currentGuess defaults to
        %       obj.ramanPumpNm
        %   'drop', dropVector: double dropVector contains grating positions that
        %       are force-excluded from the raman pump wavelength search. By
        %       default, all gratings that overlap with the search range are
        %       searched.
           
            %Default parameters
            pumpGuessNm = obj.ramanPumpNm;    %initial guess for search in nm
            thresholdNm = 10; %+/- search range around initial guess in nm
            dropGPos = [];  %grating positions to automatically exclude
            
            %parse varargin
            if nargin > 1
                for ii = 1:2:(nargin-1)
                    assert(ischar(varargin{ii}),...
                        ['Invalid argument class for name-value pair. Expected class char for name, got ' class(varargin{ii}) '.']);
                    switch varargin{ii}
                        case 'pump guess'
                            assert(isscalar(varargin{ii+1}),'Expected a scalar value for pump nm guess');
                            pumpGuessNm = varargin{ii+1};
                        case 'threshold'
                            assert(isscalar(varargin{ii+1}),'Expected a scalar value for threshold');
                            thresholdNm = varargin{ii+1};
                        case 'drop'
                            assert(isvector(varargin{ii+1}),'Expected a vector for drop grating positions');
                            dropGPos = varargin{ii+1};
                        otherwise
                            error([varargin{ii} ' is not a valid argument name.']); 
                    end
                end
            end
            
            %remember old units
            tmpUnits = cell(3,1);
            [tmpUnits{:}] = obj.getUnits();
            
            %update units to units where x-axis is in nm and raman pump scatter is positive
            obj = obj.setUnits('nm',[],'%Gain');
            
            subObj = obj.trim('wavelengths',pumpGuessNm+thresholdNm*[-1,1]);
            
            %get a single spectra by averaging over any extra dims except grating position (up to 5)
            data = permute(mean(subObj.spectra.data,[2, 3, 5]),[1,4,2,3,5]); 
            
            %sort data and wavelengths by ascending order in grating position
            [gPosSorted,gInd] = sort(subObj.gPos);
            data = data(:,gInd);
            lambda = subObj.wavelengths.data(:,gInd);
            
            %find the raman pump peak (should be strongest signal)
            [maxVal,maxInd] = max(data);
            
            %find peak base (should be close to lowest signal)
            minVal = min(data);
            
            %find pump peak fwhm while looping over sorted grating positions that do not get dropped
            pumpNm = NaN*zeros(1,subObj.sizes.nGPos);
            for ii = 1:subObj.sizes.nGPos
                %perform the fwhm routine only if the grating position is not dropped
                if ~any(dropGPos == gPosSorted(ii))
                    %find all "fwhm" intercepts in data. Sort in case lambda nm values
                    %are not sorted in ascending order.
                    intVals = sort(findGridIntercepts(lambda(:,ii),data(:,ii),minVal(ii)+(maxVal(ii)-minVal(ii))/2));

                    %Use the maximum point as an initial guess
                    pumpNm(ii) = lambda(maxInd(ii),ii);

                    %If two or more points are available, update initial guess
                    %as average of fwhm values
                    if length(intVals) > 1
                        lowPts = intVals(intVals<pumpNm(ii));  %intercepts lower in nm value than the peak
                        highPts = intVals(intVals>pumpNm(ii)); %intercepts higher in nm value than the peak
                        if ~isempty(lowPts) && ~isempty(highPts)  %if both a higher and lower point exist
                            pumpNm(ii) = (lowPts(end) + highPts(1))/2; %take the average to be the peak center
                        end
                    end
                end %dropped grating position
            end %for loop
            
            %average the results over non-dropped grating positions
            pumpNm = mean(pumpNm(~isnan(pumpNm)));
            
            %make sure pumpNm is not NaN or empty before assigning value
            assert(~isempty(pumpNm) && ~isnan(pumpNm),['Failed to find center wavelength for' obj.desc.name]);
            
            %update the raman pump wavelength
            obj.ramanPumpNm = pumpNm;
            
            %set units back to input units
            obj = obj.setUnits(tmpUnits{:});
        end
        
    end
end