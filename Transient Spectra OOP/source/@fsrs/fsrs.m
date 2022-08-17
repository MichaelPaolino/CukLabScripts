classdef fsrs < transientSpectra
    properties
        ramanPumpNm = 400;  %(double) raman pump wavelength in nm
        %implement conversions methods
    end

    % Constructor, load, get and set methods with custom implementation
    methods
    %%CONSTRUCTOR/LOAD METHODS%%

        function obj = fsrs(varargin)
        % A FSRS object contains FSRS-specific funtionality in addition to
        % everything inside the transientSpectra class. The main additional feature
        % is defining a raman pump wavelength and the raman pump unit.
        % 
        % obj = fsrs(...);
        %   Constructs a fsrs object with the same arguments as the transientSpectra
        %   constructor call.
        %
        % obj = fsrs(..., 'ramanPumpNm', nmVal);
        %   Constructs a fsrs object and additionally assigns the ramanPumpNm
        %   value. The default value is 400 nm.
        %
        % FSRS specific units are:
        %   Spectra: (%Gain) Raman gain in % or (ppmGain) in ppm
        %   Delays: none
        %   Wavelengths: (rcm-1) Raman shift away from the fundamental in cm-1
        %
        % The default units for a FSRS object are mOD, ps, and rcm-1
        %
        % See also: TRANSIENTSPECTRA, DOUBLEWITHUNITS
            
        %%--SPLIT FSRS SPECIFIC AND TRANSIENT SPECTRA SPECIFIC ARGS--%%
            %fsrs specific arguments
            %cell array of fsrs specific keywords
            fsrsKeywords = {'ramanPumpNm'};
            %double array of corresponding input argument length for each
            %keyword (1 for keyword, 2 for name-value pair, etc.)
            fsrsArgLen = [2];
            
            %Loop over varargin and if the arguemnt keyword matches a fsrs keyword:
            %1. set it as empty
            %2. add it to the fsrs args
            %3. update iterator by expected length of keyword/name-value pair
            superClassArgs = varargin;  %these will be sent to the transientSpectra constructor
            subClassArgs = cell(size(varargin));   %these will be parsed inside the fsrs constructor
            argInd = 1;
            
            while argInd <= nargin  %loop until the end of varargin is reached
               if ischar(superClassArgs{argInd})    %ensure the current argument is of type char (i.e. a keyword and not numeric/cell)
                   cmpArgInd = strcmp(superClassArgs{argInd},fsrsKeywords);
                   if any(cmpArgInd)  %if the argument matches any of the fsrs keywords
                       %Ensure that the passed argument name is unique
                       assert(sum(cmpArgInd)==1,'Multiple fsrs arguments with the same name are not allowed.');
                       
                       %Distribute arguments to sub and superclass arg cell arrays
                       subClassArgs(argInd:argInd+fsrsArgLen(cmpArgInd)-1) = superClassArgs(argInd:argInd+fsrsArgLen(cmpArgInd)-1);   %copy fsrs args into fsrs arg list
                       superClassArgs(argInd:argInd+fsrsArgLen(cmpArgInd)-1) = {[]}; %set the transientSpectra keyword and antecedent values to empty (flag for removal)
                       argInd = argInd + fsrsArgLen;    %update iterator by argument length
                   else %this is a superclass argument
                       argInd = argInd + 1;    %update iterator by 1
                   end
               else %this is a superclass argument
                   argInd = argInd + 1;    %update iterator by 1
               end
            end
            
            % remove empty cell elements
            subClassArgs = subClassArgs(~cellfun(@isempty,subClassArgs));
            superClassArgs = superClassArgs(~cellfun(@isempty,superClassArgs));
            nSubClassArgs = numel(subClassArgs);
            
        %%--BUILD SUPERCLASS OBJECT--%%
            obj = obj@transientSpectra(superClassArgs{:}); 
            
        %%--PARSE SUBCLASS ARGS--%%
           %Default options for fsrs
           doSplitSchemes = false;
        
           %To support object arrays, inputs classes will be queried as
           %elements of cell arrays. If the input is not a cell array,
           %convert it to a cell array first.
           for ii = 1:nSubClassArgs %loop over fsrs args
               if ~iscell(subClassArgs{ii}) %if input is not already a cell
                   subClassArgs{ii} = {subClassArgs{ii}}; %convert non-cell varargin elements to cell
               end
           end
           
           %get object size to be able to assign values to each element individually
           objSize = size(obj); %original object size
           argNumel = numel(obj);   %number of object elements (and cell length of expected inputs)
           obj = obj(:);    %for easy looping/assignemnt
           
           %loop over arguments and parse keywords/name-value pairs
           argInd = 1;
           while argInd <= nSubClassArgs
               %each keyword or name-value pair must start with a char argument
               assert(ischar(subClassArgs{argInd}{1}),['Expected element or cell array of chars for ',...
                      'keywords or name-value pairs. Got ' class(subClassArgs{argInd}{1}) '.']);
               
               %parse the keyword in the current index
               switch subClassArgs{argInd}{1}
                   case 'ramanPumpNm'   %assign the raman pump wavelength
                       %perform checks on value input type
                       assert(argInd+1<=nSubClassArgs,'Name-value pair ramanPumpNm requires an additional double input');
                       assert(isscalar(subClassArgs{argInd+1}{1}),'Name-value pair ramanPumpNm requires an additional',...
                              ' scalar or cell array of scalar input.');  
                       assert(any(numel(subClassArgs{argInd+1})==[1,argNumel]),['The number of ramanPumpNm must be 1 or',...
                              'match the number of objects. Expected 1 or ' num2str(argNumel) ' elements, got ',...
                               num2str(numel(subClassArgs{argInd+1})) ' elements.']); 
                           
                       %assign array of arguments directly using deal. This works for both 1 or several pump wavelengths
                       [obj(:).ramanPumpNm] = deal(subClassArgs{argInd+1}{:});
                       argInd = argInd + 2;
                   otherwise %throw error to avoid infinite loop
                       error([subClassArgs{argInd}{1} ' is an unsupported keyword or name-value pair.']);
               end
           end
           
        %%--FINAL OBJECT FORMATTING--%%
           % reshape object to match input array size
           obj = reshape(obj, objSize);
        end
                
        %%GET/SET METHODS%% 
        
        %Set the raman pump nm. This updates the raman shift unit
        %definition, which is why it requires an explicit set method.
        function obj = set.ramanPumpNm(obj,newNm)
        % Update ramanPumpNm to a new value. 
        % 
        % obj(ind).ramanPumpNm = newVal;
        %   Updates ramanPumpNm and recalculates the rcm-1 unit
            
            obj.ramanPumpNm = newNm;
            
            %todo: figure out how to do this for all properties that might
            %have unit changes. Maybe cosmetic struct has .wavelengths,
            %.delays and .spectra properties
            obj.wavelengths = obj.wavelengths.updateRule('rcm-1','Raman Shift (cm^{-1})',...
                @(f) 1e7*(1/obj.ramanPumpNm-1./f),...
                @(f) 1./(1/obj.ramanPumpNm-1e-7*f));
        end
        
        function obj = setRamanPump(obj, newNm)
        % SETRAMANPUMP sets the raman pump wavelength in nm for each element in
        % obj. Use this method when obj is an array of fsrs objects. Just like the
        % accessor method, this method updates cm-1 unit definition as well.
        % 
        % obj = obj.SETRAMANPUMP(newNm)
        %   Updates the raman pump wavelength in nm for each element in obj
        % 
        % todo: allow for newNm to be a fsrs object or array
        % See Also: findRamanPumpNm, calibrateRamanPump

            % Get object size and convert to column
            objSize = size(obj);
            objNumel = numel(obj);
            obj = obj(:);
            
            % Set ramanPumpNm by element (could also use deal?)
            for objInd = 1:objNumel
               obj(onjInd).ramanPumpNm = newNm;
            end
            
            %convert object array back to original size
            obj = reshape(obj,objSize);
        end
        
    end
    
    methods (Access = protected)
        %override superclass convertDH method to convert multi-scheme data
        %into specific scheme objects. Call this to load multiple objects
        function obj = convertDH(obj, dh_static, dh_array)
        % CONVERTDH for FSRS class calls parent class TRANSIENTSPECTRA CONVERTDH 
        % and assigns FSRS specific units to the spectra, wavelengths, and delays.
        %
        % FSRS specific units are:
        %   Spectra: (%Gain) Raman gain in % or (ppmGain) in ppm
        %   Delays: none
        %   Wavelengths: (rcm-1) Raman shift away from the fundamental in cm-1
        %
        % The default units for a FSRS object are mOD, ps, and rcm-1
        %
        % See also: DOUBLEWITHUNITS
            
            %call parent method first for generic conversion of dh_array to transientSpectra
            obj = convertDH@transientSpectra(obj, dh_static, dh_array);
            
            %add unit rules to spectra, delays, etc. todo: add to method or as a subclass of doubleWithUnits class
            %spectra
            spectraRules = doubleWithUnits([],obj.spectra);
            spectraRules = spectraRules.addRule('%Gain','Raman Gain (%)',@(f) 1e2*(10.^f-1), @(f) log10(1+1e-2*f));
            spectraRules = spectraRules.addRule('ppmGain','Raman Gain (ppm)',@(f) 1e6*(10.^f-1), @(f) log10(1+1e-6*f));
            
            %delays (currently does nothing)
            delayRules = doubleWithUnits([],obj.delays);
            
            %wavelengths
            wlRules = doubleWithUnits([],obj.wavelengths);            
            wlRules = wlRules.addRule('rcm-1','Raman Shift (cm^{-1})',...
                @(f) 1e7*(1/obj.ramanPumpNm-1./f),...
                @(f) 1./(1/obj.ramanPumpNm-1e-7*f));
            
            %update object data to include new units. In addition, force
            %data to be real-valued
            obj.spectra = doubleWithUnits(real(obj.spectra.data),spectraRules);
            obj.spectra_std = doubleWithUnits(real(obj.spectra_std.data),spectraRules);
            obj.delays = doubleWithUnits(real(obj.delays.data),delayRules);
            obj.wavelengths = doubleWithUnits(real(obj.wavelengths.data),wlRules);
            
            %set default units
            obj = obj.setUnits('rcm-1','ps','mOD');
        end
    end
    
    % Methods specific to the fsrs class that cannot be implemented in the transientSpectra class
    methods
        
        function [obj, pumpNm] = findRamanPumpNm(obj, varargin)
        % FINDRAMANPUMPNM automatically finds the raman pump wavelength in nm. This
        % is useful for a rough calibration fo the raman shift for quick data
        % preview. An internal standard is required for more accurate calibration.
        %
        % The algorithm works by searching a wavelength region defined by guess +/- 
        % threshold and searching for the maximum point. Because the fundamental
        % peak is an image of the pump scatter on the spectrometer slit, it does 
        % not have a clean shape. Instead of fitting the peak, the algroithm 
        % attempts to determine the peak center by its FWHM. If the FWHM cannot be 
        % resolved, the algorithm takes the maximum point as the center. The final 
        % pump wavelength is the average over all valid grating positions.
        %
        % By default, the guess wavelength is obj.ramanPumpNm and the threshold is
        % 10 nm. These defaults can be overriden. While these values and the search
        % is performed in nm units, but the output object will have the same units 
        % as the input.
        %
        % Note: This automatic callibration is not very accurate. For accurate
        % calibration you should manually shift/scale the wavelength axis against a
        % known standard, such as a sulfate peak. In addition, the fundamental's
        % "Rowland Ghost" sidebands can be used for a more accurate 0 cm-1
        % calibration.
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
            pumpGuessNm = obj(1).ramanPumpNm;    %initial guess for search in nm
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
            
            % Format object array dims into a column for easy looping
            objSize = size(obj);
            objNumel = numel(obj);
            obj = obj(:);
            
            for objInd = 1:objNumel
                %remember old units
                tmpUnits = cell(3,1);
                [tmpUnits{:}] = obj(objInd).getUnits();

                %update units to units where x-axis is in nm and raman pump scatter is positive
                obj(objInd) = obj(objInd).setUnits('nm',[],'%Gain');

                subObj = obj(objInd).trim('wavelengths',pumpGuessNm+thresholdNm*[-1,1]);

                %get a single spectra by averaging over any extra dims except grating position (up to 5)
                data = permute(mean(subObj.spectra.data,[2, 3, 5],'omitnan'),[1,4,2,3,5]); 

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
                assert(~isempty(pumpNm) && ~isnan(pumpNm),['Failed to find center wavelength for' obj(objInd).name]);

                %update the raman pump wavelength
                obj(objInd).ramanPumpNm = pumpNm;

                %set units back to input units
                obj(objInd) = obj(objInd).setUnits(tmpUnits{:});
            end
            
            %convert object array back to original size
            obj = reshape(obj,objSize);
        end
        
        function [obj, shiftRcm] = calibrateRamanPump(obj, pixelInd, pixelRcm)
            % Format object array dims into a column for easy looping
            objSize = size(obj);
            objNumel = numel(obj);
            pixelInd = pixelInd(:);
            pixelRcm = pixelRcm(:);
            
            %optional input to calibrate all objects in the same way
            if isscalar(pixelInd) && isscalar(pixelRcm)
                pixelInd = pixelInd*ones(1,objNumel);
                pixelRcm = pixelRcm*ones(1,objNumel);
            end
            
            shiftRcm = zeros(1,objNumel);
            obj = obj(:);
            
            for objInd = 1:objNumel
                %remember old units and set new units to rcm-1
                tmpUnits = cell(3,1);
                [tmpUnits{:}] = obj(objInd).getUnits();
                obj(objInd) = obj(objInd).setUnits('rcm-1','','');
                
                %calculate the raman shift and set as double with units object
                shiftRcm(objInd) = mean(obj(objInd).wavelengths.data(pixelInd(objInd),:,:),'all')-pixelRcm(objInd);
                shiftNm = doubleWithUnits(shiftRcm(objInd),obj(objInd).wavelengths);
                shiftNm.unit = 'nm';
                
                obj(objInd).ramanPumpNm = shiftNm.data;
                
                %set units back to input units
                obj(objInd) = obj(objInd).setUnits(tmpUnits{:});
            end
            
            %convert object array back to original size
            obj = reshape(obj,objSize);
            shiftRcm = reshape(shiftRcm,objSize);
        end
        
    end
end