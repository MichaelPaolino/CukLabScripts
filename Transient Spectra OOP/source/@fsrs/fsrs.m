classdef fsrs < transientSpectra
    properties
        ramanPumpNm = 400;  %(double) raman pump wavelength in nm
        %implement conversions methods
    end

    % Constructor, load, get and set methods with custom implementation
    methods
        function obj = fsrs(varargin)
        % A FSRS object contains FSRS-specific funtionality in addition to
        % everything inside the transientSpectra class. The main additional feature
        % is defining a raman pump wavelength and the raman pump unit, along with 
        % supporting methods to automate processing.
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
            
            % First parse input args to determine what should be passed to the superclass constructor
            if nargin == 0 %constructs default fsrs object
                superClassArgs = {};
            else %construct fsrs object from data
                % Use input parser to parse user arguments
                p = inputParser;
                p.FunctionName = 'fsrs';
                p.StructExpand = false;
                p.KeepUnmatched = true; %unmatched inputs will be passed to superclass

                % Use valVarCell to validate both cell and non-cell inputs
                % dataSource validation: handled in superclass
                p.addRequired('dataSource');
                % ramanPumpNm validation: must be a scalar number
                p.addParameter('ramanPumpNm','', @(p) valVarCell(p,@(c) isscalar(c) && ~ischar(c)));

                % parse inputs and store results in struct p.Results and unmatched inputs in p.Unmatched
                p.parse(varargin{:});

                % Build superclass arguments as cell array from unmatched results
                fn = fieldnames(p.Unmatched);   %cell array of field names
                fVal = struct2cell(p.Unmatched);    %cell array of value names
                superClassArgs = [fn(:), fVal(:)]'; %cell array of name-value pairs [2, nArgs]
                superClassArgs = [{p.Results.dataSource}; superClassArgs(:)]; %cell array of dataSource, name-value pairs [dataSource; 2*nArgs]                
            end
            
            % Pass unmatched arguments to superclass and build object
            obj@transientSpectra(superClassArgs{:}); 
            
            % Update constructed object with fsrs-specific arguments
            if nargin > 0
                % Parse subclass args
                % Convert object array to 1D array for easy looping
                objSize = size(obj);    %remember the original object size
                objNumel = numel(obj);  %for easy looping, loop over elements of arguments
                obj = obj(:); %convert object to vector array

                % convert parsed results into struct whose elements are cells
                results = ensureCellVals(p.Results);
                % ensure that the cell array size matches the object size by and convert into column vector
                results.ramanPumpNm = explicitExpand(results.ramanPumpNm, objSize);
                results.ramanPumpNm = results.ramanPumpNm(:);

                % Loop over object array elements
                for objInd = 1:objNumel
                    % Update RamanPumpNm if it exists
                    if ~isempty(results.ramanPumpNm{objInd})
                        obj(objInd).ramanPumpNm = results.ramanPumpNm{objInd};
                    end
                end
               
                %set default units
                obj = obj.setUnits('rcm-1','ps','mOD');
                
                % reshape object to match input array size
                obj = reshape(obj, objSize);
            end
        end
                
        function obj = set.ramanPumpNm(obj,newNm)
        % Update ramanPumpNm to a new value. Note: this method only works for
        % single fsrs object elements. Use setRamanPump() for object array
        % functionality.
        % 
        % obj(ind).ramanPumpNm = newVal;
        %   Updates ramanPumpNm and recalculates the rcm-1 unit
        %
        % See Also: setRamanPump
            
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
        % The new wavelength can be input as either a scalar, double array, or fsrs
        % object array. For the array cases, this method supports expansion of
        % singleton dimensions to match the size of calling fsrs object.
        % 
        % obj = obj.SETRAMANPUMP(newNm)
        %   Updates the raman pump wavelength in nm to newNm for each element in
        %   obj. If newNm is a scalar, the same value is assigned to each element.
        %   If newNm is a fsrs object, the ramanPumpNm value is copied from newNm
        %   into obj. If newNm is an array, the array's singleton dimensions are 
        %   expanded to match the size of obj. 
        % 
        % See Also: findRamanPumpNm, calibrateRamanPump

            % Get object size and convert to column
            objSize = size(obj);
            objNumel = numel(obj);
            obj = obj(:);
            
            % if input is a fsrs object array, extract ramanPumpNm and
            % convert to matrix of the same size as newNm
            if isa(newNm,'fsrs')
                newNmSz = size(newNm);
                newNm = {newNm(:).ramanPumpNm};
                newNm = reshape(cell2mat(newNm),newNmSz);
            end
            
            % Expand singleton dims of newNm to match the size of objSize
            newNm = explicitExpand(newNm,objSize);
            
            % Set ramanPumpNm by element (could also use deal?)
            for objInd = 1:objNumel
               obj(objInd).ramanPumpNm = newNm(objInd);
            end
            
            %convert object array back to original size
            obj = reshape(obj,objSize);
        end
        
    end
    
    % Protected methods that define the inner workings of the class
    methods (Access = protected)
        %override superclass assignUnits method
        function [obj,unitRules] = assignUnits(obj)
        % for FSRS class calls parent class TRANSIENTSPECTRA ASSIGNUNITS 
        % and assigns FSRS specific units to the spectra and wavelengths.
        %
        % FSRS specific units are:
        %   Spectra: (%Gain) Raman gain in % or (ppmGain) in ppm
        %   Delays: none
        %   Wavelengths: (rcm-1) Raman shift away from the fundamental in cm-1
        %
        % See also: DOUBLEWITHUNITS
            
            %call parent method first for generic conversion of dh_array to transientSpectra
            [obj,unitRules] = assignUnits@transientSpectra(obj);
            
            %add unit rules to spectra, delays, etc. 
            spectraRules = doubleWithUnits([],unitRules.spectraRules);
            spectraRules = spectraRules.addRule('%Gain','Raman Gain (%)',@(f) 1e2*(10.^f-1), @(f) log10(1+1e-2*f));
            spectraRules = spectraRules.addRule('ppmGain','Raman Gain (ppm)',@(f) 1e6*(10.^f-1), @(f) log10(1+1e-6*f));
            
            %wavelengths
            wlRules = doubleWithUnits([],unitRules.wlRules);            
            wlRules = wlRules.addRule('rcm-1','Raman Shift (cm^{-1})',...
                @(f) 1e7*(1/obj.ramanPumpNm-1./f),...
                @(f) 1./(1/obj.ramanPumpNm-1e-7*f));
            
            % Format object array dims into a column for easy looping
            objSize = size(obj);
            objNumel = numel(obj);
            obj = obj(:);
            
            % Assign unit rules to object data for each object array element, while keeping existing data
            for objInd = 1:objNumel
                obj(objInd).spectra = doubleWithUnits(real(obj(objInd).spectra.data),spectraRules);
                obj(objInd).spectra_std = doubleWithUnits(real(obj(objInd).spectra_std.data),spectraRules);
                obj(objInd).wavelengths = doubleWithUnits(real(obj(objInd).wavelengths.data),wlRules); 
            end
            
            %convert object array back to original size
            obj = reshape(obj,objSize);
            
            %package unit rules to unitrules struct
            unitRules = struct('spectraRules',spectraRules,'delayRules',unitRules.delayRules,'wlRules',wlRules);
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