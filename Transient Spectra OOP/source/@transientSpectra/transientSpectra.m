classdef transientSpectra < matlab.mixin.Heterogeneous
	properties
        %raw data
        spectra = doubleWithUnits();  %[pixels, delays, rpts, grating pos, schemes]
        spectra_std = doubleWithUnits(); %[pixels, delays, rpts, grating pos, schemes]
        wavelengths = doubleWithUnits();  %[pixels, grating pos]
        delays = doubleWithUnits();   %[delays, repeats, grating pos]
        gPos = zeros(); %[]
        
        %identifying information: todo (?): convert to calss
        desc = struct('name', '',...
                      'shortName', '', ...
                      'description', '');
                         
        schemes = {}; %list of data schemes
                         
        %cosmetic information for data display. 
        %todo: convert to class so that display method calls can auto parse
        %varargin and return an updated  cosmetics object that can 
        %auto-convert data and auto-update axis labels
        cosmetic = struct('pixelUnits','nm',...
                          'signalUnits','OD',...
                          'delayUnits','ps',...
                          'pixelLimits',[0,0],...
                          'delayLimits',[0,0],...
                          'targetScheme',1);
        
        baseUnits = struct('pixels','nm',...
                           'delay','ps',...
                           'signal','OD');
        
        %size information
        sizes = struct('nRpts', 0, ...
                       'nGPos', 0, ...
                       'nSchemes', 0, ...
                       'nDelays', 0, ...
                       'nPixels', 0);
    end
    
    %methods that children classes must implement:
%     methods (Abstract)
%         
%         %converts any inferior class to another inferior class 
%         convertTransient(obj, targetObj);
%     end
    
    %% Constructor, load, export, get and set methods
    methods

        %%**CONSTRUCTOR METHODS**%%
        function obj = transientSpectra(varargin)
        %inputs:
        %transientSpectra(); initialize default FSRS object
        %transientSpectra(path); load data holder object or structs from path and
        %   convert to 
        
           switch nargin    %number of input arguments for the constructor call
               case 0   %constructs a default object --do nothing
               case 1   %one argument is either a path or a data_holder object
                   argClass = class(varargin{1});   %determines input class
                   switch argClass
                       case 'char'  %char class, meaning a path
                           obj = loadPath(obj,varargin{1});
                       case 'data_holder'   %data_holder class
                           %call load routine for a data_holder
                       otherwise    %invalid input class
                           %return an error
                   end
                   
               case 2   %user manually input dh_static and dh_array
                   %expected input type is a pair of structs
                   if isstruct(varargin{1}) && isstruct(varargin{2})
                       %build data_holder and call load routine for a
                       %data_holder
                   else     %invalid type
                       %return an error
                   end
               otherwise
                   %return an error
           end
        end
       
        %loads a .mat file from a path and converts to a FSRS object
        function varargout = loadPath(obj,myPath)
            loaded = load(myPath);   %load the path contents 
            contents = fieldnames(loaded);   %get the variables in the loaded data
            fileName = split(myPath,'\');
            fileName = fileName{end};
            switch class(loaded.(contents{1}))   %check the variable class
                case 'struct'    %this should be a data_holder object
                    if strcmp(contents{1},'dh_static') %this is a raw data data_holder
                        %first load with acquisition data holder
                        [varargout{1:nargout}] = convertDH(obj,loaded.dh_static,loaded.dh_array);
                        
                        %update name in output objects
                        for ii = 1:nargout
                           varargout{ii}.desc.name = fileName; 
                        end
                    else
                        %return error
                    end
                case 'data_holder'
                   %this is a data_holder object
                   %call dh conversion routine
                case 'FSRS'
                   %this is a FSRS object
                otherwise
                   %error
            end
        end
        
        %convert a data holder into a transientSpectra. See convertDH.m for
        %generic implementation. Override this method in subclasses for specific
        %implementation
        obj = convertDH(obj, dh_static, dh_array);
        
        %export the object data to file
        [outputStruct, filePath] = export(obj,filePath,varargin)
        
        %%**GET-SET METHODS**%%
        %update units on all unit arrays
        function obj = setUnits(obj,wavelengthUnit,delayUnit,spectraUnit)
           if ~isempty(wavelengthUnit)
               obj.wavelengths.unit = wavelengthUnit;
           end
           
           if ~isempty(delayUnit)
               obj.delays.unit = delayUnit;
           end
           
           if ~isempty(spectraUnit)
               obj.spectra.unit = spectraUnit;
               obj.spectra_std.unit = spectraUnit;
           end
           
        end
        
        %return current units
        function [wavelengthUnit, delayUnit, spectraUnit] = getUnits(obj)
            wavelengthUnit = obj.wavelengths.unit;
            delayUnit = obj.delays.unit;
            spectraUnit = obj.spectra.unit;
        end
    end   
    
    %% Plotting/contour methods for data display
    methods
        function plotSpectra(obj, varargin)
            %%**INITIALIZE DEFAULT VALUES**%%
            %default for delay display: all delays
            delaysVal = mean(obj.delays.data,[2,3]);
            nDelays = length(delaysVal);
            delaysInd = 1:nDelays;
            showDelays = true;  %whether to display in legend
            
            %default for repeat display: average repeats
            rpts = [];
            nRpts = 1;
            showRepeats = false; %whether to display in legend
            
            %default for grating position display: all grating positions
            gPosVal = obj.gPos;
            nGPos = length(gPosVal);
            showGPos = false; %whether to display in legend
            
            %default for cosmetics
            showLegend = true;
            
            %todo: decide how to handle NaN points
            
            %%**UPDATE DEFAULT VALUES FROM USER INPUT**%%
            %defines additional arguments passed to plot function. There
            %will also be a cosmetics version of this inside the cosmetic
            %class.
            isArgParsed = false(size(varargin));
            
            %change settings from user input
            if nargin > 1
                for ii = 1:(nargin-1)   %loop over extra arguments
                    if ischar(varargin{ii}) %look for 'name' in name-value pair
                        switch varargin{ii} %switch...case over name if name is encountered
                            case 'delays'   %plot a subset of delays in data
                                %find unique values and indecies in data delays that best match user value input in name-value pair 
                                [delaysVal, delaysInd] = nearestVal(delaysVal, varargin{ii+1}(:));
                                nDelays = length(delaysVal);
                                
                                %flag the name-value pair as parsed to remove the args later below
                                isArgParsed(ii+[0 1]) = [true true];
                            case 'no legend'
                                showLegend = false;
                                
                                %flag the name as parsed to remove the args later below
                                isArgParsed(ii) = true;
                        end %switch varargin{ii}
                    end %ischar(varargin{ii})
                end %ii = 1:(nargin-1)
            end % if nargin > 1
            
            %Update extra arguments
            extraArgs = varargin(~isArgParsed);   %add unparsed arguments to extraargs
            
                     
            %generate legend
            delayStr = strtrim(cellstr([num2str(delaysVal(:)) repmat(' ps',nDelays,1)]));   %todo: replace ps with cosmetic unit and add option to specificy precision
            %rptStr = ... %todo: finish rpt formatting
            %gPosStr = ... %todo: finsih gPos formatting
            legendVal = delayStr(:);    %todo: somehow build a legend string depending on user prefs above
            
            %check hold state
            holdState = ishold();
            if ~holdState
                hold on;
            end
            
            %%**GENERATE PLOT**%%
            for ii = 1:obj.sizes.nGPos
                %add data to x, y plot
                x = obj.wavelengths.data(:,ii);
                %avg of rpts and pick the first dataScheme. todo: have case that handles whether to average rpts or not
                y = mean(obj.spectra.data(:,delaysInd,:,ii),3); %[pixels, delays, rpts, grating pos, schemes]
                %y = permute(y,[]); %change display priority

                plotArgs = [{x; y}; extraArgs{:}];  %custom generate inputs to pass to plot function
                plot(plotArgs{:});
            end
            
            %return to previous hold state
            if ~holdState
                hold off;
            end
            
            %decide on legend formattiong
            if showLegend
               legend(legendVal(:)); 
               %todo: add multi-d legend display
            end
            
            ylabel(obj.spectra.dispName);
            xlabel(obj.wavelengths.dispName);
            box on;
        end
        
        function plotTrace(obj, varargin)
            %%**INITIALIZE DEFAULT VALUES**%%
            %default for delay display: all wavelengths or wavenumbers
            xVal = obj.wavelengths.data;
            nX = length(xVal);
            waveInd = 1:nX;
            showX = true;  %whether to display in legend
            
            %default for repeat display: average repeats
            rpts = [];
            nRpts = 1;
            showRepeats = false; %whether to display in legend
            
            %default for grating position display: all grating positions
            gPosVal = obj.gPos;
            nGPos = length(gPosVal);
            showGPos = false; %whether to display in legend
            
            %default for cosmetics
            showLegend = true;
                        
            
            %todo: decide how to handle NaN points
            
            %%**UPDATE DEFAULT VALUES FROM USER INPUT**%%
            %defines additional arguments passed to plot function. There
            %will also be a cosmetics version of this inside the cosmetic
            %class.
            isArgParsed = false(size(varargin));
            
            %change settings from user input
            if nargin > 1
                for ii = 1:(nargin-1)   %loop over extra arguments
                    if ischar(varargin{ii}) %look for 'name' in name-value pair
                        switch varargin{ii} %switch...case over name if name is encountered
                            case 'wavelengths'   %plot a subset of delays in data
                                %find unique values and indecies in data delays that best match user value input in name-value pair 
                                [xVal, waveInd] = nearestVal(xVal, varargin{ii+1}(:));
                                nX = length(xVal);
                                
                                %flag the name-value pair as parsed to remove the args later below
                                isArgParsed(ii+[0 1]) = [true true];
                            case 'no legend'
                                showLegend = false;
                                
                                %flag the name as parsed to remove the args later below
                                isArgParsed(ii) = true;
                        end %switch varargin{ii}
                    end %ischar(varargin{ii})
                end %ii = 1:(nargin-1)
            end % if nargin > 1
            
            %Update extra arguments
            extraArgs = varargin(~isArgParsed);   %add unparsed arguments to extraargs
            
            %%**GENERATE PLOT**%%
            %add data to x, y plot
            x = obj.delays;
            %avg of rpts and pick the first dataScheme. todo: have case that handles whether to average rpts or not
            y = 1000*mean(obj.spectra(waveInd,:,:,:,1),4); %[pixels, delays, rpts, grating pos, schemes]
            y = permute(y,[3,2,1,4,5]); %change display priority [delays, schemes, pixels, rpts, grating pos]. todo: figure out if disply priority needs to be different for kinetic traces
            y = reshape(y,obj.sizes.nDelays,[]);    %turns y into a 2d array [delays, everything else] where left most index is most significant
            
            %generate legend
            waveStr = strtrim(cellstr([num2str(xVal(:)) repmat(' nm', nX,1)]));   %todo: replace ps with cosmetic unit and add option to specificy precision
            %rptStr = ... %todo: finish rpt formatting
            %gPosStr = ... %todo: finsih gPos formatting
            legendVal = waveStr(:);    %todo: somehow build a legend string depending on user prefs above
            
            plotArgs = [{x; y}; extraArgs{:}];  %custom generate inputs to pass to plot function
            plot(plotArgs{:});
            
            %decide on legend formattiong
            if showLegend
               legend(legendVal(:)); 
               %todo: add multi-d legend display
            end
            ylabel(obj.spectra.dispName);
            xlabel(obj.delays.dispName);
        end
    end
    
    %% Data manipulation methods that modify the spectra
    methods
        %stitches all grating positions together
        function obj = stitch(obj,varargin)
        % Stitch together available grating positions in object. This works by
        % sorting the data in ascending wavelength and stitching a pair of grating
        % positions at a time. For multiple grating positions, the previously
        % stitched grating position is treated as the first grating position in the
        % stich pair. The stitching behavior in any overlap region can be done with
        % the following strageiges: average, lower, upper, half, and linear. 
        % Stitching preserves all other dimensions, such as repeats and delays.
        % Stitching currently updates all relevant properties such as wavelengths,
        % spectra, and grating positions. todo: update spectra_std. todo: may need
        % to create a case for direct stitching
        % Inputs (varargin):
        %   strategy: Optional char array. Choose from 'average', 'lower', 'upper',
        %       'half', and 'linear'. The default is 'linear'
        % Outputs:
        %   obj: object with all grating positions stitched into one and with
        %       relevat properties updated
            
            %optional inputs: strategy for stitching. Default is linear
            if nargin==2
                strategy = varargin{1};
            else
                strategy = 'linear';
            end
            
        %%--Revert units to nm--%%
            %remember old units
            tmpUnits = cell(3,1);
            [tmpUnits{:}] = obj.getUnits();
            
            %update units to units where x-axis is in nm
            obj = obj.setUnits('nm',[],[]);
            
        %%--Sort wavelengths and grating positions--%%
            %sort everything by ascending order in terms of grating position nm, wavelength nm, delay ps
            [gPosSorted,gInd] = sort(obj.gPos);
            [wl, lInd] = sort(obj.wavelengths.data(:,gInd));    %[wavelengths, grating positions]
            t = obj.delays.data(:,:,gInd);    %[delays, repeats, grating positions]
            
            %sort by grating position first
            data = obj.spectra.data(:,:,:,gInd,:);  %[pixels, delays, rpts, GRATING POS, schemes]
                        
            %sort by wavelengths next
            for ii = 1:size(lInd,2) %loop over grating positions
                data(:,:,:,ii,:) = data(lInd(:,ii),:,:,ii,:); %[PIXELS, delays, rpts, GRATING POS, schemes]
            end
            
            %permute and reshape data so that it is easier to interpolate and stitch
            data = permute(data,[1,2,3,5,4]);   %[pixels, delays, rpts, schemes, grating pos]
            data = reshape(data,obj.sizes.nPixels,[],obj.sizes.nGPos); %[pixels, delays*rpts*schemes, grating pos]
            
        %%--Stitch wavelengths and data--%%
            %grating positions will be stitched in pairs sequentially. If
            %there are multiple grating positions, already stitched data
            %will be treated as one grating position
            tmpWl = wl(:,1);    %holder for concatonated wavelengths starting from 1st grating pos
            tmpData = data(:,:,1);  %holder for concatonated data starting from 1st grating pos
            
            for ii = 2:obj.sizes.nGPos %loop over grating positions starting from the 2nd one
                
                %First find overlap region in wavelengths
                high1 = tmpWl(end);    %highest wavelength in 1st grating position
                low2 = wl(1,ii);     %lowest wavelength in 2nd grating position
                
                [low1, low1Ind] = nearestVal(tmpWl,low2); %lowest wavelength and index in 1st grating position
                [high2, high2Ind] = nearestVal(wl(:,ii),high1); %lowest wavelength and index in 2nd grating position
                
                %Make sure indicies do not cause data/wavelength extrapolation (overlap region needs to be inclusive)
                if low1<low2   %check first grating position
                    low1Ind = low1Ind + 1;
                    low1 = tmpWl(low1Ind);
                end
                
                if high2>high1 %check second grating position
                    high2Ind = high2Ind - 1;
                    high2 = wl(high2Ind,ii);
                end
                
                %define lower and upper regions
                lowWl = tmpWl(1:low1Ind-1); %1st grating positions lower wavelengths
                highWl = wl(high2Ind+1:end,ii); %2nd grating positions higher wavelengths
                lowData = tmpData(1:low1Ind-1,:);   %1st grating positions lower data
                highData = data(high2Ind+1:end,:,ii);   %second grating positoins higher data
                
                %Next interpolate middle region
                nInd = ceil(0.5*(length(tmpWl(low1Ind:end))+length(wl(1:high2Ind,ii))));    %average number of in-between indicies
                midWl = linspace(low1,high2,nInd)';  %new wavelengths in overlap region with linear spacing between them
                mid1Data = interp1(tmpWl,tmpData,midWl,'linear');   %interpolate 1st grating position on overlap scale
                mid2Data = interp1(wl(:,ii),data(:,:,ii),midWl,'linear'); %interpolate 2nd grating position on overlap scale
                
                %Execute strategy to combine data in overlap region
                switch strategy
                    case 'average' %take the average in the overlap region
                        midData = 0.5*(mid1Data+mid2Data);  
                    case 'lower' %take the lower (1st) grating position
                        midData = mid1Data;  
                    case 'upper' %take the higher (2nd) grating position
                        midData = mid2Data;  
                    case 'half'  %stitch at half-way point
                        midData = [mid1Data(1:floor(nInd/2),:); mid2Data((floor(nInd/2)+1):end,:)];
                    case 'linear'   %do a weighted average with a linear sweep of the weights from the 1st to 2nd grating position
                        weightsFun = polyfit([midWl(1), midWl(end)],[0, 1],1);    %linear fit for weights
                        weights = polyval(weightsFun,midWl);    %calculate weights from nm values
                        midData = mid1Data.*(1-weights)+mid2Data.*weights;  %this does the weighted average
                    otherwise
                        error(['Unsupported strategy. Available stragegies are average, lower, upper, half, and linear. Got ' strategy '.']);
                end
                %concatonate the three regions across the 1st (pixel) dimension
                tmpWl = [lowWl; midWl; highWl]; %store in tmpWl so that it can be fed as 1st grating position in next iteration 
                tmpData = [lowData; midData; highData]; %store in tmpdata so that it can be fed as 1st grating position in next iteration 
            end
            
        %%--Convert wavelengths and data back to original state and update object parameters--%%
            %reshape and permute data back to original dim order
            data = reshape(tmpData,[],obj.sizes.nDelays,obj.sizes.nRpts,obj.sizes.nSchemes); %[pixels, delays, rpts, schemes, grating pos]
            data = permute(data,[1,2,3,5,4]); %[pixels, delays, rpts, grating pos, schemes]
            
            %update values for object properties
            obj.spectra.data = data;
            obj.wavelengths.data = tmpWl;
            obj.gPos = median(gPosSorted);
            obj.sizes.nPixels = length(tmpWl);
            obj.sizes.nGPos = 1;
            
            %set units back to input units
            obj = obj.setUnits(tmpUnits{:});
        end
        
        function obj = trim(obj, varargin)
            trim = struct('wls','all',...
                          'delays','all');
            
            %parse varargin
            if nargin > 1
                for ii = 1:2:(nargin-1)
                    assert(ischar(varargin{ii}),...
                        ['Invalid argument class for name-value pair. Expected class char for name, got ' class(varargin{ii}) '.']);
                    switch varargin{ii}
                        case 'wavelengths'
                            trim.wls = varargin{ii+1};
                        case 'delays'
                            trim.delays = varargin{ii+1};
                        otherwise
                            error([varargin{ii} ' is not a valid argument name.']); 
                    end
                end
            end
                      
            %trim wavelengths
            if ischar(trim.wls)  %this is a do nothing case
                % Assert correct input to ensure user isn't accidently doing something they're not aware of
                assert(strcmp(trim.wls,'all'),'Expected all keyword or wavelength range [wl1, wl2] of type double.');
                
            elseif isa(trim.wls,'double')   %this does the wavelength trim
                %ensure the user has correct input before trimming
                assert(length(trim.wls)==2, 'Expected wavelength range [wl1, wl2] of type double.');
                
                %find the wavelength range indicies
                wls = obj.wavelengths.data;  %[pixels, gPos]
                [~,wlInd] = nearestVal(wls,trim.wls,... %wlInd is [2,gPos]
                                'unique',false,'trim',false); %keep all values 
                        
                %for easy looping, do the following dim rearrangement:
                %[pixels, delays, rpts, gpos, schemes] -> [pixels, delays x rpts x schemes, gpos]
                tmpSpectra = permute(obj.spectra.data,[1,2,3,5,4]);
                tmpSpectra = reshape(tmpSpectra,obj.sizes.nPixels,[],obj.sizes.nGPos);

                %Allocate NaN double arrays and place spectra values into it 
                trimmedSpectra = nan(size(tmpSpectra)); %[pixels, delays x rpts x schemes, gpos]
                trimmedWls = nan(size(wls)); %[pixels, gpos]
                trimmedGPos = nan(size(obj.gPos));  %[1,gpos]

                %Flag remove grating positions that do not at all span the desired wavelength range.
                %There are two things to check here:
                %1. Is at least one target index is in range of the grating position?
                %2. Are the found indicies on the same, indicating the target range is above or below the grating position range
                %If 1 is true or 2 is false, then keep the grating position
                isIndInRange = any(abs(wls-permute(trim.wls(:),[3,2,1]))<0.01*(max(wls)-min(wls)),[1,3]);
                keepGPos = or(wlInd(1,:)~=wlInd(2,:),isIndInRange);
                
                %loop over grating positions to sub-select range dicated by wlInd
                for ii = find(keepGPos)
                    %copy desired subrange for each grating position into the NaN arrays starting from index 1
                    trimmedWls(1:length(wlInd(1,ii):wlInd(2,ii)),ii) = wls(wlInd(1,ii):wlInd(2,ii),ii);
                    trimmedSpectra(1:length(wlInd(1,ii):wlInd(2,ii)),:,ii) = tmpSpectra(wlInd(1,ii):wlInd(2,ii),:,ii); 

                    %remove NaN flag from new gPos only if the grating position has data in it
                    trimmedGPos(ii) = obj.gPos(ii);
                end
                
                %remove any dims that are all NaN
                trimmedWls = trimmedWls(any(~isnan(trimmedWls),2),~isnan(trimmedGPos)); %[wls, gpos]
                trimmedSpectra = trimmedSpectra(any(~isnan(trimmedWls),2),:,~isnan(trimmedGPos)); %[pixels, delays x rpts x schemes, gpos]
                trimmedGPos = trimmedGPos(~isnan(trimmedGPos));
                
                %update sizes
                obj.sizes.nPixels = size(trimmedWls,1);
                obj.sizes.nGPos = length(trimmedGPos);
                
                %convert spectra back to original dimensions and dim order
                %[pixels, delays x rpts x schemes, gpos] -> [pixels, delays, rpts, gpos, schemes]
                trimmedSpectra = reshape(trimmedSpectra,obj.sizes.nPixels,obj.sizes.nDelays,obj.sizes.nRpts,obj.sizes.nSchemes,obj.sizes.nGPos); %[pixels, delays, rpts, schemes, gpos]
                trimmedSpectra = permute(trimmedSpectra,[1,2,3,5,4]); %[pixels, delays, rpts, gpos, schemes]
                
                %add data back to object
                obj.wavelengths.data = trimmedWls;
                obj.spectra.data = trimmedSpectra;
                obj.gPos = trimmedGPos;
                
            else
                error('Expected all keyword or wavelength range [wl1, wl2] of type double.');
            end
            
            %trim delays
            if ischar(trim.delays)  %this is a do nothing case
                % Assert correct input to ensure user isn't accidently doing something they're not aware of
                assert(strcmp(trim.delays,'all'),'Expected all keyword or delay range [d1, d2] of type double.');
                
            elseif isa(trim.delays,'double')   %this does the wavelength trim
                %ensure the user has correct input before trimming
                assert(length(trim.delays)==2, 'Expected delay range [d1, d2] of type double.');
                
                %find the wavelength range indicies
                t = reshape(obj.delays.data, obj.sizes.nDelays, []);  %[delays, rpts x gPos]
                [~,tInd] = nearestVal(t,trim.delays,... %tInd is [2,rpts x gPos]
                                'unique',false,'trim',false); %keep all values 
                            
                %for easy looping, do the following dim rearrangement:
                %[pixels, delays, rpts, gpos, schemes] -> [delays, pixels x schemes, rpts x gpos]
                tmpSpectra = permute(obj.spectra.data,[2,1,5,3,4]);
                tmpSpectra = reshape(tmpSpectra,obj.sizes.nDelays,obj.sizes.nPixels*obj.sizes.nSchemes,[]);

                %Allocate NaN double arrays and place spectra values into it 
                trimmedSpectra = nan(size(tmpSpectra)); %[delays, pixels x schems, rpts x gpos]
                trimmedDelays = nan(size(t)); %[delays, rpts x gpos]
                
                %Unlike wavelengths, the delays should all be similar for all grating positions and repeats,
                %so there is no flag remove of grating positions or repeats
                
                %loop over grating positions and repeats to sub-select range dicated by wlInd
                for ii = 1:size(tInd,2)
                    %copy desired subrange for each grating position into the NaN arrays starting from index 1
                    trimmedDelays(1:length(tInd(1,ii):tInd(2,ii)),ii) = t(tInd(1,ii):tInd(2,ii),ii); %[delays, rpts x gpos]
                    trimmedSpectra(1:length(tInd(1,ii):tInd(2,ii)),:,ii) = tmpSpectra(tInd(1,ii):tInd(2,ii),:,ii); %[delays, pixels x schems, rpts x gpos]
                end
                
                %remove any dims that are all NaN
                trimmedDelays = trimmedDelays(any(~isnan(trimmedDelays),2),:); %[delays, rpts x gpos]
                trimmedSpectra = trimmedSpectra(any(~isnan(trimmedDelays),2),:,:); %[delays, pixels x schems, rpts x gpos]
                
                %update sizes
                obj.sizes.nDelays = size(trimmedDelays,1);
                
                %convert spectra back to original dimensions and dim order
                %[delays, pixels x schemes, rpts x gpos] -> %[pixels, delays, rpts, gpos, schemes]
                trimmedDelays = reshape(trimmedDelays,obj.sizes.nDelays,obj.sizes.nRpts,obj.sizes.nGPos); %[delays, rpts, gPos]
                trimmedSpectra = reshape(trimmedSpectra,obj.sizes.nDelays,obj.sizes.nPixels,obj.sizes.nSchemes,obj.sizes.nRpts,obj.sizes.nGPos); %[delays, pixels, schemes, rpts, gpos]
                trimmedSpectra = permute(trimmedSpectra,[2,1,4,5,3]); %[pixels, delays, rpts, gpos, schemes]
                
                %add data back to object
                obj.delays.data = trimmedDelays;
                obj.spectra.data = trimmedSpectra;
                
            else
                error('Expected all keyword or delay range [d1, d2] of type double.');
            end
        end
    end
    
    
end