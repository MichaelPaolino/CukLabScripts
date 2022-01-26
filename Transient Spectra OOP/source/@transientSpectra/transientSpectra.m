%thoughts on heirarchy
%data
    %transientSpectra
    %       properties: nDTable spectra, nDTable spectra std, nDTable diagnostics (add later), delays,
    %       weavlenegths, plotSettings, displaySettings, exportSettings
    %       methods: loadPath, plotSpectra, plotTrace, contour
    %   probe
    %   TR
    %   FSRS
    %       properties: ramanPumpNm, wavenumber
    %       methods: loadPath, get.wavenumber, fitRamanPump
            %loadPath should return 2 FSRS objects, a TR object
    %   tATR
    %spectra
        %spectra1D
        %spectraTransient
        %spectraD
    %electrochem
        %CV
        %it
        
        %create a nDTable class
        %create a units class that has unit conversion, base unit, and
        %   function handle for conversion
        


classdef transientSpectra < matlab.mixin.Heterogeneous
	properties
        %raw data
        spectra = zeros();  %[pixels, schemes, delays, rpts, grating pos]
        spectra_std = zeros(); %[pixels, schemes, delays, rpts, grating pos]
        wavelengths = zeros();  %[]
        delays = zeros();   %[]
        gPos = zeros(); %[]
        
        %identifying information: todo (?): convert to calss
        description = struct('name', '',...
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
            switch class(loaded.(contents{1}))   %check the variable class
                case 'struct'    %this should be a data_holder object
                    if strcmp(contents{1},'dh_static') %this is a raw data data_holder
                        %first try to load with acquisition data holder
                        [varargout{1:nargout}] = convertDH(obj,loaded.dh_static,loaded.dh_array);
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
        convertDH(obj, dh_static, dh_array);
        
        %%**GET-SET METHODS**%%

        %%**DISPLAY DATA METHODS**%%
        function plotSpectra(obj, varargin)
            %%**INITIALIZE DEFAULT VALUES**%%
            %default for delay display: all delays
            delaysVal = obj.delays;
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
                                [delaysVal, delaysInd] = nearestVal(obj.delays, varargin{ii+1}(:));
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
            
            %%**GENERATE PLOT**%%
            %add data to x, y plot
            x = obj.wavelengths;
            %avg of rpts and pick the first dataScheme. todo: have case that handles whether to average rpts or not
            y = 1000*mean(obj.spectra(:,delaysInd,:,:,1),3); %[pixels, delays, rpts, grating pos, schemes]
            %y = permute(y,[]); %change display priority
            y = reshape(y,obj.sizes.nPixels,[]);    %turns y into a 2d array [pixels, everything else] where left most index is most significant
            
            %generate legend
            delayStr = strtrim(cellstr([num2str(delaysVal(:)) repmat(' ps',nDelays,1)]));   %todo: replace ps with cosmetic unit and add option to specificy precision
            %rptStr = ... %todo: finish rpt formatting
            %gPosStr = ... %todo: finsih gPos formatting
            legendVal = delayStr(:);    %todo: somehow build a legend string depending on user prefs above
            
            plotArgs = [{x; y}; extraArgs{:}];  %custom generate inputs to pass to plot function
            plot(plotArgs{:});
            
            %decide on legend formattiong
            if showLegend
               legend(legendVal(:)); 
               %todo: add multi-d legend display
            end
            ylabel('\DeltamOD');
            xlabel('nm');
        end
        
        function plotTrace(obj, varargin)
            %%**INITIALIZE DEFAULT VALUES**%%
            %default for delay display: all wavelengths or wavenumbers
            waveVal = obj.wavelengths;
            nWave = length(waveVal);
            waveInd = 1:nWave;
            showWave = true;  %whether to display in legend
            
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
                                [waveVal, waveInd] = nearestVal(obj.wavelengths, varargin{ii+1}(:));
                                nWave = length(waveVal);
                                
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
            waveStr = strtrim(cellstr([num2str(waveVal(:)) repmat(' nm', nWave,1)]));   %todo: replace ps with cosmetic unit and add option to specificy precision
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
            ylabel('\DeltamOD');
            xlabel('ps');
        end
    end   
end


%this should be a function outside of the class
function [valsOut, ind] = nearestVal(valsIn,targetVals)
%Find the nearest values targetVals in vector vals. This is useful if a 
%user wants to select a subset of delays or wavelengths from their data. 
%The function takes the raw delays or wavelengths and returns the nearest 
%values and indicies to the target values.
%inputs:
%   valsIn: vector of values to search in for targetVals
%   targetVals: vector or values to search for in valsIn
%outputs:
%   valsOut: vector of unique values that were closest to targetVals
%   ind: indicies of valsOut, i.e. valsOut = valsIn(ind);

    %initialize new delay vectors
    nTargetVals = length(targetVals);
    ind = zeros(nTargetVals,1);
    
    %find nearest delay to user specified target delay and record the index
    for ii = 1:nTargetVals
        [~,ind(ii)] = min(abs(valsIn - targetVals(ii)));
    end

    %update delay values based on index and remove any duplicate values
    [valsOut, ia] = unique(valsIn(ind),'stable');
    ind = ind(ia);
end