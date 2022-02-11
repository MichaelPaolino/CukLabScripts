function outputStruct = export(obj,filePath,varargin)
% EXPORT object data to igor pro using an igor importable .mat file
% This method converts object spectra, wavelengths, and delays into matrix
% and wave data and saves it to a .mat file. The name of each matrix or
% column vector in the .mat file determines the name of the matrix or wave
% in the igor project. 
%
% The base name for each matrix or wave is determined by obj.desc.shortName
% and the type of matrix/wave that is being saved, e.g. matrix, spectra, 
% trace, wavelength, or delays. When multiple matricies or waves have the
% same base name, the names are made unique by automatically suffixing the 
% scheme, grating position, repeat, wavelength, and/or delay.
%
% By default, export saves matricies along with their delay and wavelength
% waves for data that contains multiple delays or a spectra wave along with
% its wavelength wave if there is only one delay. If the object contains
% multiple grating positions, each grating position gets a unique set of
% matricies or waves. If the object contains multiple repeats, the repeats
% are averaged together. These defaults can be overriden.
%
% The user also has the option to export kinetic traces and spectra using
% the 'kinetics', [wavelengths] and/or 'spectra', [delays] name-value
% pairs. 
%
% If exporting multiple objects in a loop, use the 'append', true name-
% value pair to append data to the same file.
%
% obj.EXPORT(filePath)
%   Exports the spectra, delays, and wavelengths inside obj to path
%   filePath in an igor-pro understandable format with default options,
%   as described.
%
% obj.EXPORT(filePath,'name1','value1','name2','value2',...)
%   Exports object data to filePath with options set by name-value pairs:
%   'kinetics', [wlVect]: additionaly exports kinetic traces as a function
%       of delay for wavelengths specified by wlVect (vector type double).
%       By default, kinetic traces are not exported.
%   'spectra', [delayVect]: additionaly exports spectra as a function of
%       wavelength for delays specified by delayVect (vector type double).
%       By default, spectra are not exported unless there is only one
%       delay.
%   'contours', flag: Logical flag on whether to export matrix data. The 
%       default is true unless there is only one delay.
%   'average', flag: Logical flag on whether to average together repeats.
%       The default is true.
%   'append', flag: Logical flag on whether to append saved matricies or 
%       waves to an existing .mat file. If the .mat file does not yet
%       exist, a new file is created. If the .mat file exists and has waves
%       or matricies with the same name as those being saved, the matricies
%       or waves will be replaced. The default value is false.
%
% outputStruct = obj.EXPORT(__)
%   Returns a structure containing the full (non-valid) name and data for 
%   each matrix and wave that is saved

    %set default export options for data. By default either a contour is
    %saved along with its wavelengths and delays or a spectra (wave) is
    %saved along with its wavelengths. The user can override these options.
    saveFlag = struct('delays',   obj.sizes.nDelays > 1,...  %true if multiple delays
                      'wavelengths', true,...                %by default always save wavelengths
                      'traces',   false,...                  %by default to not export traces
                      'spectra',  obj.sizes.nDelays == 1,... %by default false if there are multiple delays
                      'contours', obj.sizes.nDelays > 1);    %true if multiple delays
    
    %by default, append data to exisitng file
    appendFile = false;
    averageFlag = true;
    
    %set default naming options for exported data when there are multiple schemes and grating positions
    nameFlag = struct('shortName', true,... %show the short name
                      'scheme', false,...   %show the scheme name
                      'gPos', false,...     %show the grating position
                      'rpts', false,...     %show the rpt number
                      'delay', false,...    %show the delay number
                      'wl', false,...       %show the wavelength number
                      'spectraUnit', true); %show the data unit
    
    %define subsets
    wlSub = []; %wavelength subset to plot traces
    delaySub = 0;      %delay subset to plot spectra
    
    %parse varargin
    if nargin > 2
        for ii = 1:2:(nargin-2)
            assert(ischar(varargin{ii}),...
                ['Invalid argument class for name-value pair. Expected class char for name, got ' class(varargin{ii}) '.']);
            switch varargin{ii}
                case 'kinetics'  %User wants a subset of wavelengths to plot traces
                    wlSub = varargin{ii+1}; %update wavelength subset
                    saveFlag.traces = true; %set save traces flag to true
                    nameFlag.wl = true;
                case 'spectra' %User wants a subset of delays to plot spectra
                    delaySub = varargin{ii+1}; %update delay subset
                    saveFlag.spectra = true;   %set save spectra flag to true
                    nameFlag.delay = true;
                case 'contours' %user wants to change plot contours setting
                    saveFlag.contours = varargin{ii+1}; %update contour save flag
                case 'average'  %user wants to change default average behavior
                    averageFlag = varargin{ii+1};
                case 'append'   %user wants to change append to file setting
                    appendFile = varargin{ii+1};    
                otherwise
                    error([varargin{ii} ' is not a valid argument name.']); 
            end
        end
    end
    
    %average over repeats before exporting
    if averageFlag
        obj = obj.average();
    end
    
    %number of wavelengths and delays to save for traces and spectra
    nSubWls = length(wlSub);
    nSubDelays = length(delaySub);
    
    %get actual indicies and wavelengths from object data. These may be
    %multi-dim if there are multiple repeats and grating positions
    [wlSubVal, wlSubInds] = nearestVal(obj.wavelengths,wlSub,'trim',false,'threshold',0.01*mean(obj.wavelengths,'all'));
    [delaySubVal, delaySubInds] = nearestVal(obj.delays,delaySub,'trim',false);    
    
    %define unit strings
    unitStr = struct('spectra',obj.spectra.unit,...
                     'wl',obj.wavelengths.unit,...
                     'delay',obj.delays.unit);
    
    %generate naming rules for data. These may be further modified for
    %contours, traces, and spectra below
    %todo: figure out how to select subsets of schemes, delays, etc.
    dataName = nameRule({'shortName',  {obj.desc.shortName},                                                     0, nameFlag.shortName && ~isempty(obj.desc.shortName);...
                         'label',      {''},                                                                     0, true;...
                         'dataUnit',   {['in ' unitStr.spectra]},                                                0, nameFlag.spectraUnit;...
                         'scheme',     obj.schemes,                                                              1, (obj.sizes.nSchemes>1) || nameFlag.scheme;...
                         'gPos',       strcat({'gPos: '},strcat(strtrim(cellstr(num2str(obj.gPos(:))))), ' nm'), 2, (obj.sizes.nGPos>1) || nameFlag.gPos;...
                         'rpts',       strcat({'rpt '}, strtrim(cellstr(num2str((1:obj.sizes.nRpts)')))),        3, (obj.sizes.nRpts>1) || nameFlag.rpts;...
                         'wl',         strcat(strtrim(cellstr(num2str(wlSub(:)))), {' '}, unitStr.wl),           4, (length(wlSub)>1) || nameFlag.wl;...
                         'delay',      strcat(strtrim(cellstr(num2str(delaySub(:)))), {' '}, unitStr.delay),     4, (length(delaySub)>1) || nameFlag.delay});
    
    outputStruct = struct('name',{},'data',{});
    
    %decide what type of data to export depending on sizes
    %generate contours
    if saveFlag.contours
        %remove wavelengths and delay labels from naming rule and copy into contourName
        contourName = dataName.modify('label.values',{'Matrix'},...
                                      'wl.flag',false,...
                                      'delay.flag',false);
        
        %setup autoindexing levels for multiple name generation
        contourName = contourName.buildLevels(); 
        
        %Loop over data groups to do name conversion
        nContours = obj.sizes.nRpts*obj.sizes.nGPos*obj.sizes.nSchemes;
        contours = repmat(struct('name','','data',[]),nContours,1);
        
        %Populate countours structure with matrix and matrix name
        cc = 1;    
        for ii = 1:obj.sizes.nSchemes
            for jj = 1:obj.sizes.nGPos
                for kk = 1:obj.sizes.nRpts                   
                    %assign name and data with name autoindexing
                    [contours(cc).name, contourName] = contourName.buildName('autoIncrement',true,'delimiter',', '); 
                    contours(cc).data = obj.spectra.data(:,:,kk,jj,ii);

                    %increment contour counter
                    cc = cc + 1;
                end
            end
        end
        
        %add contour data to output structure
        outputStruct = [outputStruct; contours(~strcmp({contours.name},''))]; 
    end
    
    %generate spectra
    if saveFlag.spectra
        %remove wavelengths labels from naming rule and copy into contourName
        spectraName = dataName.modify('label.values',{'Spectra'},...
                                      'wl.flag',false);
        %setup autoindexing levels for multiple name generation
        spectraName = spectraName.buildLevels(); 
        
        %Loop over data groups to do name conversion
        nWaves = obj.sizes.nRpts*obj.sizes.nGPos*obj.sizes.nSchemes*nSubDelays;
        waves = repmat(struct('name','','data',[]),nWaves,1);

        %Populate countours structure with matrix and matrix name
        cc = 1;    
        for schemeInd = 1:obj.sizes.nSchemes
            for gPosInd = 1:obj.sizes.nGPos
                for rptInd = 1:obj.sizes.nRpts
                    %update the delay name values to the actual delay values for this repeat and grating position
                    spectraName = spectraName.modify('delay.values',strcat(strtrim(cellstr(num2str(delaySubVal(:,rptInd,gPosInd),3))), ' ', unitStr.wl));
                    for delayInd = delaySubInds(:,rptInd,gPosInd)'
                        if ~isnan(delayInd)
                            %assign name and data with name autoindexing
                            [waves(cc).name, spectraName] = spectraName.buildName('autoIncrement',true,'delimiter',', '); 
                            waves(cc).data = obj.spectra.data(:,delayInd,rptInd,gPosInd,schemeInd);
                            waves(cc).data = waves(cc).data(~isnan(waves(cc).data));
                            
                            %increment contour counter
                            cc = cc + 1;
                        else
                            spectraName = spectraName.increment();
                        end
                    end
                end
            end
        end
        
        %add spectra waves to output structure
        outputStruct = [outputStruct; waves(~strcmp({waves.name},''))]; 
    end
    
    %generate traces
    if saveFlag.traces
        %remove wavelengths labels from naming rule and copy into contourName
        traceName = dataName.modify('label.values',{'Traces'},...
                                    'delay.flag',false);
        %setup autoindexing levels for multiple name generation
        traceName = traceName.buildLevels(); 
        
        %Loop over data groups to do name conversion
        nWaves = obj.sizes.nRpts*obj.sizes.nGPos*obj.sizes.nSchemes*nSubWls;
        waves = repmat(struct('name','','data',[]),nWaves,1);

        %Populate countrous structure with matrix and matrix name
        cc = 1;    
        for schemeInd = 1:obj.sizes.nSchemes
            for gPosInd = 1:obj.sizes.nGPos
                %Update the wavelength name values for the actual available wavelengths for the grating position
                traceName = traceName.modify('wl.values',strcat(strtrim(cellstr(num2str(wlSubVal(:,gPosInd),3))), unitStr.wl));
                for rptInd = 1:obj.sizes.nRpts
                    for wlInd = wlSubInds(:,gPosInd)'
                        if ~isnan(wlInd)    %If the wavelength is outside the grating position it is flagged with NaN
                            %assign name and data with name autoindexing
                            [waves(cc).name, traceName] = traceName.buildName('autoIncrement',true,'delimiter',', '); 
                            waves(cc).data = obj.spectra.data(wlInd,:,rptInd,gPosInd,schemeInd).';

                            %increment contour counter
                            cc = cc + 1;
                        else
                            traceName = traceName.increment();
                        end
                    end
                end
            end
        end
        
        %add trace waves to output structure
        outputStruct = [outputStruct; waves(~strcmp({waves.name},''))]; 
    end
    
    %save delays
    if saveFlag.delays
        %remove wavelengths labels from naming rule and copy into contourName
        delayName = dataName.modify('label.values',strcat({'Delays in '},unitStr.delay),...
                                    'scheme.flag',false,...
                                    'wl.flag',false,...
                                    'delay.flag',false,...
                                    'dataUnit.flag',false);
        %setup autoindexing levels for multiple name generation
        delayName = delayName.buildLevels(); 
        
        %Loop over data groups to do name conversion
        nWaves = obj.sizes.nRpts*obj.sizes.nGPos;
        waves = repmat(struct('name','','data',[]),nWaves,1);

        %Populate countrous structure with matrix and matrix name
        cc = 1;    
        for gPosInd = 1:obj.sizes.nGPos
            for rptInd = 1:obj.sizes.nRpts
                %assign name and data with name autoindexing
                [waves(cc).name, delayName] = delayName.buildName('autoIncrement',true,'delimiter',', '); 
                waves(cc).data = obj.delays.data(:,rptInd,gPosInd);

                %increment contour counter
                cc = cc + 1;
            end
            
        end
        
        %add trace waves to output structure
        outputStruct = [outputStruct; waves(~strcmp({waves.name},''))]; 
    end    
    
    %save wavelengths
    if saveFlag.wavelengths
        %remove wavelengths labels from naming rule and copy into contourName
        wlName = dataName.modify('label.values',strcat({'Wavelengths in '},unitStr.wl),...
                                    'scheme.flag',false,...
                                    'rpts.flag',false,...
                                    'wl.flag',false,...
                                    'delay.flag',false,...
                                    'dataUnit.flag',false);
        %setup autoindexing levels for multiple name generation
        wlName = wlName.buildLevels(); 
        
        %Loop over data groups to do name conversion
        nWaves = obj.sizes.nGPos;
        waves = repmat(struct('name','','data',[]),nWaves,1);

        %Populate countrous structure with matrix and matrix name
        cc = 1;    
        for gPosInd = 1:obj.sizes.nGPos
                %assign name and data with name autoindexing
                [waves(cc).name, wlName] = wlName.buildName('autoIncrement',true,'delimiter',', '); 
                waves(cc).data = obj.wavelengths.data(:,gPosInd);

                %increment contour counter
                cc = cc + 1;            
        end
        
        %add trace waves to output structure
        outputStruct = [outputStruct; waves(~strcmp({waves.name},''))]; 
    end
    
    %Convert outputStruct name to structure with field names and save as a
    %file
    fullNames = {outputStruct.name}';
    validNames = matlab.lang.makeValidName(fullNames);
    validNames = matlab.lang.makeUniqueStrings(validNames);
    
    %build save structure with field names and save to file
    saveCell = [validNames, {outputStruct.data}']';
    saveStruct = struct(saveCell{:});
    
    if appendFile && isfile(filePath)   %append file if option is set to true and the file exists
        %load previous names
        saveStrings = load(filePath,'validNames','fullNames');
        
        %update name valid and full name list 
        validNames = [validNames; strsplit(saveStrings.validNames,';')'];
        fullNames = [fullNames; strsplit(saveStrings.fullNames,';')'];
        
        %remove duplicates while keeping name order
        [validNames,vi] = unique(validNames,'stable');
        fullNames = fullNames(vi);
        
        %update saveStruct and save to file
        saveStruct.validNames = strjoin(validNames,';');
        saveStruct.fullNames = strjoin(fullNames,';');
        save(filePath,'-struct','saveStruct','-append');        
    else %save without appending file
        %update saveStruct and save to file
        saveStruct.validNames = strjoin(validNames,';');
        saveStruct.fullNames = strjoin(fullNames,';');
        save(filePath,'-struct','saveStruct');
    end

    
end