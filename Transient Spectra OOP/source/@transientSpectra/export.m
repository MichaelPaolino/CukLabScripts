%export the object data to file
function [outputStruct, filePath] = export(obj,filePath,varargin)
    
    %set default export options for data. By default either a contour is
    %saved along with its wavelengths and delays or a spectra (wave) is
    %saved along with its wavelengths. The user can override these options.
    saveFlag = struct('delays',   obj.sizes.nDelays > 1,...  %true if multiple delays
                      'traces',   false,...                  %by default to not export traces
                      'spectra',  obj.sizes.nDelays == 1,... %false if there are multiple delays
                      'contours', obj.sizes.nDelays > 1);    %true if multiple delays
    
    %set default naming options for exported data when there are multiple schemes and grating positions
    nameFlag = struct('shortName', true,...
                      'scheme', false,...
                      'gPos', false,...
                      'rpts', false,...
                      'delay', false,...
                      'wl', false,...
                      'spectraUnit', true);
    
    %define subsets (todo: fix cannot be empty)
    wlSub = NaN;
    delaySub = NaN;  
    
    %Parse varargin
    %traces
    %spectra
    
    nSubWls = length(wlSub);
    nSubDelays = length(delaySub);
    
    %get actual indicies and wavelengths from object data. These may be
    %multi-dim if there are multiple repeats and grating positions
    [wlSubVal, wlSubInds] = nearestVal(obj.wavelengths,wlSub,'trim',false,'threshold',0.01*mean(obj.wavelengths,'all'));
    [delaySubVal, delaySubInds] = nearestVal(obj.delays,delaySub,'trim',false);
    
    %generate naming rules for data. These may be further modified for
    %contours, traces, and spectra below
    %todo: figure out how to select subsets of schemes, delays, etc.
    dataName = nameRule({'shortName', {obj.desc.shortName},  0, nameFlag.shortName;...
                         'scheme',     obj.schemes,     1, (length(schemes)>1) || nameFlag.scheme;...
                         'gPos',       strcat({'gPos: '},strcat(strtrim(cellstr(num2str(obj.gPos(:))))), ' nm'),         2, (length(gPos)>1) || nameFlag.gPos;...
                         'rpts',       strcat({'rpt '}, strtrim(cellstr(num2str((1:obj.sizes.nRpts)')))),                3, (length(rpts)>1) || nameFlag.rpts;...
                         'wl',         strcat(strtrim(cellstr(num2str(wlSub(:)))), unitStr.wl),             4, (length(wlSub)>1) || nameFlag.wl;...
                         'delay',      strcat(strtrim(cellstr(num2str(delaySub(:)))), unitStr.t),           4, (length(delaySub)>1) || nameFlag.delay;...
                         'dataUnit',   {['in' unitStr.sig]},                                                0, nameFlag.spectraUnit});
    
    outputStruct = struct('name','','data',[]);
    
    %decide what type of data to export depending on sizes
    %generate contours
    if saveFlag.contours
        %remove wavelengths and delay labels from naming rule and copy into contourName
        contourName = dataName.modify('wl.flag',false,'delay.flag',false);
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
        outputStruct = [outputStruct; contours]; 
    end
    
    %generate spectra
    if saveFlag.spectra
        %remove wavelengths labels from naming rule and copy into contourName
        spectraName = dataName.modify('wl.flag',false);
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
                    spectraName = spectraName.modify('delay.value',strcat(strtrim(cellstr(num2str(delaySubVal(:,rptInd,gPosInd))), unitStr.wl)));
                    for delayInd = delaySubInds
                        if ~isnan(delayInd)
                            %assign name and data with name autoindexing
                            [waves(cc).name, spectraName] = spectraName.buildName('autoIncrement',true,'delimiter',', '); 
                            waves(cc).data = obj.spectra.data(:,delayInd,rptInd,gPosInd,schemeInd);

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
        outputStruct = [outputStruct; waves]; 
    end
    
    %generate traces
    if saveFlag.traces
        %remove wavelengths labels from naming rule and copy into contourName
        traceName = dataName.modify('delay.flag',false);
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
                traceName = traceName.modify('wl.value',strcat(strtrim(cellstr(num2str(wlSubVal(:,gPosInd)))), unitStr.wl));
                for rptInd = 1:obj.sizes.nRpts
                    for wlInd = wlSubInds
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
        outputStruct = [outputStruct; waves]; 
    end
    
    %save delays
    
    
    %save wavelengths
    
    
end