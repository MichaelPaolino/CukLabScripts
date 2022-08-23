function obj = convertCList(obj,varargin)
% CONVERTCLIST converts data in the legacy/prototype conditionList database
% into a wlTR object.
    p = inputParser;
    p.FunctionName = 'convertCList';
    p.StructExpand = false;

    p.addParameter('sourceFile','');
    p.addParameter('sourceIndex',0);

    % parse inputs and store results in struct p.Results
    p.parse(varargin{:}); 

    % User can either specify which file they load or the dataList path with an index
    %if p.Results.sourceIndex == 0   %If no index, load file directly
        myData = load(p.Results.sourceFile);
    %else %if given an index, load dataList file
    %    myCList = load(p.Results.sourceFile,'conditionList');
        %tmpPath = myCList.conditionList{p.Results.sourceIndex};
    %end
    
    % Temporary code to allow conditionList data to be loaded
    obj.spectra.data = myData.dataCell{1,1}.'; %mOD [pixels, delays]
    obj.spectra_std.data = zeros(size(obj.spectra.data)); %does not contain std info, mOD [pixels, delays]
    obj.delays.data = myData.delayCell{1,1}(:); %ps [delays,1]
    obj.wavelengths.data = myData.wlCell{1,1}(:); %nm [pixels,1]
    obj.schemes = {'TR'};

    obj.sizes = struct('nRpts', 1, ...  %struct with fields type double defining the size of spectra dimensions
                       'nGPos', 1, ...
                       'nSchemes', 1, ...
                       'nDelays', size(obj.spectra.data,2), ...
                       'nPixels', size(obj.spectra.data,1));

    obj.name = p.Results.sourceFile;
end


