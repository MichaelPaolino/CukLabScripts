function obj = convertCList(obj,varargin)
% CONVERTCLIST converts data in the legacy/prototype conditionList database
% into a wlTR object. This method can load a .mat file directly or by
% lookup in the conditionList. To lookup by index, set the filepath to the
% folder containing dataList.mat and index to a non-zero value.
%
% obj = obj.convertCList(varargin)
%   Convert a conditionList entry to a wlTR object. Options for name-value
%   pair varargin are:
%       'file': char array that specifies the file path
%       'index': numeric scalar that specifies the conditionList index to
%           load.

    p = inputParser;
    p.FunctionName = 'convertCList';
    p.StructExpand = false;

    p.addParameter('file','');
    p.addParameter('index',0);

    % parse inputs and store results in struct p.Results
    p.parse(varargin{:}); 

    % User can either specify which file they load or the dataList path with an index
    if p.Results.index == 0   %If no index, load file directly
        filePath = p.Results.file;
        
    else %if given an index, load dataList file
        % Load conditionList from dataList.mat. If the path does not contain conditionList, this will throw an error
        myCList = load(fullfile(p.Results.file,'dataList.mat'),'conditionList');
        
        % List all files inside the directory
        cListDir = dir(p.Results.file);
        
        % Search for the conditionList entry. Append phonon removed if the conditionList is for phonon removed data
        if any(contains({cListDir(3:end).name},'phonon removed.mat'))   %if searching phonon removed dataList
            appendStr = ' phonon removed.mat';
        else %if searching avgCombinedOnly data list
            appendStr = '.mat';
        end
        
        % Create file name (without path)
        targetName = [myCList.conditionList{p.Results.index} appendStr];
        
        % Some download sources replace % with _, try searching for
        % replaced if not found
        if ~any(strcmp({cListDir(3:end).name},targetName))
            targetName = strrep(targetName,'%','_');    %replace % with _
        end    
        
        % load selected data to workspace
        filePath = fullfile(p.Results.file, targetName);
    end
    
    % Load data
    myData = load(filePath);
    
    % Once data is loaded into myData, copy contents into object
    obj.spectra.data = myData.dataCell{1,1}.'; %mOD [pixels, delays]
    obj.spectra_std.data = NaN(size(obj.spectra.data)); %does not contain std info, mOD [pixels, delays]
    obj.delays.data = myData.delayCell{1,1}(:); %ps [delays,1]
    obj.wavelengths.data = myData.wlCell{1,1}(:); %nm [pixels,1]
    obj.schemes = {'transient reflectance'};
    obj.t0.data = 0;
    
    % Update size struct
    obj.sizes = struct('nRpts', 1, ...  %struct with fields type double defining the size of spectra dimensions
                       'nGPos', 1, ...
                       'nSchemes', 1, ...
                       'nDelays', size(obj.spectra.data,2), ...
                       'nPixels', size(obj.spectra.data,1));
    
    % Update object name and shortName
    [~,fileName] = fileparts(filePath);
    obj.name = filePath;
    obj.shortName = fileName;
end


