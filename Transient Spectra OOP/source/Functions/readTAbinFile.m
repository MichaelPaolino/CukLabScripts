function [data, nm, T, exp, wlc, rawdata] = readTAbinFile(filePath, varargin)
% READTABINFILE extracts experimental data packaged in a TA bin file. The 
% current version will attempt to extract any usable data from from a 
% incomplete bin file. Use the additional options to control partial repeat 
% loading.
%
% The exp structure contains relavent experimental parameters, such as the
% description and acquisition metrics in the exp.labelValues field. The
% labelValues field is an ND array of size:
%   [overall repeats, repeats per delay, delay, value index]
% The elements of value index are:
%   [acq start (ms timestamp), acq stop (ms timestamp),...
%       actual delay readout (ps), x (mm), y (mm), x (mm)] 
%
% input: filePath = string for file path
%        [nRepeats] = if program crashes, last good known repeat to extract
%                     data
%        [avgRpts] = repeats range (i.e. 2:8 for 8 repepats) to include in
%        the average for data if bin file contains bad repeats.
%
% output: data = averaged data [nGpos x nDelays x nPixels]
%         nm = spectral calibration [nPixels x nGpos]
%         T = delays [nDelays x 1]
%         exp = structure that contains experimental parameters
%         wlc = white light spectrum [nRepeats, nGpos, nDelays, nPixels]
%         rawdata = unaveraged data [nRepeats, delay repeats, nGpos, nDelay, nPixels]

    % Open for reading and set to big-endian binary format
    fid = fopen(filePath,'r','ieee-be');  

    %bin file is ordered as header, data
    %1) Exp settings
    %a) Description (string)
    %b) Averages array (long)
    %c) Grating positions array (double)
    %d) Stage positions array (double)
    nStr = fread(fid,1,'uint32');
    exp.description = convertCharsToStrings(char(fread(fid,nStr,'unsigned char')));
    nAvgs = fread(fid,1,'uint32');
    exp.averages = fread(fid,nAvgs,'int32'); %frames, grating pos repeats, overall, per delay
    nGratings = fread(fid,1,'uint32');
    exp.gratings = fread(fid,nGratings,'double');
    nDelays = fread(fid,1,'uint32');
    exp.mmdelays = fread(fid,nDelays,'double');

    %2) Array calibration - saves as a 1D array
    nAxis = fread(fid,1,'uint32');
    nPixels = fread(fid,1,'uint32');
    nm = reshape(fread(fid,nAxis*nPixels,'double')',[nPixels,nAxis]);

    %overall repeat
    %stage position
    %stage position repeat
    
    %If the program keeps crashing and you want to keep a portion of the 
    % data, update the averages array as follows:
    avgRange = 0;
    if ~isempty(varargin)
        if ~isempty(varargin{1})
            exp.averages(3) = varargin{1};
        end
        if length(varargin)>1
            if ~isempty(varargin{2})
                avgRange = varargin{2};
            end
        end
    end
    
    if avgRange == 0
        avgRange = 1:exp.averages(3);
    end
    
    %labelValues dims: [overall repeats, repeats per delay, delay, value index]
    %value index corresponds to values in this vector: 
    %[acq start (ms timestamp), acq stop (ms timestamp), actual delay readout (ps), x (mm), y (mm), x (mm)] 
    nLabelValues = fread(fid,1,'uint32');
    labelValues = NaN(exp.averages(3), exp.averages(4), nDelays, nLabelValues);
    
    %spectra dim: [overall repeats, repeats per delay, grating positions, nm, delay]
    rawdata = NaN(exp.averages(3), exp.averages(4), nGratings, nDelays, nPixels);
    
    %WLCspectra dim: [overall repeats, repeats per delay, grating positions, nm, delay]
    wlc =  NaN(exp.averages(3), nGratings, nDelays, nPixels);
    
    % try loading data. fread returns an error if the file is incomplete
    try
        ind = 0;
        for ii = 1:exp.averages(3) %overall repeats
           for jj = 1:nDelays %number of time points
               for kk = 1:exp.averages(4) %repeats per delay
                    ind = ind + 1;
                    %3) Misc data
                    %a) data value array (double)
                    %b) data label array (unit16)
                    if ~all([ii,jj,kk]==1)
                        nLabelValues = fread(fid,1,'uint32');
                    end
                    labelValues(ii,kk,jj,:) = fread(fid,nLabelValues,'double');
                    nLabels = fread(fid,1,'uint32');
                    fread(fid,nLabels,'uint16');	%steps the cursor forward by uint16

                    %4) Spectra (2D array)
                    nDim1 = fread(fid,1,'uint32');
                    nDim2 = fread(fid,1,'uint32');
                    rawdata(ii,kk,:,jj,:) = permute(reshape(fread(fid,nDim1*nDim2,'double')',[nDim2,nDim1]),[2,3,1]);%[pixel, grating pos]->[grating pos, delay, pixel]
               end

               %4) Spectra (2D array)
               nDim1 = fread(fid,1,'uint32');
               nDim2 = fread(fid,1,'uint32');
               wlc(ii,:,jj,:) = permute(reshape(fread(fid,nDim1*nDim2,'double')',[nDim2,nDim1]),[2,3,1]); %[pixel, grating pos]
           end
        end
    catch
        %Warn user that the file was incomplete
        warning('TA bin file was incomplete.\nImport stopped at index [repeat, delay, rpt/delay] = [%d, %d, %d].\nFile: %s',ii,jj,kk,filePath);
    end
    
    % close disk file reference
    fclose(fid);
    
    %raw data: [nRepeats, delay repeats, nGpos, nDelay, nPixels]
    %data: [nGpos x nDelays x nPixels]
    data = permute(mean(rawdata(avgRange,:,:,:,:),[1,2],'omitnan'),[3,4,5,1,2]); %all repeats
    
    %labelValues: [overall repeats, repeats per delay, delay, value index]
    %T: [delay x 1]
    T = permute(mean(labelValues(avgRange,:,:,3),[1,2],'omitnan'),[3,1,2]); %Delays
    exp.labelValues = labelValues;
end
