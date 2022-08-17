classdef wlTR < transientSpectra
   properties
       chirpParams = 0;
       
       chirpCorrected = false;
   end
   
   properties (Access = protected)
       TRflag = struct('chirpFit',false);     
   end
   
   methods
       
       function obj = setChirp(obj, chirpFit)
        % SETCHIRP sets the chirp parameters for every element in the object array.
        % This method accepts a *.mat file with a variable called chirpFit or a
        % polyval vector of chirp parameters.
        %
        % obj = obj.SETCHIRP(chirpVector)
        % obj = obj.SETCHIRP(filePath)
        %   Sets the chirp parameters for every element of obj
        %
        % See Also: fitChirp, correctChirp, polyval
           
           % Use has an option to load chirp parameters from a .mat file
           if ischar(chirpFit)
              tmp = load(chirpFit);
              chirpFit = tmp.chirpFit;
           elseif ~isvector(chirpFit)
              error('chirpParams expected path or polyval vector.');
           end
           
           % Formtat object array dims into a column for easy looping
           objSize = size(obj);
           objNumel = numel(obj);
           obj = obj(:);
           
           %loop through each object and update chirp params
           for objInd = 1:objNumel
               obj(objInd).chirpParams = chirpFit;
           end
               
           %convert object back to original array dims
           obj = reshape(obj,objSize);
       end
       
       function obj = correctChirp(obj, varargin)
        % CORRECTCHIRP corrects the probe wavelength-dependent group delay (chirp)
        % by interpolating spectra data on a chirp-corrected delay. This
        % chirp-corrected delay is wavelength dependent and is determined by a
        % chirp correction polynomial. This polynomial is calculated by calling the
        % findChirp method. This polynomial is usually found from a seperate
        % calibration spectrum (e.g.high fluence OC with delays collected from
        % -3:0.1:3 ps). 
        %
        % Chirp correcion is applied as a wavelength-dependent offset for the
        % delay. This offset is calculated with respect to a reference wavelength,
        % which has zero offset. By default, this wavelength is the central
        % wavelength in each object element.
        %
        % This method also provides some degree of cosmetic processing to handle 
        % data clipping and extrapolation which occurs at the first and last 
        % delay. Clipping and extrapolation occurs for each wavelength that has a
        % non-zero chirp correction. By default, extrapolated values are also set 
        % to NaN.
        %
        % obj = obj.CORRECTCHIRP()
        %   Interpolates spectral data inside obj using internally assigned chirp 
        %   paramters. This call uses the default options of using the center 
        %   wavelength as the reference wavelength and without any extrapolation.
        %
        % obj = obj.CORRECTCHIRP(varargin)
        %   Interpolates spectra data inside obj with additional name-value pair
        %   options.
        %
        % Name-Value Pairs
        %   'chirpParams': (vector double) a vector of polynomial coefficients 
        %       (as defined by polyfit) to use for correction and to assign to all 
        %       object elements. The default is to use the object element member
        %       data.
        %   'extrap': (char or scalar) a constant value to assign to extrapolated
        %       points or a flag to override default interp1 extrap behavior. 
        %       Allowed flags are: 'none','extrap'. 'none' replaces extrapolated 
        %       points with NaN, and 'extrap' uses interp1 extrapolation. Default
        %       is 'none'.
        %   'wlRef': (char or scalar) the wavelength to use as the reference
        %       wavelength or a flag to calculate the reference wavelength. Allowed
        %       flags are: 'min', 'mid', or 'max' which choose the minimum, middle,
        %       or maximum wavelength in each object element. Default is 'mid'.
        %   'interp': (char) the interpolation method to use as defined by interp1.
        %       Default is 'linear'.
        %
        % See Also: FITCHIRP, INTERP1, POLYVAL
           
           % Format object array dims into a column for easy looping
           objSize = size(obj);
           objNumel = numel(obj);
           obj = obj(:);
           
           % Format varargin using input parser
           p = inputParser;
           p.FunctionName = 'correctChirp';
           p.addParameter('chirpParams', [], @(p) isvector(p));
           p.addParameter('extrap', 'none', @(p) (isscalar(p) && ~ischar(p)) || any(strcmp(p,{'none','extrap'})));
           p.addParameter('wlRef','mid');
           p.addParameter('interp','linear');
           
           % Parse arguemnts, results will be in p.Results
           p.parse(varargin{:});
           
           % Loop over object array
           for objInd = 1:objNumel
               % Temporarily set units to match chirp correction polynomial
               tmpUnits = cell(3,1);
               [tmpUnits{:}] = obj(objInd).getUnits;
               obj(objInd) = obj(objInd).setUnits('nm','ps','');
               
               % Extract data for numerical processing
               spectra = obj(objInd).spectra.data;
               t = obj(objInd).delays.data;
               wl = obj(objInd).wavelengths.data;
               
               % Assign external chirp parameters
               if ~isempty(p.Results.chirpParams)
                   obj(objInd).chirpParams = p.Results.chirpParams;
               end
               
               % Determine what wavelength is the delay reference point wlRef
               if ischar(p.Results.wlRef)
                   switch p.Results.wlRef
                       case 'min'
                           wlRef = min(wl,[],'all','omitnan');
                       case 'max'
                           wlRef = max(wl,[],'all','omitnan');
                       case 'mid'
                           wlRef = mean([min(wl,[],'all','omitnan'), max(wl,[],'all','omitnan')]);
                       otherwise
                           error(['Expected min, max, or mid for wlRef, got ' wlRef '.']);
                   end
               else
                   wlRef = p.Results.wlRef;
               end
               
               % Evalaute t0 shifts by evaluation the polynomial at the reference point wlRef
               tRef = polyval(obj(objInd).chirpParams,wlRef);
               
               %**This comment is used a reference only**
               %reshape arrays to [delays, unique, everything else]
               % spectra: [pixels, delays, rpts, grating pos, schemes]
               % delays: [delays, rpts, gPos]
               % wavelengths: [pixels, gPos]
               % tShift: [pixels, gPos]
               
               % Calculate the delay shift amount for the specific wavelengths
               tShift = polyval(obj(objInd).chirpParams, wl)-tRef;
               
               % Prepare arguments for interpolation
               args = cell(4,1);
               args{4} = p.Results.interp;
               
               % Directly pass 'extrap' to interp1
               if strcmp(p.Results.extrap,'extrap')
                   args{5} = 'extrap';
               end
               
               % Loop over all unique delay axes and interpolate spectra 
               sNearest = zeros(size(spectra));
               for gPInd = 1:obj(objInd).sizes.nGPos
                   for rptInd = 1:obj(objInd).sizes.nRpts
                       for pixelInd = 1:obj(objInd).sizes.nPixels
                           % Use linear interpolation to correct for chirp
                           args{1} = t(:,rptInd,gPInd); %old delay axis
                           args{2} = spectra(pixelInd,:,rptInd,gPInd,:); %old spectra values
                           args{3} = t(:,rptInd,gPInd)+tShift(pixelInd,gPInd); %shifted delay axis
                           
                           %interpolate spectra values on shifted delay axis and assign to old delay axis
                           spectra(pixelInd,:,rptInd,gPInd,:) = interp1(args{:});
                       end
                   end
               end
               
               % Handle extrapolation if the default extrap behavior of interp1 is not what the user wants
               if (isscalar(p.Results.extrap) && ~ischar(p.Results.extrap)) || any(strcmp(p.Results.extrap,{'none','nearest'}))
                   % find extrapolated delays (vectorized version)
                   % reshape tShift and t to match dims of spectra, and tShift for size of spectra
                   tShift2 = permute(tShift,[1,3,4,2]); %[pixels, 1, 1, gPos]
                   tShift2 = repmat(tShift2,size(spectra,1:4)./size(tShift2,1:4)); %[pixels, delays, rpts, gPos, schemes]
                   t2 = permute(t,[4,1,2,3]); %[1, delays, rpts, gPos]
                   
                   % find the minimum extrap delays
                   isExtrapMin = t2(:,1,:,:) > (t2 + tShift2);
                   isExtrapMax = t2(:,end,:,:) < (t2 + tShift2);
                   inds = repmat(or(isExtrapMin,isExtrapMax), [1,1,1,1,size(spectra,5)]);
                   
                   % handle different extrap options
                   if isscalar(p.Results.extrap) && ~ischar(p.Results.extrap)
                       spectra(inds) = p.Results.extrap;
                   else
                       switch p.Results.extrap
                           case 'none'
                               % set all extrap values to NaN
                               spectra(inds) = NaN;
    %                        case 'nearest' todo: finish writing this...
    %                            % set all extrap values to nearest
    %                            spectra(inds) = sNearset(inds);
                       end
                   end
               end

               % Assign interpolated spectra back to object data
               obj(objInd).spectra.data = spectra;
               
               %Change units back to original units
               obj(objInd) = obj(objInd).setUnits(tmpUnits{:});
           end
           
           %convert object array back to original size
           obj = reshape(obj,objSize);
       end
       
       function [obj,chirpFit] = fitChirp(obj,varargin) 
        % FITCHIRP attempts to extract the white light continuum group delay as a 
        % function of wavelength. The group delay is approximated with a polynomial
        % function (default: 5th order). The ideal data set is a high fluence OC 
        % spectra with finely spaced delay points (~100 fs) with a range of 
        % +/- 2 ps around t0.
        %
        % This method works by fitting a specified delay range (Default, -2 to 2 ps)
        % to an erf sigmoid for every wavelength in the dataset. The sigmoid's 
        % center vs. wavelength is then fit to a polynomial of specified order 
        % (default: 5th order).
        %
        % obj = obj.FITCHIRP()
        %   Extracts the WLC group delay vs. wavelength and fits the data to a 5th
        %   order polynomial.
        %
        % obj = obj.FITCHIRP(varargin)
        %   Extracts the WLC group delay vs. wavelengths and fits to a polynomial
        %   with the additional name-value pair options.
        %
        % [obj, chirpParam] = obj.FITCHIRP(__)
        %   Additionally returns array of size [polyOrder, objSize].
        %
        % Name-value pairs:
        %   'wavelengths': (double array) an array of wavelengths to subset the
        %      spectra to help speed up fitChirp. Default is to use all
        %      wavelengths.
        %   'delays': (double array) [lower, upper] a 1x2 array of delays to trim
        %      that limits the sigmoid fit range. Default is -2 to 2 ps.
        %   'order': (int) an integer > 0 that specifies the polynomial order of
        %      the group delay vs. wavelength. The default order is 5.
        %
        % See Also: CORRECTCHIRP, FINDT0, CORRECTT0, LSQFITSIGMOID,
        % POLYFIT, POLYVAL
           
           % Format object array dims into a column for easy looping
           objSize = size(obj);
           objNumel = numel(obj);
           obj = obj(:);
           
           % Format varargin using input parser
           p = inputParser;
           p.FunctionName = 'fitChirp';
           p.addParameter('wavelengths', []);
           p.addParameter('delays', [-2,2]);
           p.addParameter('order',5);
           
           % Parse arguemnts, results will be in p.Results
           p.parse(varargin{:});
           
           % Trim and subset object size to desired sub-range
           objTmp = obj.trim('delays',p.Results.delays);
           
           if ~isempty(p.Results.wavelengths)
               objTmp = objTmp.subset('wavelengths',p.Results.wavelengths);
           end
           
           % Set units to common values so that initial guess below can be
           % referenced in other data sets
           objTmp = objTmp.setUnits('nm','ps','mOD');
           
           % Average repeats and stitch grating position to make the
           % procedure below easier to run
           objTmp = objTmp.average();
           objTmp = objTmp.stitch();
           
           % Initialize chirpFit output
           chirpFit = zeros(p.Results.order+1, objNumel);
           
           % Start waitbar
           f = waitbar(0,['Fitting group delay for spectra 1 of ' num2str(objNumel)]);
           
           % Loop over individual object elements
           for objInd = 1:objNumel
               % Update waitbar
               waitbar(0,f,['Fitting group delay for spectra ' num2str(objInd) ' of ' num2str(objNumel)]);
               
               % Extract spectra data from object and use locally
               data = objTmp(objInd).spectra.data; %[wls,delays]
               wl = objTmp(objInd).wavelengths.data; %[wls,1]
               t = objTmp(objInd).delays.data;  %[delays,1]
               
               % Initialize sigmoid fit parameter matrix
               myFP = zeros(length(wl),4);
               
               % Loop over wavelengths
               for wlInd = 1:length(wl)
                   % Update waitbar
                   waitbar(wlInd/length(wl),f);
                   
                   % Do the least squares fit using the lsqsigmoid fit
                   myFP(wlInd,:) = lsqFitSigmoid(t, data(wlInd,:));
               end
               
               % Once the sigmoid fit is done, find the chirp parameters
               % and update object data
               chirpFit(:,objInd) = polyfit(wl,myFP(:,3),p.Results.order);
               obj(objInd).chirpParams = chirpFit(:,objInd);
               
               % Plot resuls
               % todo: make an option to display (and possibly log?) results
               figure;
                    contour(wl,t,data',25);
                    hold on;
                    p1 = plot(wl,myFP(:,3),'r','DisplayName','sigmoid t0');
                    p2 = plot(wl,polyval(chirpFit,wl),'k--','DisplayName','polynomial');
                    hold off;
                    
               xlabel('Wavelength (nm)');
               ylabel('Delay (ps)');
               legend([p1,p2]);
               colorbar();
           end
           
           close(f);
           
           %convert object array back to original size
           obj = reshape(obj,objSize);
           chirpFit = reshape(chirpFit, [p.Results.order+1,objSize]);
       end
       
%        function [obj, phononObj] = fitPhonons(obj, phononModel)
%            
%        end
%        
%        function obj = removePhonons(obj, varargin)
%            
%        end
   end
end