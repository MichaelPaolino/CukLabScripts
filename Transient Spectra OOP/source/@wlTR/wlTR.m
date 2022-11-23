classdef wlTR < transientSpectra
   properties
       chirpParams = 0; % (vector double) a polfit vector fitting white light group delay vs. wavelength
   end
   
   properties (Access = protected)
       TRflag = struct('chirpFit',false);     
   end
   
   % Constructor method
   methods
       function obj = wlTR(varargin)
        % A WLTR object contains white-light transient reflectance-specific 
        % funtionality in addition to everything inside the transientSpectra class.
        % The main additional features are related to chirp correction and loading 
        % legacy TA binary and conditionList data.
        % 
        % obj = wlTR(__);
        %   Constructs a wlTR object with the same arguments as the transientSpectra
        %   constructor call.
        %
        % obj = wlTR(__, 'loadType', lType);
        %   Constructs a wlTR object with data imported from a file using the
        %   import method specified with the 'loadType' name-value pair. lType is a
        %   char array or cell of char arrays. Allowed values are:
        %       'dataHolder': loads a data holder .mat file created by the LabVIEW 
        %           MATLAB API using the new acquisition program
        %       'bin': loads a legacy TA binary file acquired by the old
        %           acquisition program
        %       'cListFile': loads a .mat file that is part of the legacy
        %           conditionList database
        %       'cListIndex': loads a file specified by conditionList index, where
        %           conditionList is a cell array inside dataList.mat. When using
        %           this loadType, the file path should navigate to the directory
        %           containing the dataList.mat file. Instead of specifying a .mat
        %           file, specify the index inside the conditionList cell array. 
        %           For example: 
        %            obj = wlTR('..\Phonon Removed\5','loadType','cListIndex'); 
        %           will load the data for the 5th entry in the phonon removed
        %           conditionList. Note that this loadType only works with full
        %           paths, so the .. in the example will have to be replaced with
        %           the full path for the phonon removed folder.
        %       
        % WLTR has the same units as transientSpectra.
        %
        % The default units for a wlTR object are mOD, ps, and nm
        %
        % See also: TRANSIENTSPECTRA, DOUBLEWITHUNITS
           
           % Currently, wlTR does not define additional constructor
           % functionality. Call the transientspectra superclass
           % constructor directly.
           obj@transientSpectra(varargin{:});
       end
   end
   
   % Methods that define the inner workings of the class
   methods (Access = 'protected')
        function obj = importData(obj,myPath,loadType)
        % IMPORTDATA loads a file from a path using the import method specified by
        % a load type and populates the objects member data. This method is 
        % designed to dispatch the correct import/load method and not implement a
        % specific import/load routine, therby allowing the developer to override 
        % import implementations in subclasses. 
        %
        % The wlTR class adds implementation for loading raw TA/TR binary files
        % generated by the older version of the LabVIEW TR Acquisition program and
        % for loading processed data stored in conditionList. ConditionList is a
        % database prototype that was used to log TR data acquired from 2019 to
        % 2021.
        %
        % obj = obj.LOADPATH(path, loadType);
        %   Loads data into the object from path depending on the load type. Load
        %   types that are supported for wlTR are:
        %   'cListFile': a .mat file that corresponds to a conditionList entry
        %   'cListIndex': a condition list index specified in path by:
        %       path = '..\index'. For example, '..\Phonon Removed\1' will load the
        %       first entry in conditionList in the ..\Phonon Removed\ path 
        %
        % See Also: TRANSIENTSPECTRA
        
            switch loadType
                case 'cListFile'
                    % Call the convertCList import method to convert a conditionList file to a wlTR object
                    obj = convertCList(obj,'file',myPath);
                case 'cListIndex'
                    % Call the convertCList import method to convert a conditionList index to a wlTR object
                    [myPath,ind] = fileparts(myPath); %extract the parent 
                    obj = convertCList(obj,'file',myPath,'index',str2double(ind));
                case 'bin'
                    % Call convertTABin importmethod to load and convert a TA bin file into a wlTR object
                    obj = convertTABin(obj,myPath);
                otherwise
                    % Call superclass importData method for loading a dataholder or handle unknown loadType
                    obj = importData@transientSpectra(obj,myPath,loadType);
            end
        end
        
        % Converts a conditionList file to a wlTR object
        obj = convertCList(obj,varargin)
   end
   
   % Data correction and manipulation methods
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
           
           % User has an option to load chirp parameters from a .mat file
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
        %   'chirpParams': (vector double or char array) chirp parameters to be
        %       passed to the setChirp method. The value can either be a vector of 
        %       chirp parameters (as defined by polfit) or a path to a chirp 
        %       parameter .mat file.
        %   'wlRef': (char or scalar) the wavelength to use as the reference
        %       wavelength or a flag to calculate the reference wavelength. Allowed
        %       flags are: 'min', 'mid', or 'max' which choose the minimum, middle,
        %       or maximum wavelength in each object element. Default is 'mid'.
        %   'interp': (char) the interpolation method passed to griddedInterpolant.
        %       The default interpolation method is 'linear'.
        %   'extrap': (char) an extrapolation method passed to griddedInterpolant.
        %       The default extrapolation method is 'nearest'.
        %
        % See Also: FITCHIRP, SETCHIRP, GRIDDEDINTERPOLANT, POLYVAL
           
           % Format object array dims into a column for easy looping
           objSize = size(obj);
           objNumel = numel(obj);
           obj = obj(:);
           
           % Format varargin using input parser
           p = inputParser;
           p.FunctionName = 'correctChirp';
           p.addParameter('chirpParams', []);
           p.addParameter('wlRef','mid', @(p) (ischar(p) && any(strcmp(p,{'min','mid','max'})) || isscalar(p)));
           p.addParameter('interp','linear');
           p.addParameter('extrap','nearest');
           
           % Parse arguemnts, results will be in p.Results
           p.parse(varargin{:});
           
           % Setup a griddedInterpolant object
           F = griddedInterpolant();
           
           % Assign external chirp parameters
           if ~isempty(p.Results.chirpParams)
               obj = obj.setChirp(p.Results.chirpParams);
           end
           
           % Loop over object array
           for objInd = 1:objNumel
               % Temporarily set units to match chirp correction polynomial
               tmpUnits = cell(3,1);
               [tmpUnits{:}] = obj(objInd).getUnits;
               obj(objInd) = obj(objInd).setUnits('nm','ps','');
               
               % Extract data for numerical processing
               s = obj(objInd).spectra.data;
               t = obj(objInd).delays.data;
               l = obj(objInd).wavelengths.data;
               
               % Determine what wavelength is the delay reference point wlRef
               if ischar(p.Results.wlRef)
                   switch p.Results.wlRef
                       case 'min'
                           wlRef = min(l,[],'all','omitnan');
                       case 'max'
                           wlRef = max(l,[],'all','omitnan');
                       case 'mid'
                           wlRef = mean([min(l,[],'all','omitnan'), max(l,[],'all','omitnan')]);
                       otherwise
                           error(['Expected min, max, or mid for wlRef, got ' wlRef '.']);
                   end
               else
                   wlRef = p.Results.wlRef;
               end
               
               % Convert array sizes of t and l to match the size of s for easy vectorization
               % s old: [pixels, delays, rpts, g pos, schemes]
               % s new: [pixels, delays, rpts x g pos x schemes]
               % t old: [delays, rpts, g pos]
               % t new: [delays, rpts x g pos x schemes]
               % l old: [pixels, g pos]
               % l new: [pixels, rpts x g pos x schemes]
               t = reshape(explicitExpand(t,sizePadded(s,2:5)),obj(objInd).sizes.nDelays,[]);
               l = reshape(explicitExpand(permute(l,[1,3,2]),sizePadded(s,[1 3:5])),obj(objInd).sizes.nPixels,[]);
               s = reshape(s,obj(objInd).sizes.nPixels,obj(objInd).sizes.nDelays,[]);
               
               % Evalaute t0 shifts by evaluation the polynomial at the reference point wlRef
               tRef = polyval(obj(objInd).chirpParams,wlRef);
               
               % loop over extra dims
               for ii = 1:size(s,3)
                   % MATLAB's griddedInterpolant requires that grid vectors are sorted
                   [t(:,ii), tI] = sort(t(:,ii));
                   [l(:,ii), lI] = sort(l(:,ii));
                   s(:,:,ii) = s(lI,tI,ii);
                   
                   % Updated the griddedInterpolant object with new vectors and data
                   % Update grid vectors, values, and evaluate on new grid
                   F.GridVectors = {l(:,ii), t(:,ii)};    %set interpolant grid
                   F.Values = s(:,:,ii);                  %set interpolant values
                   
                   %Update with user defined interpolation and
                   %extrapolation (required here for MATLAB 2020a)
                   F.Method = p.Results.interp;
                   F.ExtrapolationMethod = p.Results.extrap;
                   
                   % Define the new interpolation grid and evalaute interpolant
                   
                   % Define a new delay axis for each wavelength as a matrix first
                   tShift = t(:,ii)' + polyval(obj(objInd).chirpParams, l(:,ii))-tRef;
                   % Define the corresponding evaluation wavelength matrix
                   lEval = explicitExpand(l(:,ii),size(tShift));
                   
                   % Evaluate the interpolant on the (x,y) coordinate matrix by 
                   % reshaping tShift and lEval into column vectors, concatonating,
                   % and followed by reshaping back into pixel x delay
                   s(:,:,ii) = reshape(F([lEval(:), tShift(:)]),obj(objInd).sizes.nPixels, obj(objInd).sizes.nDelays);
                   
                   % Unsort the sorted dims
                   s(lI,tI,ii) = s(:,:,ii);  
               end
               
               % After for loop, undo reshape operation
               % s old: [pixels, delays, rpts x g pos x schemes]
               % s new: [pixels, delays, rpts, g pos, schemes]
               s = reshape(s, obj(objInd).sizes.nPixels, obj(objInd).sizes.nDelays, obj(objInd).sizes.nRpts, obj(objInd).sizes.nGPos, obj(objInd).sizes.nSchemes);

               % Assign interpolated spectra back to object data
               obj(objInd).spectra.data = s;
               
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
%   'fitFun': (function handle) A custom function handle to fit the t0 with
%       of the form: t0 = f(x,y). The default is the third element of 
%       lsqFitSigmoid.
%   'showPlot': (logical) Display fit result plots. Default is true.
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
           p.addParameter('fitFun',@(x,y) sum(lsqFitSigmoid(x, y).*[0,0,1,0]), @(p) isa(p,'function_handle'));
           p.addParameter('showPlot',true,@(l) islogical(l));
           
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
                   
                   % Runs the fit function. By default this is lsqFitSigmoid, 
                   % which fits data to y = a+b*0.5*(1+erf(sqrt(0.5)*(x-x0)/s));
                   % and returns the third element of fp, [a, b, x0, s]
                   myFP(wlInd,:) = p.Results.fitFun(t, data(wlInd,:));
                   
                   % Do the least squares fit using the lsqsigmoid fit
                   %myFP(wlInd,:) = lsqFitSigmoid(t, data(wlInd,:));
               end
               
               % Once the sigmoid fit is done, find the chirp parameters
               % and update object data
               chirpFit(:,objInd) = polyfit(wl,myFP(:,3),p.Results.order);
               obj(objInd).chirpParams = chirpFit(:,objInd);
               
               % Plot resuls
               if p.Results.showPlot
                   figure;
                        contour(wl,t,data',25);
                        hold on;
                        p1 = plot(wl,myFP(:,3),'r','DisplayName','sigmoid t0');
                        p2 = plot(wl,polyval(chirpFit(:,objInd),wl),'k--','DisplayName','polynomial');
                        hold off;

                   xlabel('Wavelength (nm)');
                   ylabel('Delay (ps)');
                   legend([p1,p2]);
                   colorbar();
               end
           end
           
           close(f);
           
           %convert object array back to original size
           obj = reshape(obj,objSize);
           chirpFit = reshape(chirpFit, [p.Results.order+1,objSize]);
       end
       
       function [obj, phononFits] = removePhonons(obj, varargin)
% REMOVEPHONONS fits and removes CAWs oscillations in the spectra.
% This method stores the fit results as a STOCAWs object, which contains
% the model published in: https://doi.org/10.1021/jacs.1c04976
%
% This method works by:
% 1. Trimming each spectrum to a representative wavelength/delay range free
%   of artifact
% 2. Using SVD to remove the first two components from each spectrum,
%   thereby leaving the CAWs oscillations
% 3. Using a global nonlinear fit to fit a trimmed version of the remaining
%   spectrum to a STOCAWs model
% 4. Evaluating the best fit STOCAWs model on the untrimmed wavelength and
%   delays and subtracting the model from the full spectrum
%
% The six free parameters in the fit are strain magnitude, spatial extent, 
% formation time, an exponential coherence decay, and a 2-parameter linear
% correction to the wavelength axis.
%
% Notes: 
% -For the best fit, run this method right after pruning the data. 
% -Do not perform any interpolation (e.g. interp, stitch, correctChirp) or 
% smoothing of the data before running this method. Interpolation and 
% smoothing can destory CAWs phase information. 
% -This method fits grating positions seperately, but averages over repeats
% before fitting. The returned object will still contain its original
% repeats and grating positions.
% -This method requires an accurate chirp fit in the object data (without 
% running chirp correction). Use setChirp or fitChirp to set the parameters
% -CAWs early-time dynamics are approximated by params 'fwhm', 'xiGuess', 
% and 't0Guess'. See below and STOCAWs.M for more detail.
%
% [obj, phononFits] = obj.removePhonons(varargin)
%   Fits and subtracts the CAWs phonon contribution for each object element
%   in obj. This method returns the phonon removed spectra in obj as well
%   as their best-fit STOCAWs object, which contains all fit paramters
%   except the wavelength correction.
%   
% This method can be called with name-value pair options:
%   'wlRange' [2,1] wavelength sub-range in nm to fit. Default is [375,700]
%   'tRange' [2,1] delay sub-range in psto fit. Default is [25,1000]
%   't0Fun' a custom function handle of the form t0 = f(x,y) for
%       determining t0. The default calls lsqFitExpGrow.
%   'wlRef' (scalar) the reference wavelength in nm for determining t0 and 
%       the t = 0 ps wavelength of the chirp polynomial. Default is 540 nm.
%   'fwhm' (scalar) the pump intensity fwhm duration in ps, default is
%       0.4 ps.
%   't0Guess' (scalar) guess for the CAWs formation time, default is 1.2 ps
%   'xiGuess' (scalar) guess for the CAWs spatial extent, default is 12 nm
%   'showPlots' (logical) Display diagnostic plots. Default is true.
%   'nStarts' (scalar) Number of starts for multistart. Default is 5.
%
% See Also: STOCAWs, fitChirp, setChirp, lsqFitExpGrow

           % Format object array dims into a column for easy looping
           objSize = size(obj);
           objNumel = numel(obj);
           obj = obj(:);
           
           % parse inputs
           p = inputParser();
           p.FunctionName = 'removePhonons';
           p.addParameter('wlRange', [375,700], @(d) isvector(d) && numel(d)==2);
           p.addParameter('tRange', [25,1000], @(d) isvector(d) && numel(d)==2);
           p.addParameter('t0Fun', @(x,y) sum(lsqFitExpGrow(x, y).*[0,0,0,0,0,1]), @(p) isa(p,'function_handle'));
           p.addParameter('wlRef',540,@(d) isscalar(d));
           p.addParameter('fwhm',0.4,@(d) isscalar(d) && d>=0);
           p.addParameter('t0Guess',1.2,@(d) isscalar(d) && d>=0);
           p.addParameter('xiGuess',12,@(d) isscalar(d) && d>=0);
           p.addParameter('showPlots', true, @(l) islogical(l) && isscalar(l));
           p.addParameter('nStarts', 5, @(d) isscalar(d));
           p.parse(varargin{:});
           
           % initialize output
           phononFits = cell(objNumel,1);
           
           % assert the same scheme exists
           schInd = 1;
           
           % Start waitbar
           [~,idx] = obj.getUniqueLabels('gPos');
           nSets = sum(idx,'all');
           ctr = 0;
           f = waitbar(ctr/nSets, sprintf('Removing phonons for gPos %d of %d for object %d of %d.',0,obj(1).sizes.nGPos,0,objNumel));
           
           % Loop over objects
           for objInd = 1:objNumel
               %remember old units
               tmpUnits = cell(3,1);
               [tmpUnits{:}] = obj(objInd).getUnits();
               
               % initialize suboutput
               phononFits{objInd} = cell(obj(objInd).sizes.nGPos,1);
               
               % Correct t0 at known position to pass onto fit function
               subObj = obj(objInd).average('rpts');
               [subObj, t0Ar] = subObj.findT0('wlAr',p.Results.wlRef,'fitFun',p.Results.t0Fun);
               subObj = subObj.correctT0();
               
               % trim and average over repeats to fit a subset of data
               subTrmd = subObj.trim('wavelengths',p.Results.wlRange,'delays',p.Results.tRange);
               
               % loop over grating positions
               for gInd = 1:subObj.sizes.nGPos
                   % Update waitbar
                   waitbar(ctr/nSets, f, sprintf('Removing phonons for gPos %d of %d for object %d of %d.', gInd,subObj.sizes.nGPos, objInd, objNumel));
                   ctr = ctr+1;
                   
                   % get numeric spectra trimmed data for the grating position
                   % getNumeric also removes NaN values, important for SVD
                   [s,l,t] = getNumeric(subTrmd.getLabel('gPos',gInd));
                                     
                   % Use SVD to remove 1st two components--removes emissive/absorptive dynamics and 
                   % keeps phonon oscillations
                   [U,S,V] = svd(s);
                   M = s-U(:,1:2)*S(1:2,1:2)*V(:,1:2)';
                   
                   % Phonon model OOP constructor
                   fm = STOCAWs(l,t,subTrmd.chirpParams,p.Results.wlRef,p.Results.fwhm);
                   
                   % Set guess parameters for spatial extent and formation time
                   xg = p.Results.xiGuess;
                   t0g = p.Results.t0Guess;
                   
                   % Set early-time CAWs formation parameters as guess parameters
                   fm.xiGuess = xg;
                   fm.t0Guess = t0g;
                   
                   % Setup initial guess fit parameters and bounds. l0 and l1 are a 
                   % linear correction to the spectrometer calibration. eta0 is the 
                   % strain amplitude, xi is the spatial extent, t0 is the formation time,
                   % and tau is a CAWs coherence decay constant modelling spectrometer
                   % resolution and acoustic damping.
                             %l0,  l1,  eta0,  xi,  t0,  tau
                   myGuess = [0,   1,   1e-3,  xg,  t0g, 1000];
                   typicalX= [1,   1,   1e-3,  5,   1,   500];
                   myLB =    [-10, 0.9, 0,     0,   0,   0];
                   myUB =    [10,  1.1, 0.1,   100, 10,  1e6]; 
                   
                   % Setup MultiStart optimization with lsqnonlin and tight convergence criteria
                   opts = optimoptions(@lsqnonlin,'FunctionTolerance',1e-9,'MaxFunctionEvaluations',2000,'MaxIterations',2000,'Display','off','TypicalX',typicalX);
                    
                   problem = createOptimProblem('lsqnonlin',...
                       'objective', @(fp) fm.evalMslp(fp(1)+fp(2)*l, fp(3), fp(4), fp(5), fp(6))-M,...
                       'x0',myGuess,'lb',myLB,'ub',myUB,'options',opts);
                   
                   % Run MultiStart solver
                   myFP = run(MultiStart,problem,p.Results.nStarts);
                   
                   % Update fitModel to the best solution
                   [~,fm] = fm.evalMslp(myFP(1)+myFP(2)*l, myFP(3), myFP(4), myFP(5), myFP(6));
                   
                   % Update phononFits
                   phononFits{objInd}{gInd} = fm;
                   
                   % Remove phonons and update object data
                   for rInd = 1:obj(objInd).sizes.nRpts
                      % Correct both wavelengths and t0 in obj
                      obj(objInd).spectra.data(:,:,rInd,gInd,schInd) = ...
                                obj(objInd).spectra.data(:,:,rInd,gInd,schInd) - ...
                                fm.evalMlt(myFP(1)+myFP(2)*obj(objInd).wavelengths.data(:,gInd),...
                                    obj(objInd).delays.data(:,rInd,gInd)-t0Ar{1}(:,:,gInd));
                   end
                   
                   % Plot results
                   if p.Results.showPlots
                       % get numeric spectra untrimmed data for current grating position
                       sFull = subObj.spectra.data(:,:,1,gInd,1);
                       lFull = subObj.wavelengths.data(:,gInd);
                       tFull = subObj.delays.data(:,1,gInd);

                       % Create full fitModel with full wavelengths and delays
                       [~,fmFull] = fm.evalMlt(myFP(1)+myFP(2)*lFull,tFull);                  
                       MFull = fmFull.M;
                       
                       % Plot contour of results                     
                       figure('Name',sprintf('%s, gPos: %d',subObj.shortName,subObj.gPos(gInd)));
                       subplot(1,4,1)
                           [~,c] = contour(lFull,tFull,sFull',23);
                           title('Raw');
                           xlabel(subObj.wavelengths.dispName);
                           ylabel(subObj.delays.dispName);
                           
                       subplot(1,4,3)
                           contour(lFull,tFull,sFull'-MFull',c.LevelList);
                           title('Raw - Phonons');
                           xlabel(subObj.wavelengths.dispName);
                           
                      subplot(1,4,2)
                           [~,c] = contour(lFull,tFull,MFull',23);
                           title('Phonons');
                           xlabel(subObj.wavelengths.dispName);
                           
                       subplot(1,4,4)
                           contour(l,t,fm.M'-M',c.LevelList);
                           title('Residuals');
                           xlabel(subObj.wavelengths.dispName);
                           ylabel(subObj.delays.dispName);
                   end
                   
                   % update waitbar
                   waitbar(ctr/nSets, f);
               end
               
               % set units back to input units
               obj(objInd) = obj(objInd).setUnits(tmpUnits{:});
           end
           
           % close waitbar
           close(f);
           
           %convert object array back to original size
           obj = reshape(obj,objSize);
       end

   end
end