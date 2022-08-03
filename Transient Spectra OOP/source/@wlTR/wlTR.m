classdef wlTR < transientSpectra
   properties
       t0Shift = 0;
       chirpParams = 0;
       chirpCorrected = false;
   end
   
   properties (Access = protected)
       TRflag = struct('chirpFit',false);     
   end
   
   methods
       
       function obj = correctChirp(obj, chirpPoly)
           % Format object array dims into a column for easy looping
           objSize = size(obj);
           objNumel = numel(obj);
           obj = obj(:);
           
           % Loop over object array
           for objInd = 1:objNumel
               
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
% obj = obj.CHIRPFIT()
%   Extracts the WLC group delay vs. wavelength and fits the data to a 5th
%   order polynomial.
%
% obj = obj.CHIRPFIT(varargin)
%   Extracts the WLC group delay vs. wavelengths and fits to a polynomail
%   with the additional name-value pair options:
%
% Name-value pairs:
%   'wavelengths': (double array) an array of wavelengths to subset the
%      spectra to help speed up fitChirp. Default is to use all
%      wavelengths.
%   'delays': (double array) [lower, upper] a 1x2 array of delays to trim
%      that limits the sigmoid fit range. Default is -2 to 2 ps.
%   'order': (int) an integer > 0 that specifies the polynomial order of
%      the group delay vs. wavelength. The default order is 5.
           
           %--PREFORMAT DATA--
           % Format object array dims into a column for easy looping
           objSize = size(obj);
           objNumel = numel(obj);
           obj = obj(:);
           
           % Format varargin using input parser
           p = inputParser;
           p.FunctionName = 'correctT0';
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
                   myFP(wlInd,:) =  lsqFitSigmoid(t, data(wlInd,:));
               end
               
               % Once the sigmoid fit is done, find the chirp parameters
               % and update object data
               chirpFit(:,objInd) = polyfit(wl,myFP(:,3),p.Results.order);
               obj(objInd).chirpParams = chirpFit(:,objInd);
               
               % Plot resuls
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
   end
end