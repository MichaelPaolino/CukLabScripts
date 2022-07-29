classdef wlTR < transientSpectra
   properties
       t0Shift = 0;
       chirpParams = 0;
       chirpCorrected = false;
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
       
       function [obj,chirpFit] = fitChirp(obj,wlRange,tRange,polyOrder)
           %--PREFORMAT DATA--
           % Format object array dims into a column for easy looping
           objSize = size(obj);
           objNumel = numel(obj);
           obj = obj(:);
           
           % Trim object size to desired sub-range
           objTmp = obj.trim('wavelengths',wlRange,'delays',tRange);
           
           % Set units to common values so that initial guess below would be valid
           objTmp = objTmp.setUnits('nm','ps','mOD');
           
           % Average repeats and stitch grating position
           objTmp = objTmp.average();
           objTmp = objTmp.stitch();
                      
           % Define sigmoid fit function that starts at a, rises to b, with
           % std s centered at x0
           mySigmoid = @(a,b,x0,s,x) a+b*0.5*(1+erf(sqrt(0.5)*(x-x0)/s));           
           
           % Define lsqnonlin options for tight convergance and quiet operation
           opts = optimoptions(@lsqnonlin,'TolFun',1e-9,'TypicalX',[0.2 5 0.5 0.2]','Display','off');
           
           % Initialize chirpFit output
           chirpFit = zeros(polyOrder+1, objNumel);
           
           % Start waitbar
           f = waitbar(0,['Fitting group delay for spectra 1 of ' num2str(objNumel)]);
           
           % Loop over individual object elements
           for objInd = 1:objNumel
               % Update waitbar
               waitbar(0,f,['Fitting group delay for spectra ' num2str(objInd) ' of ' num2str(objNumel)]);
               
               % Extract spectra data from object and use locally
               data = objTmp(objInd).spectra.data;
               wl = objTmp(objInd).wavelengths.data;
               t = objTmp(objInd).delays.data;
               
               % Initialize sigmoid fit parameter matrix
               myFP = zeros(length(wl),4);
               
               % Loop over wavelengths
               for wlInd = 1:length(wl)
                   % Update waitbar
                   waitbar(wlInd/length(wl),f);
                   
                   TF = isnan(data(wlInd,:));
                   dataTmp = data(wlInd,~TF);
                   tTmp = t(~TF);
                   
                   % determine a good initial guess for b and its bounds
                   minVal = min(dataTmp);
                   maxVal = max(dataTmp);
                   
                   if maxVal > abs(minVal)
                       bG = maxVal;
                       bGL = 0;
                       bGU = 1.5*maxVal;
                   else
                       bG = minVal;
                       bGL = 1.5*minVal;
                       bGU = 0;
                   end
                   
                   % Setup guess for least squares search
                   %          a  b    x0  s
                   myGuess = [0  bG   0   0.2]; %initial guess
                   myLB =    [-1 bGL  -2  0];   %lower bound
                   myUB =    [1  bGU  2   1];   %upper bound
                   
                   %do the least squares fit
                   myFP(wlInd,:) = lsqnonlin(@(fp) mySigmoid(fp(1),fp(2),fp(3),fp(4),tTmp(:))-dataTmp(:),myGuess,myLB,myUB,opts);
               end
               
               % Once the sigmoid fit is done, find the chirp parameters
               % and update object data
               chirpFit(:,objInd) = polyfit(wl,myFP(:,3),polyOrder);
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
           chirpFit = reshape(chirpFit, [polyOrder+1,objSize]);
       end
   end
end