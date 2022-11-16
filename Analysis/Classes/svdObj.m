classdef svdObj   
    properties (SetAccess=protected, GetAccess=public)
        Ur = 1; % raw spectral U matrix [nU, nU, nSets]
        Sr = 1; % raw weights S matrix [nU, nV, nSets]
        Vr = 1; % raw kinetic V matrix [nV, nV, nSets]
        
        Ut = 1; % constraining basis on U [nU, nC]
        Xr = 1; % raw constraining matrix for constrained SVD [nC, nC, nSets]
        nC = 2; % working number of components for constrained SVD and display
        
        flC = false(); % 'display' flag to flip U and V components [nU, nSets]
        pOrd = 0; % 'display' permute U and V component display order [nU, nSets]
    end
    
    properties (Dependent, SetAccess=protected, GetAccess=public)
        U;  % [nU, nC, nSets] 'display' U = Ur*T components 
        S;  % [nC, nC, nSets] 'display' S = T'*Sr*T weights 
        V;  % [nV, nC, nSets] 'display' V = Ur*T components 
        X;  % [nC, nC, nSets] 'display' X = T'*X constraining matrix 
        T;  % [nC, nC, nSets] 'display' transformation matrix 
        
        Uc;  % [nU, nC, nSets] constrained U*X components 
        SVc; % [nV, nC, nSets] constrained (inv(X)*S*V')' components 
        
        Sd;  % [nC, nSets] 'display' S weights diagonals 
        Srd; % [max(nU, nV), nSets] raw S weights diagonals 
        
        M;   % [nU, nV, nSets] reconstructed data matrix from first nC components 
        Mc;  % [nU, nV, nSets] reconstructed data matrix from constrained basis and first nC components 
        
        nU; % number of rows in M and components in U
        nS; % number of diagonal entries in S, min(nU,nV)
        nV; % number of cols in M and components in V
        nSets; % number of spectra sets in data
    end
    
    % Constructor and get/set methods
    methods
        function obj = svdObj(varargin)
% SVDOBJ performs a singular value decomposition (SVD) on an input spectra 
% data set and contains methods and dependent properties to display the 
% results in various forms. 
%
% obj = svdObj(data)
%   Builds a svdObj by performing SVD on numeric array data. Data is a
%   [m, n, nSets] matrix that calcualtes SVD over the m and n dimensions
%   for each spectrum over nSets. By default, assigns 2 components for
%   display and constraining.
%
% obj = svdObj(data, nC)
%   Same call as above except assigns nC components for display and
%   constraining.
%
% Properties:
%
% The main properties can be subdivided into three types:
%   1. raw components (suffixed with r)-Ur, Sr, Srd, and Vr that return 
%       all components as generated by MATLAB's [Ur,Sr,Vr] = svd(data) call
%   2. display components (no suffix)-U, S, Sd, and V that return the 
%       working number of components (default 2) modified by basis flip and
%       permute transformations (matrix T)
%   3. constrained components (suffixed with c)-Uc and SVc that return the
%       working number of components constrained to a target basis (default
%        is the 1st spectrum) with constraint transformation matrix X
%
% Additional properties:
%
% The class contains properties that return various transformation 
% matricies related to component display and constraining:
%   X, Xr, Ut (constraining basis), and T
%
% The class also reconstructs the initial data set from raw, display, and
% constrained components with methods and properties:
%   M, Mr(nC), Mc
%
% Finally, properties prefixed with n indicate the number of components, 
% sizes of decomposed matricies, and the number of spectra in the data set:
%   nC, nU, nS, nV, and nSets
%
% Note: all properties are read-only. To change the number of components,
% call the obj.setC(newNC) method.
%
% Methods:
%
% The class contains methods to sign flip and permute the columns of 
% U, S, and V for the purpose of comparing components across the spectra:
%   obj.flipC(__), obj.permuteC(__), and obj.resetC()
%
% The class also contains preview methods for plotting SVD components that
% mirror the various display forms:
%   obj.previewr(__), obj.preview(__), and obj.previewc(__)
%
% Finally, the class contains analysis methods that perform the SVD
% rotation analysis and constrained SVD:
%   obj.rotAnalysis() and obj.constrain(basis)
%
% See Also: svd

            % Parse inputs
            p = inputParser();
            p.addRequired('data');
            p.addOptional('nC',2);
            p.parse(varargin{:});
           
            % SVD matrix sizes
            nU = size(p.Results.data, 1); %number of wavelengths
            nV = size(p.Results.data, 2); %number of delays
            nSets = size(p.Results.data, 3); %number of spectra in set
            
            % Initialize svd component matricies
            obj.Ur = zeros(nU,nU,nSets);
            obj.Sr = zeros(nU,nV,nSets);
            obj.Vr = zeros(nV,nV,nSets);
            
            % Start waitbar
            w = waitbar(0,['Calculating SVD for spectra 1 of ' num2str(nSets)]);
            
            % Loop over sets and calculate SVD
            for sInd = 1:nSets
                % Update waitbar
                waitbar((sInd-1)/nSets,w,['Calculating SVD for spectra ' num2str(sInd) ' of ' num2str(nSets)]);
                
                % do raw SVD
                [obj.Ur(:,:,sInd),obj.Sr(:,:,sInd),obj.Vr(:,:,sInd)] = svd(p.Results.data(:,:,sInd));
                
                % Update waitbar
                waitbar(sInd/nSets,w);
            end
            
            % Close waitbar
            close(w);
            
            % Set the number of components
            obj.nC = p.Results.nC;
            
            % Set flip and permute orders
            obj = obj.resetC();
            
            % Calculate constrained SVD using 1st spectrum as constraining set
            obj = obj.constrain(obj.Ur(:,1:obj.nC,1));
            
        end
        
        function obj = setC(obj, nC)
% SETC updates the display and constrain component number and triggers an 
% update of the constraining basis. If the display permute order contains
% indicies outside of the component number, this method returns a warning
% and resets the permute order.
%
% obj = obj.setC(nC)
%   Updates the number of components in obj.
%
% See Also: constrain, flipC, permuteC, preview
            
            % Make sure the display component number does not exceed the
            % number of singular values
            assert(nC <= obj.nS, 'The component number cannot exceed the number of singular values, obj.nS = %d.', obj.nS);
            
            % Update component number and remember old component number
            oldnC = obj.nC;
            obj.nC = nC;
            
            % Update constraining basis only if the constraining basis does not match the size. 
            if nC ~= size(obj.Ut,2)
                if nC > oldnC  % when there are more new components than old components
                    newUt = obj.Ur(:,1:nC,1); % new cosntraining basis components are from the 1st nC cols of 1st spectrum Ur 
                    newUt(:,1:oldnC,1) = obj.Ut; % keep the old constraining basis
                else % when there are less or the same abount of components compared to old components
                    newUt = obj.Ut(:,1:nC); % keep a subset of the old constraining basis
                end
                
                % Update constrain
                obj = obj.constrain(newUt);
            end
            
            % Check to make sure pOrd does not contain permute indicies outside of obj.nC
            % This may be okay, just did not have time to test this use case
            inds = obj.pOrd(1:nC,:) > nC;   % permute indicies outside of obj.nC
            if any(inds, 'all')
                warning('The permute order obj.pOrd contains indicies outside of the range of obj.nC = %d. The permute order will be reset for offending spectra.', nC)
                setInds = any(inds,1);  %spectra that are offending
                obj.pOrd(:,setInds) = repmat((1:obj.nU)',1,sum(setInds)); %initialize offending permute orders to consecutive
            end
        end
        
        function U = get.U(obj)
        % Returns a display U, where U includes flipping and swapping components
            U = obj.Ur(:,1:obj.nC,:);
            T = obj.T;
            for i = 1:obj.nSets
               U(:,:,i) = U(:,:,i)*T(:,:,i); 
            end
        end
        
        function S = get.S(obj)
        % Returns a display S, where S includes flipping and swapping components
            S = obj.Sr(1:obj.nC,1:obj.nC,:);
            T = obj.T;
            for i = 1:obj.nSets
               S(:,:,i) = T(:,:,i)'*S(:,:,i)*T(:,:,i); 
            end
        end
        
        function V = get.V(obj)
        % Returns a display V, where V includes flipping and swapping components
            V = obj.Vr(:,1:obj.nC,:);
            T = obj.T;
            for i = 1:obj.nSets
               V(:,:,i) = V(:,:,i)*T(:,:,i); 
            end
        end
        
        function X = get.X(obj)
        % Returns a display X, where X includes flipping and swapping components
            X = obj.Xr;
            T = obj.T;
            for i = 1:obj.nSets
                X(:,:,i) = T(:,:,i)'*X(:,:,i);
            end
        end
        
        function T = get.T(obj)
        % Returns the display transformation matrix that flips and swaps components
            T = zeros(obj.nC,obj.nC,obj.nSets);
            for i = 1:obj.nSets
                T(:,:,i) = tMat(obj.flC(1:obj.nC,i),obj.pOrd(1:obj.nC,i));
            end
        end
        
        function Uc = get.Uc(obj)
        % Returns a constrained U, where the constraints were set by obj.constrain(basis)    
            Uc = zeros(obj.nU,obj.nC,obj.nSets);
            for i = 1:obj.nSets
                Uc(:,:,i) = obj.Ur(:,1:obj.nC,i)*obj.Xr(:,:,i);
            end
        end
        
        function SVc = get.SVc(obj)
        % Returns a constrained U, where the constraints were set by obj.constrain(basis)    
            SVc = zeros(obj.nV,obj.nC,obj.nSets);
            for i = 1:obj.nSets
                SVc(:,:,i) = (obj.Xr(:,:,i)\obj.Sr(1:obj.nC,1:obj.nC,i)*obj.Vr(:,1:obj.nC,i)')';
            end
        end
        
        function Sd = get.Sd(obj)
        % Returns the diagonal entries of display S
            Sd = zeros(obj.nC, obj.nSets);
            tmp = obj.S;    % to avoid excessive dependent property calls
            for i = 1:obj.nSets
               Sd(:,i) = diag(tmp(:,:,i)); 
            end
        end
        
        function Srd = get.Srd(obj)
        % Returns the diagonal entries of raw S
            Srd = zeros(obj.nS, obj.nSets);
            for i = 1:obj.nSets
               Srd(:,i) = diag(obj.Sr(:,:,i)); 
            end
        end
        
        function Mr = Mr(obj, nC)
% MR calculates the SVD reconstructed spectra from the specified number of 
% raw components.
%
% Mr = obj.Mr(nC)
%   Returns a spectra set Mr [nU, nV, nSets] reconstructed from the first
%   nC number of raw SVD components.

            % Make sure nC is within the available number of components
            assert(nC > 0 && nC < max(obj.nU, obj.nV),'The number of reconstruction components must be between 1 and the maximum number of row or column dims, %f.', max(obj.nU, obj.nV));
            
            % initialize variables
            Mr = zeros(obj.nU, obj.nV, obj.nSets);
            
            % Loop over spectral sets and calculate Mc
            for i = 1:obj.nSets
                Mr(:,:,i) = obj.Ur(:,1:min(nC,obj.nU),i)*obj.Sr(1:min(nC,obj.nU),1:min(nC,obj.nV),i)*obj.Vr(:,1:min(nC,obj.nV),i)';
            end
            
        end
        
        function M = get.M(obj)
        % Reconstruct the spectral data set from the display SVD components.

            % initialize variables
            M = zeros(obj.nU, obj.nV, obj.nSets);
            
             % copy dependent variables to avoid multiple calculations
            U = obj.U;
            S = obj.S;
            V = obj.V;
            
            % Loop over spectral sets and calculate Mc
            for i = 1:obj.nSets
                M(:,:,i) = U(:,:,i)*S(:,:,i)*V(:,:,i)';
            end
        end
        
        function Mc = get.Mc(obj)
        % Reconstruct the spectral data set from the constrained SVD components

            % initialize variables
            Mc = zeros(obj.nU, obj.nV, obj.nSets);
            
            % copy dependent variables to avoid multiple calculations
            X = obj.X;
            S = obj.S;
            V = obj.V;
            
            % Loop over spectral sets and calculate Mc
            for i = 1:obj.nSets
                Mc(:,:,i) = obj.Ut*(X(:,:,i)\S(:,:,i))*V(:,:,i)';
            end
        end
        
        function nU = get.nU(obj)
        % Returns the number of number wavelengths or elements in a U component
           nU = size(obj.Ur,1); 
        end
        
        function nS = get.nS(obj)
        % Returns the number of number diagonal elements in Sr
            nS = min(obj.nU, obj.nV);
        end
        
        function nV = get.nV(obj)
        % Returns the number of number delays or elements in a V component
           nV = size(obj.Vr,1); 
        end
        
        function nSets = get.nSets(obj)
        % Returns the number of spectra sets in obj
           nSets = size(obj.Ur,3);
        end
    end
    
    % Methods to manipulate display properties U, S, and V
    methods
        function obj = flipC(obj, cInd, dataSetInds, compLogical)
% FLIPC flags SVD spectral and kinetic components for flipping. This method
% does not affect the raw SVD components, only the display of components 
% retrieved by obj.U, obj.V, obj.subU(), and obj.subV().
%
% obj = obj.flipC(compInd, dataSetInds, compLogical)
%   Assigns a flip flag to the subset of components specified by compInd
%   and dataSetInds. The flag compLogical values are false = no flip and
%   true = flip.
%
% Example:
%   obj = obj.flipC(2,[1,3,5],true); flags the 2nd spectral component of 
%   the 1st, 2nd, and 3rd spectrum for flipping (multiply by -1) upon 
%   display. This is equivalent of multiplying the 1st two components by 
%   the expanded identity in both US and SV':
%   M = (U12*[1  0])*([1  0]*S12*[1  0])*([1  0]*V12')
%            [0 -1]   [0 -1]     [0 -1]   [0 -1]
%
% See Also: permuteC, preview, resetC, setC
            
            if islogical(dataSetInds)
                assert(numel(dataSetInds) == obj.nSets, 'The number of data set indicies must match the number of spectra, obj.nSets = %d.', obj.nSets);
            else
                assert(all(dataSetInds<=obj.nSets) && all(dataSetInds>0), 'Data set indicies must be within 1 and %d.', obj.nSets);
            end
            assert(cInd>0 && cInd<=obj.nU, 'Component index cInd must be within 1 and %d.', obj.nU);
            
            obj.flC(cInd,dataSetInds) = compLogical;
        end
        
        function obj = permuteC(obj, dataSetInds, permuteOrder)
% PERMUTEC flags SVD spectral and kinetic components for order permutation.
% This method does not affect the raw SVD components, only the display of
% components retrieved by obj.U, obj.V, obj.subU(), and obj.subV().
%
% obj = obj.permuteC(dataSetInds, permuteOrder)
%   Assigns a permute order to the first subset of components specified by 
%   permuteOrder for spectrum numbers given by dataSetInds.
%
% Example:
%   obj = obj.permuteC(1,[2,1]); flags the 2nd component to display first
%   and 1st component to display 2nd for the 1st spectrum in the data set.
%   This is equivalent of multiplying the 1st two components by the
%   expanded identity in both US and SV':
%   M = (U12*[0 1])*([0 1]*S12*[0 1])*([0 1]*V12')
%            [1 0]   [1 0]     [1 0]   [1 0]
% 
% See Also: flipC, resetC, preview

            n = numel(permuteOrder); %number of permute inds to change
            
            % Make sure that the permute order does not include invalid indicies 
            assert(all(sort(permuteOrder(:))==(1:n)'),'Invalid permute order.');
            assert(max(permuteOrder) <= obj.nC,'Permute order cannot include indicies outside of the number of display components, obj.nC = %d.', obj.nC);
            assert(all(dataSetInds<=obj.nSets) && all(dataSetInds>0),'Data set indicies must be within 1 and %d.', obj.nSets);
            
            % Assign permute order
            for i = 1:numel(dataSetInds)
                obj.pOrd(1:n,dataSetInds(i)) = permuteOrder(:);
            end
        end
        
        function obj = resetC(obj)
% RESETC removes any user flip and permute assignments. This method does not
% affect the raw SVD components, only the display of components retrieved 
% by obj.U, obj.V, obj.subU(), and obj.subV().
%
% obj = obj.resetC();
%   Reset flip and permute assignemtns.
%
% See Also: flipC, permuteC
            
            obj.flC = false(obj.nU, obj.nSets);    %initialize flip flag to false
            obj.pOrd = repmat((1:obj.nU)',1,obj.nSets); %initialize permute order to consecutive
        end
    end
    
    % Methods to display/preview components
    methods
        function previewc(obj, varargin)
% PREVIEWC generates a figure with constrained components.
% The components are plotted as a 2 x nC subplot, where the top row 
% contains the first nC components of U*X and the bottom row contains the
% first nC components of inv(X)*S*V'.
%
% [] = obj.previewc()
%   Generates a plot starting from the current figure to display the 
%   first nC constrained components of U*X and inv(X)*S*V'. By default, 
%   the x-axis values are indicies and the legend values are indicies.
%
% [] = obj.previewc('name1',value1,'name2',value2,...)
%   The same call as above but with additional name-value pairs:
%   
%    'xLabel' (char) x-axis label for the U*X components 
%    'yLabel' (char) x-axis label for the inv(X)*S*V' components
%    'xVals', (vector) x-axis values for the U*X components
%    'yVals', (vector) x-axis values for the inv(X)*S*V' components
%    'legend', (cell of char) legend labels for the spectra in the data set
%
% See Also: constrain, preview, previewr, rotAnalysis
            
            % generate cosmetic display parameters inside an input parser object
            defLegendVals = cellfun(@(n) num2str(n),num2cell(1:obj.nSets),'uniformOutput',false);
            p = parsePlotInputs(1:obj.nU, 1:obj.nV, defLegendVals, varargin{:});
            
            % generate figure with subplot
            plotSVD(obj.Uc, obj.SVc, p.Results.xVals, p.Results.yVals, p.Results.xLabel, p.Results.yLabel, p.Results.legend,'');
            
            % Add basis components to U plots
            for i = 1:obj.nC
                subplot(2,obj.nC,i)
                hold on;
                plot(p.Results.xVals,obj.Ut(:,i),'--k','DisplayName','basis');
                hold off;
            end
        end
        
        function preview(obj, varargin)
% PREVIEWC generates a figure with display components of U and V.
% The components are plotted as a 2 x nC subplot, where the top row 
% contains the first nC components of U and the bottom row contains the
% first nC components of V.
%
% [] = obj.preview()
%   Generates a plot starting from the current figure to display the 
%   first nC constrained components of U and V. By default, the x-axis 
%   values are indicies and the legend values are indicies.
%
% [] = obj.preview('name1',value1,'name2',value2,...)
%   The same call as above but with additional name-value pairs:
%   
%    'xLabel' (char) x-axis label for the U components 
%    'yLabel' (char) x-axis label for the V components
%    'xVals', (vector) x-axis values for the U components
%    'yVals', (vector) x-axis values for the V components
%    'legend', (cell of char) legend labels for the spectra in the data set
%
% See Also: flipC, permuteC, previewc, previewr, resetC, setnC

            % generate cosmetic display parameters inside an input parser object
            defLegendVals = cellfun(@(n) num2str(n),num2cell(1:obj.nSets),'uniformOutput',false);
            p = parsePlotInputs(1:obj.nU, 1:obj.nV, defLegendVals, varargin{:});
            
            % generate figure with subplot
            plotSVD(obj.U, obj.V, p.Results.xVals, p.Results.yVals, p.Results.xLabel, p.Results.yLabel, p.Results.legend,'');
        end
        
        function previewr(obj, nC, varargin)
% PREVIEWR generates a figure with constrained components.
% The components are plotted as a 2 x nC subplot, where the top row 
% contains the first nC components of Ur and the bottom row contains the
% first nC components of Vr.
%
% [] = obj.previewr(nC)
%   Generates a plot starting from the current figure to display the 
%   first nC constrained components of Ur and Vr. By default, the x-axis 
%   values are indicies and the legend values are indicies.
%
% [] = obj.previewr(__,'name1',value1,'name2',value2,...)
%   The same call as above but with additional name-value pairs:
%   
%    'xLabel' (char) x-axis label for the Ur components 
%    'yLabel' (char) x-axis label for the Vr components
%    'xVals', (vector) x-axis values for the Ur components
%    'yVals', (vector) x-axis values for the Vr components
%    'legend', (cell of char) legend labels for the spectra in the data set
%
% See Also: preview, previewc

            % generate cosmetic display parameters inside an input parser object
            defLegendVals = cellfun(@(n) num2str(n),num2cell(1:obj.nSets),'uniformOutput',false);
            p = parsePlotInputs(1:obj.nU, 1:obj.nV, defLegendVals, varargin{:});
            
            % generate figure with subplot
            plotSVD(obj.Ur(:,1:nC,:), obj.Vr(:,1:nC,:), p.Results.xVals, p.Results.yVals, p.Results.xLabel, p.Results.yLabel, p.Results.legend,'');
        end
    end
    
    % Data analysis methods
    methods
        function obj = constrain(obj, Ut)
% CONSTRAIN performs a constrained SVD of the full SVD data set onto a
% specified basis. This method calculates the transformation matricies
% (obj.Xr) to transform the raw basis into the target constrained basis. 
% Use dependent properties obj.cU and obj.cSV to return the constrained 
% U and V components, respectively. Use obj.X to return the display 
% constraint matrix, given by M = U(T'X)(inv(X)T)SV'), where T is the 
% display transform matrix.
%
% obj = constrain(obj, Ut)
%   Constrains the Ur components in obj to Ut by calculating Xr. Ut needs
%   to be [nU, nC]. If the number of cols in Ut is different then Ut, this
%   method will display a warning and update obj.nC by calling 
%   obj.setC(nC).
%
% The Nature Materials paper (https://doi.org/10.1038/s41563-021-01118-9) 
% uses the following conventions:
%   i (n in the paper) represents the U components of the ith data set
%   t represents the constraining basis
%   M is a m x n spectrum matrix 
%   U is a m x nC matrix
%   S is a nC x nC matrix
%   X is the constraint matrix nC x nC  
% 
% The constrained matrix Xi can be calculated according to:
%   Mi = UiSiVi'
%   Mi = (UiXi)(inv(X)SiVi') ~~ Ut(inv(X)SiVi')
%   Xi = (Ui)'Ut
% where Ui represents the raw SVD U components, Ut represents the 
% corresponding constraining basis, and ~~ is approx equal. 
%
% Note: For the purpose of display, this method does not normalize Xi as
% described in the Nature Materials SI. To get normalized Xi, normalize the
% basis vectors according to their vector norm, Ut --> Ut./vecnorm(Ut). 
% When the basis vectors are normalized, the columns of Xi should also be
% normalized. In practice, the columns of Xi are close-to-normalized
% because the constraining basis vectors Ut cannot be exactly derived from
% a linear combination of Ui: Ut ~~ Ui*Xi. diag(Xi'*Xi) provides a measure 
% of consistency between the constraining and spectral bases, where the
% values should be close to 1.
%
% Also See: previewc, setC, vecnorm
            
            % Ensure that the constraining basis matches the size requirements for the raw basis
            assert(size(Ut,1)==obj.nU,'The number of rows in the constraining basis must match the number of rows in obj.U');
            
            % Update basis in obj
            obj.Ut = Ut;
            
            % Update the number of display components to match the size of the constraining basis
            ntC = size(Ut,2);
            if ntC~=obj.nC
                warning('The number of display components, obj.nC = %d, was updated to match the number of basis components, Ut = %d',obj.nC,ntC);
                obj = obj.setC(ntC);
            end
                
            % initialize transformation matrix size
            obj.Xr = zeros(ntC,ntC,obj.nSets);
            
            % calculate transformation matrix for each data set
            for i = 1:obj.nSets
                obj.Xr(:,:,i) = obj.Ur(:,1:ntC,i)'*Ut;
            end
        end
        
        function [th, R] = rotAnalysis(obj)
% ROTANALYSIS performs a rotation analysis on the full SVD data set.
% This analysis calculates the transformation matrix from one spectrum to 
% another in the first 2 components and returns the rotation angle 
% represented by the transformation matrix. 
% 
% [th, R] = obj.rotAnalysis()
%   Performs a rotation analysis on the first two 'display' components 
%   inside obj. Returns th as the rotation angles in Rij (in radians) as a 
%   [nSets, nSets] matrix where i and j represet row and col, respectively.
%   Optionally returns the transformation matrix Rij as a [2, 2, nSets,
%   nSets] matrix. See below for more details on Rij.
%
% Note: To get the rotationally constrained bases, use the constrain method
% along with the desired basis. For example:
% obj = obj.constrain(obj.U(:,1:2,k)); constrains to the kth basis of U.
% Use obj.Uc and obj.SVc to return the constrained components and
% obj.previewc(__) to preview the constrained components.
%
% Note: this method uses the obj.U, obj.S, and obj.V display properties, 
% meaning that any flips or permutations are included in the outputs.
%
% The Nature Materials paper (https://doi.org/10.1038/s41563-021-01118-9) 
% uses the following conventions:
%   i and j represent two different spectra from one data set
%   i represents the reference basis
%   j represents the target basis
%   M is a m x n matrix
%   U is a m x 2 matrix
%   S is a 2 x 2 matrix
%   V is a n x 2 matrix
%   R is a 2 x 2 matrix
%
% The rotation matrix R transforms a target spectrum in terms of the
% reference basis as (~~ means aprox. equal):
%   Mj ~~ (Ui*inv(Rij))*Sj*Vj'
%   Ui ~~ Uj*Rij
%   Rij = Uj'*Ui ~~ [cos(thij) -sin(thij)]
%                   [sin(thij)  cos(thij)]
%
% If i and j exist on the same 2D plane in hyperspace, R is a rotation 
% matrix and can be interpreted as the amount of rotation required to 
% rotate the jth basis into the ith basis and spectra i and j can be said 
% to share a common basis.
%
% See Also: constrain, previewc

            obj = obj.setC(2); % set the number of components to 2
            
            % initialize outputs
            R = zeros(2,2,obj.nSets,obj.nSets);
            th = zeros(obj.nSets,obj.nSets);
            
            % loop over spectra sets 
            for i = 1:obj.nSets
                for j = 1:obj.nSets
                    % Calcualte R matrix
                    R(:,:,i,j) = obj.U(:,:,j)'*obj.U(:,:,i);

                    % Extract angle from R matrix
                    if sign(R(1,1,i,j))==R(2,2,i,j)
                        th(i,j) = atan2(R(2,1,i,j),R(1,1,i,j));
                    else 
                        Rtmp = R(:,:,i,j)*[1 0; 0 -1];
                        th(i,j) = atan2(Rtmp(2,1),Rtmp(1,1));
                    end
                end
            end
        end
    end
end

% Reuse routines in the class
function T = tMat(subFlC, subPOrd)
% Calculates the transformation matrix based on array subsets of flip and 
% permute order.
%
% T = tMat(subFlC, subPOrd)
%   Logical vector subFlC specifies if the component needs to be flipped.
%   Double vector subPOrd specifies the permute order of the columns
%   Transformations go as M = USV' = (UT)(T'ST)(VT)'
%   For example: if subFLc is [0,1] and subPOrd is [2,1]:
%   T = [0  1]
%       [-1 0]
%    
    % generate a diagonal matrix of 1 or -1
    f = ones(1,numel(subFlC));
    f(subFlC) = -1;
    T = diag(f);
    
    % permute col order
    T = T(:,subPOrd);
end

function plotSVD(U,V,xVals,yVals,xLabel,yLabel,legendVals,linespec)
% PLOTSVD generates a subplot of components of U (x) and V (y) along with
% plot display options like the U and V x-axis values and labels and the
% data set labels.

    nC = size(U,2);
    for i = 1:nC
        subplot(2,nC,i);   % U components
            p = plot(xVals,squeeze(U(:,i,:)),linespec);
            % update DisplayName property for legend
            for j = 1:numel(p)
               p(j).DisplayName = legendVals{j}; 
            end
            xlabel(xLabel);
            xlim([min(xVals),max(xVals)]);
            title(sprintf('U Comp. %s',num2str(i)));
        subplot(2,nC,nC+i); % V components
            p = plot(yVals,squeeze(V(:,i,:)),linespec);
            % update DisplayName property for legend
            for j = 1:numel(p)
               p(j).DisplayName = legendVals{j}; 
            end
            title(sprintf('V Comp. %s',num2str(i))); 
            xlabel(yLabel);
            xlim([min(yVals),max(yVals)]);
    end

end

function p = parsePlotInputs(defXVal, defYVal, defLegendVals, varargin)
% Parse display options for plotting SVD components. Input arguments are
% default values and extra arguments that will be parsed by an input parser
% scheme.

    p = inputParser();
    p.addParameter('xLabel','U elem.');
    p.addParameter('yLabel','V elem.');
    p.addParameter('legend',defLegendVals);
    p.addParameter('xVals',defXVal);
    p.addParameter('yVals',defYVal);
    p.parse(varargin{:});

end