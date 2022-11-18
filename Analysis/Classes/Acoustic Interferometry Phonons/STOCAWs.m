classdef STOCAWs
    properties (Constant)
        B = 170.57; % Bulk modulus in GPa
        G = 108.46; % Shear modulus in GPa
        R = 0.928;  % Acoustic amplitude reflection coefficient for STO/water interface
        rho = 5.11; % density in g/cm3
        cOptical = 299792.458; % speed of light in nm/ps
        cAcoustic = sqrt((STOCAWs.B+4/3*STOCAWs.G)*1e10/STOCAWs.rho)*1e-5; % acoustic velocity in nm/ps calcualted from elasticity constants
        
        p1122 = 0.095;   % elasto-optic coefficient at 633 nm
        th0 = 45/180*pi; % incident angle onto sample cell in radians
    end
    
    % Non-Static properties--Update when changing wavelength, read-only, visible to the user
    properties
        l;          % wavelengths [nl,1]--set by constructor
        t;          % delays [nt,1]--set by constructor
        cp;         % chirp polynomial [1,np]--set by constructor
        lt0;        % wavelength scalar that defines t0 for t and cp
        tau = inf;    % CAWs exponential attenuation caused by acoustic damping and spectrometer resolution
        
        xi = 15;    % spatial extent of transformation strain in nm--set by user
        t0 = 1e-3;  % formation time of transformation strain in ps--set by user
        fwhm = 0;   % fwhm duration of pump pulse in ps--set by user
        eta0 = 0;   % amplitude of transformation strain--set by user
        
        t0Guess = 1.2;  % Guess at t0 used to reconstruct early-time dynamics in obj.M
        xiGuess = 12.2; % Guess at xi used to reconstruct early-time dynamics in obj.M
        
        nH2O;  % refractive index of H2O evaluated on l [nl,1]
        nSTO;  % complex refractive index of STO evaluated on l [nl,1]
        
        thi;   % incident angle in water in radians evaluated on l [nl,1]
        tht;   % transmitted angle in STO in radians evaluated on l [nl,1]
        
        rs1;   % complex amplitude reflection water/STO evaluated on l [nl,1]
        ts1;   % complex amplitude transmission water/STO evaluated on l [nl,1]
        rs2;   % complex amplitude reflection STO/water evaluated on l [nl,1]
        ts2;   % complex amplitude transmistion STO/water evaluated on l [nl,1]
        
        deps;  % elasto-optic disperson dPerm/dStrain evaluated on l [nl,1]
        freq;  % acoustic frequency in GHz evaluated on l [nl,1]
    end
    
    % Protected Properties--Update before varying parameters, internal, not visible to the user
    properties (Access=protected)
       nElips;  % Actual 0.1% STO real refractive index measured by ellipsometry and loaded by constructor
       kUVVIS;  % Actual 0.1% STO kappa measured by UVVIS absorption and loaded by constructor
       
       % Sellmeier parameterization of H2O from: https://refractiveindex.info/?shelf=main&book=H2O&page=Daimon-21.5C
       nSelmH2O = @(lum) sqrt(1+5.689093832E-1./(1-5.110301794E-3./lum.^2)+...
                              1.719708856E-1./(1-1.825180155E-2./lum.^2)+...
                              2.062501582E-2./(1-2.624158904E-2./lum.^2)+...
                              1.123965424E-1./(1-1.067505178E1./lum.^2));
       
       etaConv = 3*STOCAWs.B*1e10/(STOCAWs.rho*(STOCAWs.cAcoustic*1e5)^2); %Conversion factor as defined in JACS paper in eq. 6, but without -eps0*xi/vSTO
       rdeps; % complex acoustic scattering prefactor evaluated on l [nl,1]
       depsl; % Acoustic pulse spectrum evaluated on acoustic wavelength [nl,1]
    end
    
    %Dependent Properties--Calculate when called, visible to the user
    properties (Dependent)
        amp;    % phonon oscillation amplitude in mOD evaluated on l [nl,1]
        phase;  % phonon oscillation phase in radians evaluated on l [nl,1]
        Ms;     % [l,t] phonon mOD spectrum evaluated on l and t without t0/causality [nl,nt]
        M;      % [l,t] phonon mOD spectrum considering pump fwhm and causality [nl,nt]
    end
        
    % Constructor and get/set methods
    methods
        function obj = STOCAWs(l, t, cp, lt0, fwhm)
% STOCAWs is a model object that calculates the CAWs contribution to TR
% spectra. Use this class to fit various aspects of the CAWs spectrum to
% experimental data. This class follows the model described in the 2021
% JACS paper: https://doi.org/10.1021/jacs.1c04976
%
% obj = STOCAWs(l, t, cp, lt0, fwhm)
%   Returns a STOCAWs object with wavelengths set to l, delay set to t,
%   chirp paramters set to cp, the chirp reference wavelength lt0, and pump
%   pulse duration intensity fwhm to fwhm
%
% This model contains the following properties:
% 1. Properties that describe elastic and optical constants of 0.1% SrTiO3
%    (B, G, rho, cOptical, cAcoustic, and p1122) as well as the 
%    experimental geometry (th0)
% 2. Properties that are free parameters for the model (eta0, xi, t0)
% 3. Properties that are required to reconstruct actual TR CAWs data
%    (l, t, sp, lt0, tau, and fwhm)
% 4. Intermediate properties that are calculated on l, including
%    interpolated refractive indicies (nH2O, nSTO), Snell's law angles (thi
%    and tht), Fresnel reflection/transmission coefficients (rs1, ts1, rs2,
%    ts2), elasto-optic effect (deps), and SBS phase matching frequency 
%   (freq)
% 5. Non-accessible protected properties
% 6. Dependent properties that are calculated on-the-fly, including the
%    acoustic amplitude spectrum evaluated on l, the acoustic phase
%    spectrum evaluated on l, and the TR CAWs spectra Ms and M evaluated 
%    on l and t.
%
% Note: this class uses a heuristic to calculate early-time formation of
% the CAWs signal in obj.M. The formation is parameterized by obj.fwhm,
% obj.t0guess and obj.xiguess. See obj.M documentation for more detail.
%
% This model includes methods that can be grouped into three classes:
% 1. Constructor + get/set methods to calculate dependent properties.
% 2. Methods that update and recalculate properties. These methods start
%    with set and are suffixed with the property they update, such as l for
%    obj.l and p for any object pulse parameters (eta0, xi, t0, tau, fwhm).
% 3. Methods that are designed to be used in fitting routines as function
%    handles and return calculated spectra. These methods start with eval
%    and are suffixed with the property that they return (e.g. M and Ms) 
%    and change (e.g. l, t, and p for pulse parameters). 
%    Todo: add evalAmpl, evalPhasel, evalAmplp, evalPhaselp, etc.
%
% See Also: lsqnonlin

            obj.l = l(:);    % set wavelengths
            obj.t = t(:);    % set delays
            obj.cp = cp;     % set chirp parameters
            obj.lt0 = lt0;   % set chirp t0 wavelength
            obj.fwhm = fwhm; % set pump pulse duration
            
            % Load STO refractive index data into object
            obj.nElips = load('0p1_ellipsometry.mat', 'el_n_fit');
            obj.kUVVIS = load('Abs_UVVis_0deg_aoi.mat');
            
            % Update any params that are l-dependent but do not need to be recalc'd
            obj = obj.setl(l);
        end
        
        function amp = get.amp(obj)
            amp = -1000*log10(1-2*abs(obj.ts1.*obj.ts2.*obj.rdeps.*obj.deps.*obj.depsl)./abs(obj.rs1));
        end
        
        function phase = get.phase(obj)
            phase = pi-angle(obj.rs1.*conj(obj.ts1.*obj.ts2.*obj.rdeps.*obj.deps.*obj.depsl));
        end
        
        function Ms = get.Ms(obj)
            % Chirp corrected wavelength dependent time [nl, nt]
            tChirp = obj.t'-polyval(obj.cp,obj.l)+polyval(obj.cp,obj.lt0);
            
            % Calcualted sub spectrum (without t0/causality)
            Ms = obj.amp.*exp(-(0.002*pi*imag(obj.freq)+1/obj.tau).*tChirp).*cos(0.002*pi*tChirp.*real(obj.freq)+obj.phase);
        end
        
        function M = get.M(obj)
% Returns calcualted CAWs spectrum (with t0/causality) [nl, nt]
%
% In reality, the true function should include early-time acoustic 
% interference and optical thin-layer interference effects that occur when 
% the acoustic pulse has not yet left the surface strain's spatial extent.
% The growth of the oscillations appears to be approx an exponential 
% growth in formation time t0 times a growth in the amount of time it takes
% the acoustic pulse to travel the spatial extent xi:
%   (1-exp(-t/t0)).*(1-exp(-t/xi*c_acoustic))
%
% Note: this is just a heuristic, this approx may not be valid in all cases
% See acousticPulse_numerical_full.m for more detail.
%
% Since the values of xi and t0 may not be accurate when fit, use 
% obj.xiGuess and obj.t0Guess to set the early-time CAWs behavior.
%
% If obj.fwhm is set to 0 or there are not enough points near t = 0 ps to
% calculate convolution, this function returns M without pump convolution
% effects, but with early-time CAWs approx effects.
            
            % Chirp corrected wavelength dependent time [nl, nt]
            tChirp = obj.t'-polyval(obj.cp,obj.l)+polyval(obj.cp,obj.lt0);
            
            % Pre-calculate Ms with pump convolution
            ht = hvsd(tChirp);
            M = ht.*(1-exp(-ht.*tChirp/obj.t0Guess)).*(1-exp(-ht.*tChirp.*obj.xiGuess/obj.cAcoustic)).*obj.Ms;
            
            % Check for non-zero or invalid pump duration
            if obj.fwhm > 0
                % Gaussian filter function
                f = 2.355/obj.fwhm/sqrt(2*pi)*exp(-(5:0.5:5).^2/2);
                dt = 0.5*obj.fwhm/2.355;

                % Select a sufficiently early sub-range for non-uniform convolution
                maxt = max(polyval(obj.cp,obj.l)-polyval(obj.cp,obj.lt0))+2*obj.fwhm;
                [~,maxti] = min(abs(obj.t-maxt));
                
                mint = min(polyval(obj.cp,obj.l)-polyval(obj.cp,obj.lt0))-2*obj.fwhm;
                [~,minti] = min(abs(obj.t-mint));
                
                % At least three points are required for convolution
                if (maxti-minti) >= 3
                    % Set obj fwhm to 0 because these effects will be included in nonuniform conv
                    obj.fwhm = 0;
                    
                    % re-calculate early-time M without pump convolution
                    tChirp = tChirp(:,minti:maxti);
                    ht = ht(:,minti:maxti);
                    M(:,minti:maxti) = ht.*(1-exp(-ht.*tChirp/obj.t0Guess)).*(1-exp(-ht.*tChirp.*obj.xiGuess/obj.cAcoustic)).*obj.Ms(:,minti:maxti);

                    % Loop over wavelengths and convolve the early-time data with pump pulse duration
                    for lInd = 1:size(M,1)
                        M(lInd,minti:maxti) = nonuniformConv(M(lInd,minti:maxti),tChirp(lInd,:),f,dt);
                    end
                end
            end
        end
    end
    
    % Set parameters and recalculate properties
    methods
        function obj = setl(obj,l)
% SETL updates the wavelengths in the object and triggers an update of
% properties dependent on l.
%
% obj = obj.setl(l)
%   Update l (nl,1) in nm and recalculate related properties.
%
% % See Also: SETP, SETLP

            % Update l in obj member data
            obj.l = l(:);
            
            % Real index of refraction of water
            obj.nH2O = obj.nSelmH2O(obj.l/1000);

            % Complex index of refraction of STO
            n1 = polyval(obj.nElips.el_n_fit, obj.l);
            n2 = interp1(obj.kUVVIS.nm, obj.kUVVIS.k, obj.l);
            obj.nSTO = n1+1i*n2;
            
            % Real AOI in water
            obj.thi = asin(1./obj.nH2O.*sin(obj.th0));
            
            % Calculate fresnel factors (reflection/transmission amplitudes)
            [obj.rs1, obj.ts1, ~, ~, obj.tht] = fresnel(obj.nH2O, obj.nSTO, obj.thi); %H2O to STO
            [obj.rs2, obj.ts2]                = fresnel(obj.nSTO, obj.nH2O, obj.tht); %STO to H2O
            
            % Calculate acsoutic scattering prefactors and elasto-optic contribution
            obj.rdeps = obj.cAcoustic*1i*pi./(obj.l.*obj.nSTO.*cos(obj.tht)); % prefactor in eq. 4 in JACS paper
            obj.deps = -obj.p1122.*(obj.nSTO.^4);   % ignore elasto-optic dispersion and calculate dPermxx/dStrain in eq. 4 in JACS paper
            
            % Calculate phase-matched phonon frequency
            obj.freq = omega(obj.cAcoustic, obj.tht, obj.nSTO, obj.l);
      
            % Calculate acoustic pulse shape ignoring optical attenuation
            obj.depsl = pulsef(obj.xi,obj.t0,obj.fwhm,obj.eta0,obj.etaConv,obj.R,obj.cAcoustic,real(obj.freq));
        end
        
        function obj = setp(obj, eta0, xi, t0, fwhm)
% SETP sets acoustic pulse parameters and triggers a recalculation of the
% pulse shape.
%
% obj = obj.setp(eta0, xi, t0, fwhm)
%   Update pulse parameters and recalculate pulse spectrum
%   eta0: (scalar) transformation strain amplitude 0-1
%   xi: (scalar) transformation strain 1/e spatial extent in nm
%   t0: (scalar) transformation strain 1/2 formation time in ps
%   fwhm: (scalar) pump pulse intensity duration fwhm in ps
%
% See Also: SETL, SETLP
           
            obj.eta0 = eta0;
            obj.xi = xi;
            obj.t0 = t0;
            obj.fwhm = fwhm;
            
            % Calculate acoustic pulse shape ignoring optical attenuation
            obj.depsl = pulsef(obj.xi,obj.t0,obj.fwhm,obj.eta0,obj.etaConv,obj.R,obj.cAcoustic,real(obj.freq));
        end

        function obj = setlp(obj, l, eta0, xi, t0, fwhm)
% SETLP sets acoustic pulse parameters and updates the wavelengths which 
% triggers a recalculation of most parameters
%
% obj = obj.setlp(l, eta0, xi, t0, fwhm)
%   Update pulse parameters and recalculate pulse spectrum
%   l: ([nl,1]) wavelengths at which to evaluate spectrum
%   eta0: (scalar) transformation strain amplitude 0-1
%   xi: (scalar) transformation strain 1/e spatial extent in nm
%   t0: (scalar) transformation strain 1/2 formation time in ps
%   fwhm: (scalar) pump pulse intensity duration fwhm in ps
%
% See Also: SETL, SETP

            obj.eta0 = eta0;
            obj.xi = xi;
            obj.t0 = t0;
            obj.fwhm = fwhm;
            obj = obj.setl(l);
        end
    end
    
    % Fast evaluation of phonon model for use with function handles
    methods
        function [Ms, obj] = evalMslp(obj, l, eta0, xi, t0, tau)
% EVALMSLP returns a fast evaluation of the sub CAWs spectrum based on 
% input parameters.
%
% [Ms, obj] = evalMslp(obj, l, eta0, xi, t0, tau)
%   Returns the sub CAWs spectrum [nl, nt] based on parameters:
%   l: ([nl,1]) wavelengths in nm at which to evaluate spectrum
%   eta0: (scalar) transformation strain amplitude 0-1
%   xi: (scalar) transformation strain 1/e spatial extent in nm
%   t0: (scalar) transformation strain 1/2 formation time in ps
%   tau: (scalar) CAWs 1/e decay caused by spectrometer resolution and
%       acoustic damping
%
% See Also: EVALMSP, EVALMLP, EVALMP, lsqnonlin
            
            obj.tau = tau;
            obj = obj.setlp(l,eta0,xi,t0,obj.fwhm);
            Ms = obj.Ms;
        end
        
        function [Ms, obj] = evalMsp(obj, eta0, xi, t0, tau)
% EVALMSP returns a fast evaluation of the sub CAWs spectrum based on input
% parameters.
%
% [Ms, obj] = evalMsp(obj, eta0, xi, t0, tau)
%   Returns the sub CAWs spectrum [nl, nt] based on parameters:
%   eta0: (scalar) transformation strain amplitude 0-1
%   xi: (scalar) transformation strain 1/e spatial extent in nm
%   t0: (scalar) transformation strain 1/2 formation time in ps
%   tau: (scalar) CAWs 1/e decay caused by spectrometer resolution and
%       acoustic damping
%
% See Also: EVALMSLP, EVALMLP, EVALMP, lsqnonlin
            
            obj.tau = tau;
            obj = obj.setp(eta0,xi,t0,obj.fwhm);
            Ms = obj.Ms;
        end
        
        function [Ms, obj] = evalMslt(obj, l, t)
% EVALMSLT returns a fast evaluation of the sub CAWs spectrum based on 
% input parameters.
%
% [Ms, obj] = evalMslt(obj, l, t)
%   Returns the sub CAWs spectrum [nl, nt] based on parameters:
%   l: ([nl,1]) wavelengths in nm at which to evaluate spectrum
%   t: ([nt,1]) wavelengths in nm at which to evaluate spectrum
%
% See Also: EVALMLT, lsqnonlin

           obj.t = t;
           obj = obj.setl(l);
           Ms = obj.M;
        end
        
        function [M, obj] = evalMlp(obj, l, eta0, xi, t0, tau)
% EVALMLP returns a fast evaluation of the CAWs spectrum based on input
% parameters.
%
% [M, obj] = evalMlp(obj, l, eta0, xi, t0, tau)
%   Returns the CAWs spectrum [nl, nt] based on parameters:
%   l: ([nl,1]) wavelengths in nm at which to evaluate spectrum
%   eta0: (scalar) transformation strain amplitude 0-1
%   xi: (scalar) transformation strain 1/e spatial extent in nm
%   t0: (scalar) transformation strain 1/2 formation time in ps
%   tau: (scalar) CAWs 1/e decay caused by spectrometer resolution and
%       acoustic damping
%
% See Also: EVALMSLP, EVALMSP, EVALMP, lsqnonlin
            
            obj.tau = tau;
            obj = obj.setlp(l,eta0,xi,t0,obj.fwhm);
            M = obj.M;
        end
        
        function [M, obj] = evalMp(obj, eta0, xi, t0, tau)
% EVALMP returns a fast evaluation of the CAWs spectrum based on input
% parameters.
%
% [M, obj] = evalMp(obj, eta0, xi, t0, tau)
%   Returns the CAWs spectrum [nl, nt] based on parameters:
%   eta0: (scalar) transformation strain amplitude 0-1
%   xi: (scalar) transformation strain 1/e spatial extent in nm
%   t0: (scalar) transformation strain 1/2 formation time in ps
%   tau: (scalar) CAWs 1/e decay caused by spectrometer resolution and
%       acoustic damping
%
% See Also: EVALMSLP, EVALMSP, EVALMLP, lsqnonlin
            
            obj.tau = tau;
            obj = obj.setp(eta0,xi,t0,obj.fwhm);
            M = obj.M;
        end
        
        function [M, obj] = evalMlt(obj, l, t)
% EVALMLT returns a fast evaluation of the CAWs spectrum based on input
% parameters.
%
% [M, obj] = evalMlt(obj, l, t)
%   Returns the CAWs spectrum [nl, nt] based on parameters:
%   l: ([nl,1]) wavelengths in nm at which to evaluate spectrum
%   t: ([nt,1]) wavelengths in nm at which to evaluate spectrum
%
% See Also: EVALMSLT, lsqnonlin

           obj.t = t;
           obj = obj.setl(l);
           M = obj.M;
        end
        
     end
    
end

function f = omega(v,th,n,l)
% Acoustic scattering frequency as determined by Brillouin Scattering
% phase matching conditions in GHz
%
% f = omega(v,th,n,l)
%   Acoustic frequency l in GHz for acoustic velocity v in nm/ps, incidence
%   angle th in radians, refractive index n, and optical wavelength l in
%   nm. v, th, n, and l can be vectors of the same size and will return f
%   as a vector of the same size.

    f = 2000*v*n.*cos(th)./l;
end

function depsf = pulsef(xi,t0,fwhm,eta0,etaConv,R,c,f)
% PULSEF calculates the frequency domain acoustic spectrum for a pulse
% generated by a transformation strain initiated by a pump pulse. 
%
% depsf = pulsef(xi,t0,fwhm,R,c,f)
%   Returns complex change in strain, deps, generated by a transformation
%   strain that decays in space exponentially with distance xi (nm), grows 
%   exponentially with time t0 (ps), and is initiated with a pulse fwhm
%   intensity duration of fwhm (ps). R is the surface reflection
%   coefficient, c acoustic velocity (nm/ps), and f is a vector of acoustic
%   frequencies (in GHz) to evaluate on.

    spatial_right = @(R,omega) -R./(1+1i*omega); % Fourier transform of eta > 0 for spatial part
    spatial_left = @(omega) 1./(1-1i*omega);     % Fourier transform of eta < 0 for spatial part
    temporal = @(omega) 1./(1+1i*omega);         % Foureir transform of eta > 0 for temporal part
    pump = @(omega) exp(-omega.^2/2);            % Fourier transform of pump temporal intensity profile

    w_spatial = xi*2*pi*f/1000/c;   % Convert spatial frequency to unitless form in Eq. 6
    w_temporal = t0*2*pi*f/1000;    % Convert temporal frequency to unitless form in Eq. 6
    w_pump = fwhm/2.355*2*pi*f/1000; % Convert pump intensity profile frequency to unitless form in Eq. 6

    T_w = xi/c*(spatial_right(R,w_spatial) + spatial_left(w_spatial));  % Fourier transform of spatial pulse shape
    g_w = temporal(w_temporal); % Fourier transform of rate of transformation strain growth
    P_w = pump(w_pump);         % Fourier transform of pump temporal intensity profile
    
    % Complex acousic pulse spectrum in the frequency domain, Eq. 6
    depsf = -eta0*etaConv/2*T_w.*g_w.*P_w; 
end

function y = hvsd(x)
    y = 0.5*(x == 0) + (x > 0);
end