function varargout = lsqFitExpGrow(x, y)
% LSQFITEXPGROWD fits vectors x, y to an exponential growth using lsqnonlin
% while determining a good initial guess for the data for fast convergence.
%
% The fit function is:
% ht = heaviside(t-t0)
% y = (a+b*ht-c*(1-exp(-ht*t/tau)))**g(t,fwhm)
% where t0 is time zero, a is a baseline, b is the amplitude of the 
% absorptive component, c is the amplitude of the emissive component, and
% tau is the emissive time constant
%  
% The outputs of this function are the outputs of lsqnonlin:
% [x,resnorm,residual,exitflag,output,lambda,jacobian]
% where x is an array of best-fit values, [a, b, c, tau, fwhm, t0]
%
% [varargout] = lsqFitExpGrow(x, y)
%   Fits x, y data to an exp growth as a wrapper to lsqnonlin
%
% See Also: LSQNONLIN
    
    %Define exp growth fit function
    f = 1/sqrt(2*pi)*exp(-(5:0.5:5).^2/2);
    myExpGrow = @(a,b,c,tau,fwhm,t0,t) nonuniformConv(a+hvsd(t-t0).*(b-c*(1-exp(-hvsd(t-t0).*(t-t0)/tau))),t,2.355*f/fwhm,0.5*fwhm/2.355);
    
    % convert x-y data to column vectors
    y = y(:);
    x = x(:);

    % remove any NaN values 
    TF = ~isnan(y) & ~isnan(x);
    y = y(TF);
    x = x(TF);

    % determine a good initial guess for b and its bounds by determining
    % the sign of the sigmoid
    bg = max(y);
    cg = bg-min(y);
    pp = find(y > (bg/2));
    tg = x(pp(1));
    
    % Setup guess for least squares search
    %           a       b     c     tau    fwhm   t0
    guess =    [0,      bg,   cg,   0.4,   0.3,   tg]; 
    typicalX = [0.1,    1,    1,    1,     0.3,   0.5];
    lb =       [-inf,   0,    0,    0.1,   0.1,  min(x)];
    ub =       [inf,    inf,  inf,  10,    5,    max(x)];
    
    %Define least squares options
    opts = optimoptions(@lsqnonlin,'TolFun',1e-9,'TypicalX',typicalX,'Display','off');
    
    % do the least squares fit
    [varargout{1:nargout}] = lsqnonlin(@(fp) myExpGrow(fp(1),fp(2),fp(3),fp(4),fp(5),fp(6),x)-y,guess,lb,ub,opts);
end

function y = hvsd(x)
    y = 0.5*(x == 0) + (x > 0);
end