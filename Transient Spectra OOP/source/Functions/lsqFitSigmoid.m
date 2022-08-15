function varargout = lsqFitSigmoid(x, y)
% LSQFITSIGMOID fits vectors x, y to a erf sigmoid function using lsqnonlin
% while determining a good initial guess for the data for fast convergence.
%
% The fit function is:
% y = a+b*0.5*(1+erf(sqrt(0.5)*(x-x0)/s));
% where a is a baseline, b is the amplitude, x0 is the half-rise point, and
% s is the standard-deviation.
%  
% The outputs of this function are the outputs of lsqnonlin:
% [x,resnorm,residual,exitflag,output,lambda,jacobian]
% where x is an array of best-fit values, [a, b, x0, s]
%
% [varargout] = LSQFITSIGMOID(x, y)
%   Fits x, y data to a erf sigmoid as a wrapper to lsqnonlin
%
% See Also: LSQNONLIN
    
    %Define sigmoid fit function
    mySigmoid = @(a,b,x0,s,x) a+b*0.5*(1+erf(sqrt(0.5)*(x-x0)/s));
    
    %Define least squares options
    opts = optimoptions(@lsqnonlin,'TolFun',1e-9,'TypicalX',[0.2 5 0.5 0.2]','Display','off');

    % convert x-y data to column vectors
    y = y(:);
    x = x(:);

    % remove any NaN values 
    TF = ~isnan(y) & ~isnan(x);
    y = y(TF);
    x = x(TF);

    % determine a good initial guess for b and its bounds by determining
    % the sign of the sigmoid
    minVal = min(y);
    maxVal = max(y);

    if maxVal > abs(minVal)
        bG = maxVal;
        x0 = x((y > mean([maxVal, minVal])),1);
        x0 = x0(1);
    else
        bG = minVal;
        x0 = x((y < mean([maxVal, minVal])),1);
        x0 = x0(1);
    end

    % Setup guess for least squares search
    %        a      b    x0      s
    guess = [0,     bG,  x0,     0.2]; %initial guess
    lb =    [-inf, -inf, min(x), 0];
    ub =    [inf, inf,   max(x), inf];
    
    % do the least squares fit
    [varargout{1:nargout}] = lsqnonlin(@(fp) mySigmoid(fp(1),fp(2),fp(3),fp(4),x)-y,guess,lb,ub,opts);
end