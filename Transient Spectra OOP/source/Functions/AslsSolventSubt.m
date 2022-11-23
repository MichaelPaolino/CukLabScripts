function x = AslsSolventSubt(varargin)
% ASLSSOLVENTSUBT estimates a smooth subtractor for solvent subtraction 
% using asymmetric least squares (ASLS) method. The goal is to find
% x (subtractor) in y = kx, where y is sample spectrum and k is the solvent
% spectrum such that a baseline is generated and negative peaks are
% avoided.
%
% The problem can be stated as a weighted least squares minimization with
% smoothing penalites on the final spectrum and subtractor:
% L=argmin((y-kx)'W(y-kx)+l1*Abs(D(x))+l2*Abs(D(y-kx)),x)
% Here, W is a diagonal wieght matrix, l1 and l2 are smoothing penalties,
% and D is the second order difference operator. 
%
% Diagonal elements of the weight matrix are set to 1 when the subtracted 
% data point is positive, i.e. if y-kx>0, and to a large number h if the
% subtracted point is negative , i.e. if y-kx<0. The algorithm runs
% iteratively starting from W = I, the identity matrix, and updates W in
% each iteration until the changes in x converge or the maximum number of
% iterations are reached.
% 
% Typical values for h and the penalites are:
%   h = 1e6, penalty on negative points, h<1e9 (numericaly usntable above)
%   l1 = 1e7, penalty on smoothness of the subtraction vector x and l1>=h
%   l2 = 1e5, penalty on smoothness of the  final spectra y-k*x and l2<h
%
% x = ASLSSOLVENTSUBT(y, k)
%   Calculates the smooth subtractor x for sample spectrum y and background
%   spectrum k using the default penalties above. The algorithm stops at 10
%   iterations and if the change in x between iterations is < 0.5
%
% x = ASLSSOLVENTSUBT(y, k, 'name', value)
%   The same call as above with additional name-value parameters. Allowed
%   names are:
%       'params': 1x3 vector for custom penalties [l1, l2, h]
%       'maxRec': int for the maximum number of recursions. Default is 10
%       'convCrit': scalar for the convergance criteria. Default is 0.5
%       'silent': logical on whether to display information for each
%           iteration. The default is true (do not display).
%
% Written 11/12/2014 Ilya Vinogradov
% Modified 11/13/2014 Updated to exact solution to minimization problem
% Modified 11/21/2014 Update to include penalty for negative values of x
% Modified 9/21/2022 Updated to use inputparser for optional inputs
% See paper for more: Yang, C.; Peng, S; Xie, Q.; Peng, J; Wei, J.; Hu, Y., 
%   FTIR Spectral Subtraction Based on Asymmetric Least Squares.
%   International Conference on BMEI. 2011, 4, 1724.
    
    % Parse optional inputs
    p = inputParser();
    p.addRequired('y');
    p.addRequired('k');
    p.addParameter('params',[1e7,1e5,1e6]);
    p.addParameter('maxRec',10);
    p.addParameter('convCrit',0.05);
    p.addParameter('silent',true);
    
    p.parse(varargin{:});

    y = p.Results.y(:);    % sample spectra
    k = p.Results.k(:);    % solvent spectra
    
    l1 = p.Results.params(1);   %smoothnes of subtraction vector
    l2 = p.Results.params(2);   %smoothnes of subtracted spectra
    h1 = p.Results.params(3);   %penalty on negative values of subtracted spectra
    h2 = 0;   %penalty on negative values of subtraction vector
    
    maxRec = p.Results.maxRec;
    convCrit = p.Results.convCrit;
    talk = ~p.Results.silent;
    
    %quadratic difference matrix
    matsize=size(y,1);
    d=diff(speye(matsize),2)'*diff(speye(matsize),2);

    %will use solvent spectra as sparse diagonal matrix (makes math easier)
    k=sparse(diag(k));
    W=speye(matsize);  %initial guess for W is a sparse diagonal matrix
    V=speye(matsize);  %initial guess for V

    x0=ones(matsize,1);    %initial guess -> linear subtraction
    x=x0;                  %use initial guess in optimization
    
    for j=1:maxRec
        %we want solution to:
        %L=argmin((y-kx)'W(y-kx)+l1*Abs(Dx) + l2*Abs(D(y-kx)),x)
        W=sparse(diag(1+(y<spdiags(k).*x)*h1));  %Negative peak penalty
        V=sparse(diag((x<0)*h2));   %Negative values of x penalty
        H=2*(W*k^2+l1*d+l2*d*k^2+V);  %calculate quadratic term with smoothing
        H=(H+H')/2; %Ensures hessian is symmetric (numerical error?)
        f=-2*(W*k*y+d*k*y);     %calculate linear term with negative penalty
        x=-inv(H)*f;    %exact solution when H is pos-def
        L=0.5*x'*H*x+f'*x+y'*W*y+y'*d*y;  %Value of minimized answer
        convergence=sum((x-x0).^2); %check for convergence using least squares
        if convergence<=convCrit
            if convergence == 0 && talk
                disp(['Problem converged. Exact solution found after ',...
                    num2str(j),' passes.']);
            elseif talk
                disp(['Convergence criteria met after ', num2str(j),...
                    ' passes. Convergence = ', num2str(convergence)]);
            end
            if talk
                disp(['Minimized solution L = ', num2str(L)]);
            end
            break;
        end
        x0=x;   %use last solution as next guess
        if talk
            disp(['Pass ', num2str(j), ' complete! Convergence = ',...
                num2str(convergence)]);
            disp(['Minimized solution L = ', num2str(L)]);
        end
    end

    if convergence>convCrit
        warning(['Warning: Convergence criteria not met after %d passes.\n',...
                 'Convergence = %g > %g\n',...
                 'Minimized solution L = %g'], maxRec, convergence, convCrit, L)
    end
end