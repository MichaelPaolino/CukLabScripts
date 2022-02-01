function xVals = findGridIntercepts(x, y, yTarget)
%Finds all intercept at value yTarget in x, y data. Intercepts are linearly
%interpolated on x. Note: this function only finds intercepts that cross
%yTarget. Constant-values at yTarget and local extrema may not be found.
%Inputs:
%   x: 1d array of grid values that define y, [nx, 1]
%   y: 1d array of lookup values of x. Must be the same length as y [nx, 1] 
%   yTarget:  scalar that specifies intercept value to search for (e.g.
%       x-intercept is yTarget = 0)
%Output:
%   xVals: interpolated x values for the intercept. If no intercepts were
%       found, xVals is empty. 
%

    %pre-formatting of inputs
    x = x(:);
    y = y(:);
    nx = length(x);
    assert(nx==length(y), ['Length of x must match length of y.',...
                           '\nlength(x) = ' num2str(nx) '\nlength(y) = ' num2str(length(y))]); 
    
    %find all intercepts using threshold difference of adjacent grid points.
    %Any difference that is not zero is an intercept on the grid
    int = diff(y>yTarget)~=0;
    intInd1 = find(int);
    
    %ensure that intercept points were found
    if ~isempty(intInd1)
        %determine the adjacent point around the intercept
        direc = int(intInd1);
        intInd2 = intInd1+direc;

        %take care of out-of-index points
        highClip = intInd2 > nx;
        lowClip = intInd2 < 1;
        intInd2(highClip) = intInd2(highClip) - 2;
        intInd2(lowClip) = intInd2(lowClip) + 2;

        %construct x value pairs and y value pairs for intercept interpolation
        xPairs = [x(intInd1), x(intInd2)];
        yPairs = [y(intInd1), y(intInd2)];

        %loop through each intercept and fit a line through the points
        %adjacent to the intercept
        nInt = length(intInd1);
        linInterp = zeros(nInt,2);
        for ii = 1:nInt
            linInterp(ii,:) = polyfit(xPairs(ii,:),yPairs(ii,:),1);
        end

        xVals = (yTarget - linInterp(:,2))./linInterp(:,1);
    else  %no intercept points found
        xVals = [];
    end
    
end