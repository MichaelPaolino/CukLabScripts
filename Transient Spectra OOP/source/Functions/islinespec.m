function isls = islinespec(x)
% ISLINESPEC checks if input x is a linespec input for the plot function.
% https://stackoverflow.com/questions/8486904/passing-varargs-to-plot-function
%
% isls = ISLINESPEC(x)
%   Returns logical true if char array x is a valid linespec. Returns false
%   otherwise.
%
% See also PLOT

    isls = false;

    if ~ischar(x)
        return;
    end

    lineStyleSpecifiers = {'--','-.','-',':'};
    markerSpecifiers    = {'square','diamond','pentagram','hexagram','+','o','*','.','x','s','d','^','v','>','<','p','h'};
    colorSpecifiers     = {'r','g','b','c','m','y','k','w'};

    for oo=1:length(lineStyleSpecifiers)
        k = strfind(x,lineStyleSpecifiers{oo});
        if ~isempty(k)
            x(k:(k+length(lineStyleSpecifiers{oo})-1)) = [];
            break;
        end
    end

    for oo=1:length(markerSpecifiers)
        k = strfind(x,markerSpecifiers{oo});
        if ~isempty(k)
            x(k:(k+length(markerSpecifiers{oo})-1)) = [];
            break;
        end
    end

    for oo=1:length(colorSpecifiers)
        k = strfind(x,colorSpecifiers{oo});
        if ~isempty(k)
            x(k:(k+length(colorSpecifiers{oo})-1)) = [];
            break;
        end
    end

    if isempty(x)
        isls = true;
    end
end