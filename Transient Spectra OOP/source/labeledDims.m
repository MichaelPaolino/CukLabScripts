classdef labeledDims < double
    properties
        dims = struct('name',{''},...
                      'size',[]);
    end
    
    methods
        function obj = labeledDims(varargin)
            %handles adding double data by superclass constructor           %constr
            if nargin == 0  %default constructor, empty array
                array = [];
            else %pass first argument as double array to double superclass
                array = varargin{1};
            end
            obj = obj@double(array);    %call double superclass constructor
            
            %handle additional inputs for units definitions
            switch nargin
                case 0  %default constructor, handled above, do nothing
                case 1  %double constructor, handled above, do nothing
                otherwise
                    error(['labeledDims expects up to one input. Got ' num2str(nargin) '.']);
            end
            
            
        end
    end
    
   %Override implementations from double class
    methods
        %indexed subset of array, e.g. myArray(1:10)
        function sref = subsref(obj,s)
            switch s(1).type  
                case '.' %dot indexing, e.g. myArray.property
                    switch s(1).subs
                        case 'data' %myArray.data returns type double
                            d = double(obj);
                            if length(s)<2 %myArray.data
                                sref = d;
                            elseif length(s)>1 && strcmp(s(2).type,'()') %myArray.data(ind)
                                sref = subsref(d,s(2:end));
                            end
                        otherwise
                            sref = subsref@double(obj,s);   %perform standard dot-notation call from the double superclass
                   end
                case '()' %array indexing, e.g. myArray(ind) returns type doubleWithUnits
                    d = double(obj);
                    newd = subsref(d,s(1));
                    obj = doubleWithUnits(newd,obj);
                    if length(s)<2 %myArray(ind)
                        sref = obj;
                    elseif length(s)>1 && strcmp(s(2).type,'.') %myArray(ind).unit
                        sref = subsref(obj,s(2:end));
                    else
                        error('Not a supported indexing expression');
                    end
                case '{}'
                    error('Not a supported indexing expression');
            end
        end
        
        %assign subset of an array, e.g. myArray(1:10) = 1:10;
        function obj = subsasgn(obj,s,b)
            switch s(1).type
                case '.'
                    switch s(1).subs
                        case 'data'
                            if length(s)<2
                                obj = doubleWithUnits(b,obj);
                            elseif length(s)>1 && strcmp(s(2).type,'()')
                                d = double(obj);
                                newd = subsasgn(d,s(2:end),b);
                                obj = doubleWithUnits(newd,obj);
                            end
                        otherwise
                            obj = subsasgn@double(obj,s,b);   %perform standard dot-notation call from the double superclass
                    end
                case '()'
                   d = double(obj);
                   newd = subsasgn(d,s(1),b);
                   obj = doubleWithUnits(newd,obj);
                case '{}'
                    error('Not a supported indexing expression')
            end
        end
        
        %horizontal array concatonation, e.g. [myArray1, myArray2]
        function obj = horzcat(varargin)
            %format varargin to be ready for array concatonation
            [arrayCell, obj, oldUnitInd] = formatArrayCat(varargin{:});
            
            %superclass horizontal concatonation
            array = horzcat(arrayCell{:});
            
            %rebuild concatonated array object and convert to old current
            %unit of the first element
            obj = doubleWithUnits(array,obj);
            obj = convert(obj,oldUnitInd);
        end
        
        %vertical array concatonation, e.g. [myArray1; myArray2]
        function obj = vertcat(varargin)
            %format varargin to be ready for array concatonation
            [arrayCell, obj, oldUnitInd] = formatArrayCat(varargin{:});
            
            %double class vertical concatonation on double cell array
            array = vertcat(arrayCell{:});
            
            %rebuild concatonated array object and convert to old current
            %unit of the first element
            obj = doubleWithUnits(array,obj);
            obj = convert(obj,oldUnitInd);
        end
        
        %formats doubleWithUnits argument input for array concatonation methods
        function [arrayCell, obj, oldUnitInd] = formatArrayCat(varargin)
            %todo: check that objects have the same base unit to avoid unexpected results
            
            %remember current unit index of the first doubleWithUnits before conversion to base unit 
            oldUnitInd = varargin{1}.currentUnitInd;
            
            %convert doubleWithUnits cell array to base unit
            arrayCell = cellfun(@(fObj) convert(fObj,fObj.baseUnitInd),varargin,'UniformOutput',false);
            obj = arrayCell{1}; %keep the first doubleWithUnits properties
            
            %convert doubleWithUnits cell array to double cell array
            arrayCell = cellfun(@double,arrayCell,'UniformOutput',false);
        end
    end
end