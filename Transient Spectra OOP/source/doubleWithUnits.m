classdef doubleWithUnits < double
    properties
        %array data is stored in double superclass

        unitRules = struct('unit','none',...
                           'displayName','none',...
                           'base2unit', @(f) f,...
                           'unit2base', @(f) f);
    end
    
    properties (Access = private)
                %todo: implement type check
        baseUnitInd = 1;   %index of base unit rule
        currentUnitInd = 1;    %index of current unit rule
    end
    
    properties (Dependent)
        unit
        base
        dispName
    end
    
    %doubleWithUnits constructor and get/set methods
    methods
        %constructor
        function obj = doubleWithUnits(varargin)
            %handles adding unitless double array data by caling double superclass constructor           %constr
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
                case 2  %add new array but copy all other object properties
                    copyObj = varargin{2};
                    obj.baseUnitInd = copyObj.baseUnitInd;
                    obj.currentUnitInd = copyObj.currentUnitInd;
                    obj.unitRules = copyObj.unitRules;
                case 3  %add array, assign base unit and current unit, assign short name, no additional rules
                    obj = obj.changeBaseName(varargin{2},  varargin{3});
                otherwise
                    error(['doubleWithUnits expects up to three inputs. Got ' num2str(nargin) '.']);
            end
            
            
        end
        
        %returns name of the current unit
        function val = get.unit(obj)
        	val = obj.unitRules(obj.currentUnitInd).unit; 
        end
        
        %returns base unit
        function val = get.base(obj)
        	val = obj.unitRules(obj.baseUnitInd).unit; 
        end
        
        %returns display name of current unit
        function val = get.dispName(obj)
        	val = obj.unitRules(obj.currentUnitInd).displayName; 
        end
        
        %sets the current unit by calling convert
        function obj = set.unit(obj,targetUnit)
            obj = convert(obj, targetUnit); 
        end
        
        %changes the base unit. Did not use a set method since both the unit
        %and display name need to be changed and set requires exactly 2 inputs
        function obj = changeBaseName(obj,newUnit,newDisplayName)
            obj.unitRules(obj.baseUnitInd).unit = newUnit;
            obj.unitRules(obj.baseUnitInd).displayName = newDisplayName;
        end
        
    end
    
    %doubleWithUnits methods
    methods
        %add a unit and its conversion rules
        function obj = addRule(obj,unitName,displayName,base2unit,unit2base)
            rule = struct('unit',unitName,...
                          'displayName',displayName,...
                          'base2unit', base2unit,...
                          'unit2base', unit2base);
            obj.unitRules = [obj.unitRules rule];
        end

        %update an existing unit and its conversion rules
        function obj = updateRule(obj,unitName,displayName,base2unit,unit2base)
            %get index of unit that will be updated
            ruleInd = obj.getConversion(unitName);
            
            %if array is in units that are being updated, convert back to
            %the base unit before updating rule
            isUpdateUnit = false;
            if ruleInd == obj.currentUnitInd
                obj = convert(obj, obj.baseUnitInd);
                isUpdateUnit = true;
            end
            
            %update the unit rule
            obj.unitRules(ruleInd) = struct('unit',unitName,...
                          'displayName',displayName,...
                          'base2unit', base2unit,...
                          'unit2base', unit2base);
            
            %if array is in units that are being updated, convert from
            %base unit back to unit that was updated
            if isUpdateUnit
                obj = convert(obj, ruleInd);
            end

        end
        
        %converts the array to the target unit
        function obj = convert(obj, targetUnit)          
            %get the conversion rule index for the targetUnit if it is a char
            if ischar(targetUnit)
                targetUnitInd = obj.getConversion(targetUnit);
            else
                targetUnitInd = targetUnit;
            end
            
            %get array from double superclass
            array = double(obj);
            
            %if the array is not in the base unit, convert it back to base unit first
            if obj.baseUnitInd ~= obj.currentUnitInd
                array = obj.unitRules(obj.currentUnitInd).unit2base(array);
            end
            
            %convert to to target unit. Conversion fails for empty arrays
            if ~isempty(array)
                array = obj.unitRules(targetUnitInd).base2unit(array);
            end
            obj.currentUnitInd = targetUnitInd;
            
            %reconstruct object with new array data
            obj = doubleWithUnits(array,obj);
        end
        
        %list available unit rules
        function unitsOut = listUnits(obj)
            nUnits = length(obj.unitRules);
            unitsOut = cell(1,nUnits);
            for ii = 1:nUnits
                unitsOut(ii) = {obj.unitRules(ii).unit};
            end
        end
        
        %return the unit rule index for targetUnit
        function unitInd = getConversion(obj,targetUnit)
        %Finds the index for a specific unit
        %Inputs:
        %   targetUnit: character array for the unit to find, i.e. 'nm'
        %Outputs:
        %   unitInd: index of obj.unitRules that corresponds to unitStruct
            units = listUnits(obj);
            unitInd = find(strcmp(units,targetUnit));
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

