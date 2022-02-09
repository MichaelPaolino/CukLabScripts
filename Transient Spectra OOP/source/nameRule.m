classdef nameRule
   properties
       % Struct array that keeps track of rules.
       % The array must be a column with specified fields/field types and the 
       % name field must have unique names
       rules(:,1){mustBeRuleStruct, mustHaveUniqueName}...
             = struct('name', '',...        %char - lookup name
                      'values', {{''}},...  %cell array of char - allowed name values
                      'level', 0,...        %double - nested autoincrement level. Lower number is less significant 
                      'flag', false(0));    %logical - display the rule when generating a name
       
       % Struct array the defines nest levels for autoincrement
       % todo: add type checks
       levels(:,1) = struct('ind',1,...
                            'max',1);
   end
   
   properties (Dependent, Access = private)
       %this is the index of the next free rule position in the rules struct
       ruleInd; 
   end
   
   methods
       function obj = nameRule(varargin)
       % Creates a nameRule object.
       %
       % obj = nameRule()
       %   Default constructor
       %
       % obj = nameRule(rulesCell)
       %   Constructs object and calls obj.add(rulesCell)
       %
       % obj = nameRule(rulesStruct)
       %    Constructs object and calls obj.add(rulesStruct)
       %
       
           switch nargin
               case 0   %default constructor
               case 1   %convert a rules structure or cell to nameRule                  
                   obj = obj.add(varargin{1});
               otherwise
                   error('Too many inputs. Allowed calls: nameRule() or nameRule(ruleStruct)');
           end
       end
       
       function val = get.ruleInd(obj)
       % Return the index of the first empty rule or the last rule + 1
       %
       % val = obj.ruleInd()
       
           if isempty(obj.rules(1).name)
               val = 0;
           else
               val = length(obj.rules); 
           end
       end
       
       function names = list(obj)
       % Return a cell array of all rule names
       %
       % names = list(obj) 
       
          names = {obj.rules.name}; 
       end
       
       function ind = search(obj,ruleName)
       % Returns the index of rule with name ruleName. Throws an error if
       % the rule was not found
       %
       % ind = search(obj,ruleName) 
       
           names = list(obj);
           ind = find(strcmp(names,ruleName));
           assert(~isempty(ind),...
               ['Rule ' ruleName ' not found. Available rules are ' strjoin(names,', ') '.']);
       end
       
       function obj = add(obj,varargin)
       % Adds rules to the end of the nameRule object rules structure. 
       %
       % obj = add(obj,ruleCell)
       %   Rule cell is a cell array {nRules, 4} that has rules as rows and 'name'
       %   'values', 'level', and 'flag' entries as columns in that order
       %
       % obj = add(obj,ruleStruct)
       %   Adds a valid ruleStruct (nRules,1) with fields 'name','values','level',
       %   and 'flag'
       %
       % obj = add(obj,name,values,level,flag)
       %   Adds a single rule with char name, cell values, double int and logical
       %   flag
       
           switch nargin
               case 2   %add cell array or struct array of inputs
                   if iscell(varargin{1})   %for cell array input
                       %Convert cell array to struct. Cell array columns should 
                       %correspond to different fields and rows to rules
                       ruleStruct = struct('name', varargin{1}(:,1),...
                                           'values', varargin{1}(:,2),...  
                                           'level', varargin{1}(:,3),...
                                           'flag', varargin{1}(:,4));
                       obj = obj.add(ruleStruct);   %call add on struct
                   elseif isstruct(varargin{1})
                       nRules = length(varargin{1});    %number of rules in struct
                       obj.rules(obj.ruleInd + (1:nRules)) = varargin{1};   %add rules to the end of obj.rules
                   else
                       error(['Unsupported input type for add. Expected cell or struct. Got ' class(varargin{1}) '.']);
                   end    
               case 5
                   %directly add a single rule to obj.rules
                   obj.rules(obj.ruleInd+1) = struct('name', varargin{1},...
                                             'values', {varargin{2}},...  
                                             'level', varargin{3},...
                                             'flag', varargin{4});
               otherwise
                   error('Invalid number of inputs. Expected add(ruleStruct) or add(ruleName, ruleValues, rulelevel, ruleFlag).');
           end
       end
       
       function obj = modify(obj,varargin)
       % Modifies a rule field specified by name-value pairs 
       % The rule and field is selected by the name in the name-value pair
       % as 'ruleName.field'. If a rule's name is changed, subsequent
       % name-value pairs should refer to the old name.
       %
       % obj = modify(obj,'nameA.fieldName1',newValue1,...
       %                  'nameB.fieldName2',newValue2....)
       %   Modifies rule nameA's fieldName1, rule nameB's fieldName2, etc. 
       %   to values newValue1, newValue2, etc. 
       %
       % Allowed field names are:
       %   'name', 'values', 'level', and 'flag'. Corresponding value types are
       %    char,   cell,     double, and  logical.
       %
           
           % loop through varargin to update rule fields
           oldObj = obj;    %make a copy of the old object to confusion if a rule's name is updated
           for ii = 1:2:(nargin-2)
               dotInd = strfind(varargin{ii},'.');  %find the . separator
               
               ruleName = varargin{ii}(1:dotInd-1); %get the ruleName
               fieldName = varargin{ii}(dotInd+1:end); %get the fieldName
               
               ind = oldObj.search(ruleName);  %get the name index from the old object
               assert(any(strcmp(fieldName,fields(obj.rules))),[fieldName ' is not a valid field name. Valid fields are ' strjoin(fields(obj.rules),', ') '.']);
               obj.rules(ind).(fieldName) = varargin{ii+1}; %update the rule    
           end
       end
       
       function obj = remove(obj,ruleName)
       % Removes a from the nameRule object
       %
       % obj = remove(obj,ruleName)
       %    Removes the rule named ruleName from the nameRule object.
       
           ind = obj.search(ruleName);
           obj.rules(ind) = [];
       end
       
       function obj = rearrange(obj, varargin)
       % Rearranges the rule order in the rules structure
       %
       % obj = rearrange(obj, newIndOrder)
       %    returns an object with rules ordered by the double array newIndOrder. 
       %    Equivalent to calling obj.rules = obj.rules(newIndOrder)
       %
       % obj = rearrange(obj, 'name1', 'name2',...)
       %    returns an object with rules orderd by char 'name1', 'name2', etc.
       
           if isa(varargin{1},'double') %update rule order by index
               obj.rules = obj.rules(varargin{1});
           elseif isa(varargin{1},'char') %update rule order by name order
               oldObj = obj;
               obj = nameRule();
               for ii = 1:nargin-1  %loop through each argument
                  ind = oldObj.search(varargin{ii});    %find name index in old rule
                  obj = obj.add(oldObj.rules(ind));     %add old rule to new rule object
               end
           else
               error('Unexpected input class for obj.rearrange');
           end
       end
       
       function obj = buildLevels(obj)
       % Builds autoincrement level struct from entries in the rules struct.
       % The level number corresponds to the struct array row index. The current 
       % level index, level.ind is automatically set to 1. For a given level, the
       % maximum index is set to the length of the shortest rules.value array at
       % that level. Rules that have their flags set to false are ignored in the
       % max calculation. If a level's max index cannot be calculated, it is set
       % to one, which effectively skips the level when running generate.
       %
       % obj = buildLevels(obj)
       %   returns a nameRule object with updated level struct array.
           
           %get available levels and flags in rules
           levelVals = [obj.rules.level];
           flagVals = [obj.rules.flag];
           
           %remove duplicate levels and ignore rules that have flag = false
           uniqueLevels = unique(levelVals.*flagVals);
           %remove level 0, which is a flag for no autoincrement
           uniqueLevels = uniqueLevels(uniqueLevels>0);
           
           %get number of elements of value cell arrays in rules
           valueSizes = zeros(1,obj.ruleInd);
           for ii = 1:obj.ruleInd
               valueSizes(ii) = length(obj.rules(ii).values);
           end
           
           %reset values for levels struct array
           obj = obj.setLevelMax(ones(1,max(uniqueLevels)));
           obj = obj.reset;
           
           %populate new levelVals
           for ii = uniqueLevels
               %take the smallest length to avoid out-of-index conditions
               obj.levels(ii).max = min(valueSizes((levelVals.*flagVals)==ii)); 
           end
       end
       
       function obj = reset(obj)
       % Resets all autoincrement level indecies back to 1
       %
       % obj = reset(obj)
       %    Returns an object with all elements of obj.levels.ind = 1
           
           obj = obj.setLevelIndex(ones(1,length(obj.levels)));
       end
       
       function obj = setLevelIndex(obj,indList)
       % Sets the index for all levels while preserving old max values.
       % The length of levels is automatically resized to the length of
       % indList.
       %
       % obj = setLevelIndex(obj,indList)
       %   set structure array elements of obj.level.ind to elements of 
       %   double array indList.
       
          nElem = length(indList);
          newMax = ones(1,nElem);
          newMax(1:length(obj.levels)) = [obj.levels.max];
          newMax = newMax(1:nElem);
            
          obj.levels = struct('ind',num2cell(indList),...
                              'max',num2cell(newMax));
       end
       
       function obj = setLevelMax(obj,maxList)
       % Sets the max index for all levels while preserving old ind values.
       % The length of levels is automatically resized to the length of
       % maxList.
       %
       % obj = setLevelMax(obj,maxList)
       %   set structure array elements of obj.level.max to elements of 
       %   double array maxList.
          
          nElem = length(maxList);
          newInd = ones(1,nElem);
          newInd(1:length(obj.levels)) = [obj.levels.ind];
          newInd = newInd(1:nElem);
            
          obj.levels = struct('ind',num2cell(newInd),...
                              'max',num2cell(maxList));
       end
       
       function obj = increment(obj)
       % Increments the most significant (last) index in obj.levels by one.
       % If any index exceeds its max value, its value is reset to 1 and the less
       % significant index is incremented instead. If all indicies have been
       % incremented to max, the values for all indicies are reset to 1.
       %
       % obj = increment(obj)
       %   Returns an object with incremented indecies.
          
          %check to make sure there are levels to increment
          if ~isempty(obj.levels)
              %start by incrementing the most significant index
              obj.levels(end).ind = obj.levels(end).ind+1;

              %check, starting from the most significant index, that all
              %indicies are at or below their max value
              for ii = length(obj.levels):-1:1
                  if obj.levels(ii).ind > obj.levels(ii).max
                      if ii==1  %required to avoid subindexing 0
                          obj = obj.reset;
                      else %reset the current index value to 1 and increment the less significant index
                          obj.levels(ii).ind = 1;
                          obj.levels(ii-1).ind = obj.levels(ii-1).ind+1;
                      end
                  end
              end
          end
       end
       
       function [nameOut, obj] = buildName(obj,varargin)
       % Generates a name based on the nameRule rules and autoincrement levels.
       % This method works by looping over and executing each rule and
       % concatonating the result to the nameOut char array. The method supports
       % autoincrementing the level indicies for use with nested for loops to
       % generate incremental sequences of names. By default, autoincrement is not
       % enabled, but can be enabled with the 'autoIncrement', true name-value 
       % pair. In addition, the method skips rules whose flag is set to false. The
       % flags can be set by executing obj = obj.modify(ruleName,'flag',condition)
       % before running this method. 
       % 
       % [nameOut, obj] = buildName(obj)
       %   Executes the name rules without autoincrement and returns char array
       %   nameOut. The returned obj is a copy if the input obj.
       %
       % [nameOut, obj] = buildName(obj,'name1',value1,'name2',value2,...)
       %   Executes the name rules with the options set the name-value pairs. This
       %   call may return a modified object depending on the options. 
       %
       % name-value pairs:
       %   'delimiter': adds a delimiter between every rule with a true flag. The
       %       delimiter is specified by a char array value.
       %   'autoIncrement': a logical true value turns autoincrement on and 
       %       returns an object with incremented level indicies. When a level 
       %       index reaches levels.max, it resets to 1 and increments the less 
       %       significant index (lower row or level number)
           
           %Set defaults
           autoInc = false;
           delimiter = '';
           
           %update name-value pairs
           for ii = 1:2:(nargin-1)
               assert(ischar(varargin{ii}),'Name argument must be of type char.');
               switch varargin{ii}
                   case 'delimiter'
                       delimiter = varargin{ii+1};
                   case 'autoIncrement'
                       autoInc = varargin{ii+1};
                   otherwise
                       error([varargin{ii} ' is not a valid name-value pair']);
               end
           end
           
           %verify that each rule has a level defined in the level struct
           assert(max([obj.rules.level].*[obj.rules.flag])<=length(obj.levels),...
               ['The rules struct references more levels than are defined in obj.levels.',...
                'Run obj.buildLevels() before obj.execute and any increment loops to update the levels struct.']);
           
           %verify that each level index is less than or equal to the max
           assert(all([obj.levels.ind]<=[obj.levels.max]),'A levels struct ind is set to a value higher than max.');
            
           %loop through and execute each rule
           nameOut = '';
           for ii = 1:obj.ruleInd
               if obj.rules(ii).flag
                   %update the index of which value the rule uses
                   valueInd = 1;
                   levelInd = 0;
                   if obj.rules(ii).level ~= 0
                       levelInd = obj.rules(ii).level;
                       valueInd = obj.levels(levelInd).ind;
                   end
                   
                   %assert that the value index in the current rule corresponds to a value in the cell array
                   assert(valueInd <= length(obj.rules(ii).values),...
                       ['The level ' num2str(levelInd) ' index = ' num2str(valueInd),...
                        ' is larger than the available values in rule ' obj.rules(ii).name,...
                        '. The length of values are ' num2str(length(obj.rules(ii).values)) '.']);
                   
                   %concatonate name
                   nameOut = [nameOut obj.rules(ii).values{valueInd} delimiter];
               end
           end
           
           %remove last delimiter
           nameOut = nameOut(1:end-length(delimiter));
           
           %autoincrement if flagged
           if autoInc
               obj = obj.increment;
           end
       end
   end
end

% Validation functions
function mustBeRuleStruct(a)
%validate subtypes of rule struct fields

   assert(all(strcmp(fields(a),{'name';'values';'level';'flag'})),...
        ['Input structure must have field names: name, values, level, flag. Got: ' strjoin(fields(a),', ') '.']);
   assert(isa(a(1).name,'char'),...
        ['Name field must be a char array. Got ' class(a(1).name) '.']);
   assert(isa(a(1).values,'cell'),...
        ['Values field must be a cell array of char arrays. Got ' class(a(1).values) '.']);
   assert(isa(a(1).values{1},'char'),...
        ['Values field must be a cell array of char arrays. Got ' class(a(1).values) '.']);
   assert(isa(a(1).level,'double'),...
        ['level field must be of type double. Got ' class(a(1).level) '.']);
   assert(isa(a(1).flag,'logical'),...
        ['Flag field must be of type logical. Got ' class(a(1).flag) '.']);
end

function mustHaveUniqueName(a)
%the name field in rule struct must have unique names

    assert(length(unique({a.name}))==length(a),...
        ['Rule names must be unique. Rule names were ' strjoin({a.name},', ') '.']);
end