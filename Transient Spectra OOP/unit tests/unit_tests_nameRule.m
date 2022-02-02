% Unit tests for nameRule class

%% test constructors
myRule = nameRule();   %default constructor

%constructor 
nameStruct = struct('name', 'asdf',...
                     'values', {{'1','2','3'}},...  
                     'level', 1,...
                     'flag', false());
                 
myRule = nameRule(nameStruct);
assert(length(myRule.rules) == 1, 'Failed to add or update struc with nameRule constructor');

%% test adding rules
myRule = nameRule();   %default constructor

%add using cell array. This also tests adding using struct because of the recursive call
myRule = myRule.add({'asdf',{'1','2','3'},1,true;...
                   'ghjk',{'4'},5,false});
assert(all(strcmp({myRule.rules.name},{'asdf','ghjk'})),'addRule failed to add name rules using cell input');

%add using argument input
myRule = myRule.add('qwer',{'1','2','3'},1,true);
assert(all(strcmp({myRule.rules.name},{'asdf','ghjk','qwer'})),'addRule failed to add name rules using argument input');

try
    myRule = myRule.add('qwer',{'1','2','3'},1,true);
    error('Failed to stop adding non-unique rule name');
catch ME
    assert(contains(ME.message,'Rule names must be unique'),'Failed to stop adding non-unique rule name');
end

%% test modifying rules
myRule = testRule1();    %build myRule object with test data

%change a single entry
myRule = myRule.modify('asdf.flag',false);
assert(~myRule.rules(1).flag,'modifyRule failed to update flag');

%change two entries
myRule = myRule.modify('ghjk.name','qwer','ghjk.flag',true);
assert(strcmp(myRule.rules(2).name,'qwer') && myRule.rules(2).flag,'modifyRule failed to update name and flag');

%% test removing a rule
myRule = testRule1();    %build myRule object with test data
               
myRule = myRule.remove('asdf');
assert(length(myRule.rules) == 1,'Failed to remove rule');

%% test rearranging rules
myRule = testRule1();    %build myRule object with test data

%rearrange rule by index
myRule = myRule.rearrange([2,1]);
assert(all(strcmp({myRule.rules.name},{'ghjk','asdf'})),'Failed to rearrange rules by index');

%rearrange rule by name
myRule = myRule.rearrange('asdf','ghjk');
assert(all(strcmp({myRule.rules.name},{'asdf','ghjk'})),'Failed to rearrange rules by name');

%% test building autoincrement levels
%basic level test
[myRule,wl,rpt,schemes,gPos,desc] = testRule2();
myRule = myRule.buildLevels;
assert(all([myRule.levels.max]==[length(schemes),length(gPos),length(rpt),length(wl)]),...
    'Failed to build levels')

%test building levels when there are cell arrays of different length at each level
myRule = myRule.add('test',[schemes(:); schemes(:)],1,true);
myRule = myRule.rearrange('desc','scheme','test','gpos','rpts','wl');
myRule = myRule.buildLevels;
assert(all([myRule.levels.max]==[length(schemes),length(gPos),length(rpt),length(wl)]),...
    'Failed to build levels');

%test building levels when rules have a false flag
myRule = myRule.modify('rpts.flag',false);
myRule = myRule.buildLevels;
assert(myRule.levels(myRule.rules(myRule.search('rpts')).level).max==1,...
    'Failed to correctly build levels with false flags');

%test building levels when rules have a false flag when multiple rules reference the same level
%This is useful if you want to toggle between two different arrays of different length
myRule = myRule.modify('scheme.flag',false);
myRule = myRule.buildLevels;
testInd = myRule.search('test');
assert(myRule.levels(myRule.rules(testInd).level).max==length(myRule.rules(testInd).values),...
    'Failed to correctly build levels with false flags and rules referencing more than one level');

%% test level increment
[myRule,wl,rpt,schemes,gPos,desc] = testRule2();
myRule = myRule.modify('gpos.flag',false,'wl.flag',false);
myRule = myRule.buildLevels;

for ii = 1:length(schemes)
    for jj = 1:length(rpt)
        assert(all([myRule.levels.ind] == [ii,1,jj]),'Failed to increment inside loop');
        [~, myRule] = myRule.buildName('autoIncrement',true,'delimiter',', ');
    end
end
assert(all([myRule.levels.ind] == [1,1,1]),'Failed to resent increment values');

%% test name generation
[myRule,wl,rpt,schemes,gPos,desc] = testRule2();

%test call for 1st index
myRule = myRule.buildLevels;
[myName, myRule] = myRule.buildName('autoIncrement',true,'delimiter',', ');
assert(strcmp(myName,'sample A pH 13, scheme A, 400 nm, rpt 1, 400 nm'),'Failed to build correct name');
assert(all([myRule.levels.ind] == [1,1,1,2]), 'Failed to autoincrement');

%test call after setting different index and with rules set to false
myRule = myRule.modify('rpts.flag',false,'scheme.flag',false);
myRule = myRule.buildLevels;
myRule = myRule.setLevelIndex([1,2,1,2]);
[myName, myRule] = myRule.buildName('autoIncrement',true,'delimiter',', ');
assert(strcmp(myName,'sample A pH 13, 425 nm, 410 nm'),'Failed to build correct name');
assert(all([myRule.levels.ind] == [1,2,1,3]), 'Failed to autoincrement');

%% All tests pass
disp('All tests passed!');

%% functions specific to test script
function myRule = testRule1()
%generate a nameRule for testing methods
    myRule = nameRule();   %default constructor

    %add using cell array. This also tests adding using struct because of the recursive call
    myRule = myRule.add({'asdf',{'1','2','3'},1,true;...
                       'ghjk',{'4'},5,false});
end

function [myRule,wl,rpts,schemes,gPos,desc] = testRule2()
%generate a nameRule for testing execution
    wl = 400:10:450;
    rpts = 1:10;
    schemes = {'scheme A','scheme B','scheme C'};
    gPos = [400, 425, 450];
    desc = 'sample A pH 13';
    
    myRule = nameRule({'desc',  {desc},     0,true;...
                       'scheme',schemes(:), 1,true;...
                       'gpos',  strcat(strtrim(cellstr(num2str(gPos(:)))), ' nm'), 2,true;...
                       'rpts',  strcat({'rpt '}, strtrim(cellstr(num2str(rpts(:))))), 3,true;...
                       'wl',    strcat(strtrim(cellstr(num2str(wl(:)))), ' nm'),   4,true});
                   
                   
end