%this file contains unit tests for the doubleWUnits class

%% CONSTRUCTOR METHODS
    %% default constructor
    myArray = doubleWithUnits();

    %% construct an arrayWithUnits without any unit rules
    myArray = doubleWithUnits(400:450,'nm','Wavelength (nm)');

%% ACCESSOR METHODS (DEPENDENT)
    myArray = doubleWithUnits(400:450,'nm','Wavelength (nm)');
    assert(strcmp(myArray.unit,'nm'),'get.unit failed');
    assert(strcmp(myArray.base,'nm'),'get.base failed');
    
    myArray = myArray.changeBaseName('nanometer','Wavelength (nanometer)');
    assert(strcmp(myArray.base,'nanometer'),'changeBaseName failed');

%% SUBCLASS METHODS
    %% test adding unit rules and converting to a unit other than the base unit
    l_0 = 400; 
    myArray = testArray();
    myArray = myArray.convert('cm-1');
    testVal = double(myArray);
    assert(testVal(1) == 0, 'conversion from nm to cm-1 failed');
    
    %% test alternative convert unit method set.unit
    l_0 = 400; 
    myArray = testArray();
    myArray.unit = 'cm-1';
    testVal = double(myArray);
    assert(testVal(1) == 0, 'conversion from nm to cm-1 failed');
    
    %% test conversion of a non-base unit to another non-base unit
    %Note: run previous test
    myArray = myArray.convert('um');
    testVal = double(myArray);
    assert(testVal(1) == 0.4, 'conversion from cm-1 to um failed');

    %% test changing a rule that is not the current unit (i.e. change raman pump wavelength)
    %Note: run previous test
    l_0 = 375;
    myArray = myArray.updateRule('cm-1','Raman Shift (cm^{-1})',... %wavelength to raman shift (nm to cm-1) rule 
                                 @(f) 1e7*(1/l_0-1./f),...  %nm to cm-1
                                 @(f) 1./(1/l_0-1e-7*f));   %cm-1 to nm
    myArray = myArray.convert('cm-1');
    testVal = double(myArray);
    assert(testVal(1) == 1e7*(1/375-1/400), 'updating cm-1 failed when current unit is different than updated unit');

    %% test changing a rule that is the current unit (i.e. change raman pump wavelength)
    %Note: run previous test
    myArray = myArray.convert('cm-1');
    l_0 = 400;
    myArray = myArray.updateRule('cm-1','Raman Shift (cm^{-1})',... %wavelength to raman shift (nm to cm-1) rule 
                                 @(f) 1e7*(1/l_0-1./f),...  %nm to cm-1
                                 @(f) 1./(1/l_0-1e-7*f));   %cm-1 to nm
    testVal = double(myArray);
    assert(testVal(1) == 0, 'updating cm-1 failed when current unit is the same as updated unit');

%% DOUBLE SUPERCLASS METHODS
    %% SUBSREF AND ?? METHODS

        %% test superclass subref calls for subclass properties
        myArray = testArray();
        assert(strcmp(myArray.unit,'nm'),'Failed to make superclass subsref property call');
        assert(strcmp(myArray(1).unit,'nm'),'Failed to make superclass subsref property call');

        %% test superclass subref calls for subclass properties that are arrays
        myArray = testArray();
        myUnitRules = myArray.unitRules;
        assert(isa(myUnitRules,'struct'),'Failed to return unitRules struct');
        assert(length(myUnitRules)==3,'Failed to return all elements of unitRules struct');

        myUnitRules = myArray.unitRules(1);
        assert(isa(myUnitRules,'struct'),'Failed to return unitRules struct on a subset call');
        assert(length(myUnitRules)==1,'Failed to return all elements of unitRules struct');

        %% test array subset calls with different return types
        myArray = testArray();
        assert(isa(myArray(1:10),'doubleWithUnits'),'Array subset call did not return doubleWithUnits class');
        assert(isa(myArray.data,'double'),'Array data call did not return double class');
        assert(isa(myArray.data(1:10),'double'),'Array data subset call did not return double class');
        
        %% test superclass subref calls for subclass methods
        myArray = testArray();
        myArray = myArray(1:10).convert('cm-1');
        assert(length(myArray) == 10,'Array subset selection failed');
        assert(myArray(1) == 0,'Array conversion on array subset failed');
        
        %% test superclass subref assignment calls
        myArray = testArray();
        myArray.data = 400:5:450;
        myArray.data(1) = 0;
        myArray(1:2) = [1 2];
        %myArray(1:2) = myArray(3:4);    %!!!fails!!!

    %% ARRAY CONCATONATION
        %% test array concatonation methods
        myArray = doubleWithUnits(400:450,'nm','Wavelength (nm)');
        myArray = [myArray myArray];
        assert(all(size(myArray) == [1 102]),'horzcat failed');
        myArray = [myArray; myArray];
        assert(all(size(myArray) == [2 102]),'vercat failed');

        %% test array concatonation methods with different units
        myArray = testArray();

        myArray = [myArray.convert('cm-1'); myArray.convert('um')];
        assert(all(size(myArray) == [2 51]),'vercat failed with different units');
        testVal = double(myArray);
        assert((testVal(1)) < 1e-6,'vercat failed to convert units');
        
        %% test array concatonation method with mixed doubleWithUnits and double classes
        myArray = testArray();
        %myArray = [0 myArray];  %!!!fails!!!

%% TEST ARRAY FUNCTIONS
function myArray = testArray()
    l_0 = 400; 
    myArray = doubleWithUnits(400:450,'nm','Wavelength (nm)');
    myArray = myArray.addRule('cm-1','Raman Shift (cm^{-1})',@(f) 1e7*(1/l_0-1./f),@(f) 1./(1/l_0-1e-7*f));    %wavelength to raman shift (nm to cm-1) rule 
    myArray = myArray.addRule('um','Wavelength (\mm)',@(f) 1e-3*f,@(f) 1e3*f);    %wavelength to raman shift (nm to cm-1) rule 
end