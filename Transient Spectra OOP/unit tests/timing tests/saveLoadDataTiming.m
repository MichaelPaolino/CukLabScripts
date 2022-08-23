% Conclusions: Matlab build-in file saving is slow because of compression.
% Compression does not offer a huge size advantage for random numeric
% data. Normally compression would offer a speed boost because read/write
% is slow for hard disk data, but if compression does not reduce the size,
% it translates to a time cost. The tools available for chunked saving are 
% also fairly slow for structs--may require developing custom hdf5 format 
% for optimal performace. 
%
% For loading in the TR class, allow user to specify custom load string, 
% otherwise attempt to auto-determine load format retain compatibility with
% existing scripts. Note: hdf5 is Igor compatible and contains APIs in many
% languages.

%% Benchmark matfile save -v7.3
myStruct = genStructData(1000);
tic;
myMFile = matfile('testWriteMatFile.mat','Writable',true);
myMFile.myStruct = myStruct; %78.7 MB
toc %~3.3s

%% Benchmark save -v7.3 compression
myStruct = genStructData(1000);
tic;
save('testWriteCompression.mat','myStruct','-v7.3'); %78.7 MB
toc %~3.3s

%% Benchmark save -v7.3 no compression
myStruct = genStructData(1000);
tic;
save('testWriteNoCompression.mat','myStruct','-v7.3','-nocompression'); %82.6 MB
toc %~0.33s

%% Benchmark chunksave with matfile -v7.3
chunkSize = 10;
numChunks = 10;
pauseTime = 2; %in seconds
myStruct = genStructData(chunkSize*numChunks);

tic;
myMFile = matfile('testWriteMatFileChunked.mat','Writable',true);

for chunkInd = 1:chunkSize:(numChunks*chunkSize)
    myMFile.myStruct(chunkInd:(chunkInd+chunkSize-1),1) = myStruct(chunkInd:(chunkInd+chunkSize-1),1);
    pause(pauseTime);   %Allows background file processes to complete?
end

totalTime = toc; %1.6 s
procTime = totalTime - pauseTime*numChunks;
disp(['Benchmark chunksave with matfile -v7.3: ' num2str(procTime) ' s']);

%% Compare against loop of saving full dataset with matfile -v7.3
% Surprisingly, this is much slower than one would anticipate from the first
% three benchmarks... It is not clear if there is dead time due to the OS
% doing file processing

chunkSize = 10;
numChunks = 10;
pauseTime = 0; %in seconds

myStruct = genStructData(chunkSize*numChunks);
myBlankStruct(size(myStruct,1),1) = myStruct(end,1);

tic;

for chunkInd = 1:chunkSize:(numChunks*chunkSize)
    save('testWriteNoCompression.mat','myStruct','-v7.3','-nocompression');
    pause(pauseTime);   %Allows background file processes to complete?
end

totalTime = toc; %4.6 s !?!?
procTime = totalTime - pauseTime*numChunks;
disp(['Benchmark chunksave with matfile -v7.3: ' num2str(procTime) ' s']); 

%% functions
function dataOut = genStructData(len)
    dataOut = struct('a',rand(100,100),'b',rand(20,20));
    dataOut = repmat(dataOut,len,1);
end