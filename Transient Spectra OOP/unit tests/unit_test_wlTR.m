%% Load test data set


%% findT0 tests
testData = wlTR('22-07-28_15h03m26s_OC_chirp_calibration_W_Water.mat','short name','chirp');

%flip sign
testData.spectra.data = -testData.spectra.data;

%test without stitching
[~,t0] = testData.findT0;
assert(all(size(t0)==[1,1,2]),'Incorrect t0 size');
assert(all(t0 > 0.5) & all(t0 < 0.7),'Could not find correct t0');

%test with stitching
testData = testData.stitch;

[~,t0] = testData.findT0;
assert(all(size(t0)==[1,1]),'Incorrect t0 size');
assert((t0 > 0.5) & (t0 < 0.7),'Could not find correct t0');

%test actual multi-d data with repeats
testData = wlTR('22-07-22_20h36m26s_Sample_W_pH_12100mMNa.mat','short name','pH 12');
[~,t0] = testData.findT0;
assert(all(size(t0)==[1,4,2]),'Incorrect t0 size');
assert(all(t0 > 0.5,'all') & all(t0 < 0.7,'all'),'Could not find correct t0');

%test with stitch and average
testData = testData.average;
testData = testData.stitch;
[~,t0] = testData.findT0;
assert(all(size(t0)==[1,1]),'Incorrect t0 size');
assert((t0 > 0.5) && (t0 < 0.7),'Could not find correct t0');

%% correctT0 tests
testData = wlTR('22-07-22_20h36m26s_Sample_W_pH_12100mMNa.mat','short name','chirp');
testData = testData.trim('delays',[-2,2]);

figure;
testData.plotKinetics('wavelengths',540);
hold on;

testData = testData.correctT0('average',false);
testData.plotKinetics('wavelengths',540);
hold off;

assert(all(size(testData.t0.data)==[1,4,2]),'Incorrect obj.t0 size');

testData = wlTR('22-07-22_20h36m26s_Sample_W_pH_12100mMNa.mat','short name','chirp');
testData = testData.trim('delays',[-2,2]);

figure;
testData.plotKinetics('wavelengths',540);
hold on;

testData = testData.correctT0();
testData.plotKinetics('wavelengths',540);
hold off;

assert(all(size(testData.t0.data)==[1,1]),'Incorrect obj.t0 size');

%% fit chirp tests
testData = wlTR('22-07-28_15h03m26s_OC_chirp_calibration_W_Water.mat','short name','chirp');

testData = testData.stitch();
testData = testData.subset('wavelengths',375:5:725);

[testData, chirpParam] = testData.fitChirp();

testData = wlTR('22-07-28_15h03m26s_OC_chirp_calibration_W_Water.mat','short name','chirp');

[testData, chirpParam] = testData.fitChirp('order',7,'wavelengths',365:5:735,'delays',[-3,3]);