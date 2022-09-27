testData = wlTR('22-07-28_15h03m26s_OC_chirp_calibration_W_Water.mat','shortName','chirp');
testData = testData.fitChirp('order',7,'wavelengths',365:5:735,'delays',[-3,3]);

testData2 = wlTR('22-07-22_20h36m26s_Sample_W_pH_12100mMNa.mat','shortName','chirp');
testData2 = testData2.trim('delays',[-3,3]);
testData2.chirpParams = testData.chirpParams;

testData3 = testData2.correctChirp('interp','spline');