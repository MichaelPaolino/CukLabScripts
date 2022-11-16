%% Test class
l = 375:700;
t = -5:120;
load(fullfile(repoPath(),'CukResearchGroup','TR-VIS','chirpFitParam22Jul28.mat'));
myCAWs = STOCAWs(l,t,chirpFit,540);
myCAWs = myCAWs.setp(0.05/100,15,1.22,0.4);

% Make sure these params are reasonable
figure;
subplot(1,3,1)
    plot(l,myCAWs.amp);
    xlim([475,700]);
    ylim([0,0.16]);
subplot(1,3,2)
    plot(l,180*myCAWs.phase/pi);
    ylim([120,180]);
subplot(1,3,3)
    plot(l,myCAWs.freq);

figure;
contourf(l,t,myCAWs.M',23);
