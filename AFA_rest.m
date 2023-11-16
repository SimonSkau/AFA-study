
clear all
close all

patcon = {'AFA_pat_rest_','AFA_con_rest_'};
patconLength = [20, 20];


for subType = 1:2;
for subID = 1:patconLength(subType); 
load([patcon{subType} num2str(subID) '.nirs'],'-mat');

[subType subID ]
ncomp = 0;
%Run basic processing stream
trange = [-5 300]; 
fs = length(t)/max(t);
%trange = [-2 12];
fs = length(t)/max(t);
tIncMan = ones(size(d(:,1))); 
procResult.dod = hmrIntensity2OD(d);

procResult.SD = enPruneChannels(d,SD,ones(size(t)),[1e-4  1e+0 ],15,[0.0 40.0],3);

procResult.SD.MeasListAct(end/2+1:end) = procResult.SD.MeasListAct(1:end/2); %This shouldn't be necessary, not sure why it is

procResult.dod = hmrBandpassFilt(procResult.dod,fs,0.01,0.5);

procResult.tInc = ones(size(t));
%Try removing 1st principal component from whole dataset

procResult.dod = enPCAFilter(procResult.dod,procResult.SD,procResult.tInc,[ncomp ncomp]);

STDEVthresh = 30;
AMPthresh =  1;
[procResult.tInc, procResult.tIncCh] = hmrMotionArtifactByChannel(procResult.dod, fs, procResult.SD, ones(length(d),1), 1, 5, STDEVthresh, AMPthresh); %procResult.tInc was before procResult.SD
%[procResult.tIncPCA, procResult.tIncChPCA] = hmrMotionArtifactByChannel(dod_PCA, fs, procResult.SD, ones(length(d),1), 1, 5, STDEVthresh, AMPthresh); %procResult.tInc was before procResult.SD

p = 0.99; %recommended in the littreture, -1 to not use spline
procResult.dod = hmrMotionCorrectSpline(procResult.dod, t, procResult.SD, [procResult.tInc, procResult.tIncCh], p); %procResult.tIncCh
%procResult.dodSplinePCA = hmrMotionCorrectSpline(dod_PCA, t, procResult.SD, [procResult.tInc, procResult.tIncCh], p); %procResult.tIncCh

%[procResult.dodPCA,svs,nSV] = hmrMotionCorrectPCA(procResult.SD, procResult.dodSpline, procResult.tInc, nSV );

%procResult.dodPCA = hmrBandpassFilt(procResult.dodPCArecursePCA,fs,0.01,0.2);

XX = []; % this is based on age for each individual inorder to get the pathlanget factor. Values are not added here due to anonymity. This is take from Scholkmann F, Wolf M.
%General equation for the differential pathlength factor of the 499 frontal human head depending on wavelength and age. Journal
%of biomedical optics.500 2013;18(10):105004- 

%Convert to conc
procResult.dc = hmrOD2Conc(procResult.dod,procResult.SD, [XX(subID,:)]); 
%procResult.dcPCA = hmrOD2Conc(procResult.dodPCA,procResult.SD,[6 6]);

% y_reg = hmrSSRegressionByChannel(y, SD, rhoSD_ssThresh, flagSSmethod)
rhoSD_ssThresh = 16;
flagSSmethod = 1; 
procResult.dc  = hmrSSRegressionByChannel(procResult.dc, procResult.SD, rhoSD_ssThresh, flagSSmethod);
%dc_PCA  = hmrSSRegressionByChannel(dc_PCA, SD, rhoSD_ssThresh, flagSSmethod);

%Calc HRFs using simple block average;
[yHRF_BA, ~, tHRF] = hmrBlockAvg(procResult.dc, s, t, trange);
%[yHRF_BA_Regressed, ~, tHRF] = hmrBlockAvg(dc_Regressed, s, t, trange);
%yHRF_BA_PCA = enPCAFilter(squeeze(yHRF_BA(:,:,:,1)),procResult.SD,ones(size(tHRF)),[1 1 1]);


% hh1 = figure; 
% plotProbe_old(squeeze(yHRF_BA(:,:,:)), tHRF, SD, hh1, [], [1 1]);



%channels of interest n = 44
coi = [1,2,3,4,5,6,9,10,11,12,13,14,15,16,17,18,19,20,21,23,24,25,26,27,28,29,30,31,32,33,34,35,38,39,40,41,42,43,44,45,46,47,50,51];

all_oxy{subType}{subID} = yHRF_BA(:,1,coi);

end
end


%%
