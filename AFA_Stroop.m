%Stroop
clc
close all
clear all
patcon = {'AFA_pat_Stroop_','AFA_con_Stroop_'};
patconLength = [40];
for groupID = [1:2]; %1 pateints, 2 controls
for subID =  [1:20];%[number of particpants in each group i.e. N=20]
snirf = SnirfLoad([patcon{groupID} sprintf('%02d',subID) '.snirf']);

[patcon{groupID} sprintf('%02d',subID) ]

procResult.mlActMan{1} = ones(size(snirf.GetMeasList,1),1);
procResult.tInc{1}= ones(size(snirf.Get_t,1),1);

mlActAuto = hmrR_PruneChannels(snirf.data, snirf.probe ,procResult.mlActMan , procResult.tInc, [1e-03  1e+0 ], 15, [0.0 40.0]);
% mlActAuto = hmrR_PruneChannels(data, probe, mlActMan, tInc, dRange, SNRthresh, SDrange)

nSV = 0;
dod = hmrR_Intensity2OD( snirf.data );
[dod3, svs, nSV] = hmrR_PCAFilter( dod, mlActAuto, procResult.tInc, nSV );
dod2 = hmrR_BandpassFilt(dod3, 0.05, 0.5);

tMotion = 0.5;% Units of seconds.    Typical value ranges from 0.1 to 0.5.
tMask = 1; %Units of seconds. Typical value ranges from 0.5 to 1
STDEVthresh = 30;
AMPthresh =  1;
procResult.tIncAuto = hmrR_MotionArtifact(snirf.data, snirf.probe, procResult.mlActMan, mlActAuto, procResult.tInc, tMotion, tMask, STDEVthresh, AMPthresh);

tRange =[-5 5];
[procResult.stim, tRange] = hmrR_StimRejection(snirf.data, snirf.stim, procResult.tIncAuto, procResult.tInc, tRange);

XX = []; %this is based on age for each individual inorder to get the pathlanget factor. This is take from Scholkmann F, Wolf M.
%General equation for the differential pathlength factor of the 499 frontal human head depending on wavelength and age. Journal
%of biomedical optics.500 2013;18(10):105004- 
 dc = hmrR_OD2Conc( dod2, snirf.probe, [XX(subID,:)] );

% data - this is the concentration data with dimensions #time points x [HbO/HbR/HbT] x #channels
% stim - stimulation vector (# time points x #conditions)=1 at stim onset otherwise =0
% probe - source detector stucture (units should be consistent with rhoSD_ssThresh)
procResult.mlActAuto = mlActAuto;
procResult.Aaux = [];
rcMap = [];%snirf.aux;
trange = [-2 12 ]; 
glmSolveMethod = 1;
idxBasis = 1;
paramsBasis = [0.5 0.5];%[0.1 3.0 1.8 3.0];%[0.1 3.0 1.8 3.0];
rhoSD_ssThresh = 15;
flagNuisanceRMethod = 2;
driftOrder = 3;
c_vector = 0;
[dcAvg, dcAvgStd, nTrials, dcNew, dcResid, dcSum2, beta, R, hmrstats] = hmrR_GLM(dc, procResult.stim, snirf.probe, procResult.mlActAuto, procResult.Aaux, procResult.tInc, rcMap, trange, glmSolveMethod, idxBasis, paramsBasis, rhoSD_ssThresh, flagNuisanceRMethod, driftOrder, c_vector);
% [data_yavg, data_yavgstd, nTrials, data_ynew, data_yresid, data_ysum2, beta_blks, yR_blks, hmrstats] = ...
% hmrR_GLM(data_y, stim, probe, mlActAuto, Aaux, tIncAuto, rcMap, trange, glmSolveMethod, idxBasis, paramsBasis, rhoSD_ssThresh, flagNuisanceRMethod, driftOrder, c_vector)

result = dcAvg.GetDataTimeSeries('reshape');

oxy1{groupID}{subID} = squeeze(result(:,1,:,1));
oxy2{groupID}{subID} = squeeze(result(:,1,:,2));


% hh1 = figure; %Array figure
% plotProbe_old(squeeze(result(:,:,:,1)), dcAvg.time, snirf.Get_SD, hh1, [], [1 1]);


%ROI is Regions of intrestes
F1 = [1 2 3 4]; %R-FPA
F2 = [9 10 11 38]; %L-FPA 
F3 = [30 31 32 34]; %R-aVL
F4 = [13 14 15 24]; %L-aVL 
F5 = [33 35 46 47]; %R-pVL 
F6 = [23 25 27 28]; %L-pVLnfr
F7= [5 6 16 17 19 29]; %R-aDL
F8 = [12 39 40 41 42 44]; %L-aDL
F9 = [18 20 21 45]; %R-pDL
F10 = [26 43 50 51]; %L-pDL


F1_A = squeeze(result(:,1,F1,1)) * 1000000;
F1_A(F1_A==0)=NaN;
F1_A= nanmean(F1_A')' ;
time_index = find(F1_A==nanmax(F1_A(40:130)));
F1_A = nanmean(F1_A(time_index-10:time_index+10));

F2_A = squeeze(result(:,1,F2,1)) * 1000000;
F2_A(F2_A==0)=NaN;
F2_A= nanmean(F2_A')' ;
time_index = find(F2_A==nanmax(F2_A(40:130)));
F2_A = nanmean(F2_A(time_index-10:time_index+10));

F3_A = squeeze(result(:,1,F3,1)) * 1000000;
F3_A(F3_A==0)=NaN;
F3_A= nanmean(F3_A')' ;
time_index = find(F3_A==nanmax(F3_A(40:130)));
F3_A = nanmean(F3_A(time_index-10:time_index+10));

F4_A = squeeze(result(:,1,F4,1)) * 1000000;
F4_A(F4_A==0)=NaN;
F4_A= nanmean(F4_A')' ;
time_index = find(F4_A==nanmax(F4_A(40:130)));
F4_A = nanmean(F4_A(time_index-10:time_index+10));

F5_A = squeeze(result(:,1,F5,1)) * 1000000;
F5_A(F5_A==0)=NaN;
F5_A= nanmean(F5_A')' ;
time_index = find(F5_A==nanmax(F5_A(40:130)));
F5_A = nanmean(F5_A(time_index-10:time_index+10));

F6_A = squeeze(result(:,1,F6,1)) * 1000000;
F6_A(F6_A==0)=NaN;
F6_A= nanmean(F6_A')' ;
time_index = find(F6_A==nanmax(F6_A(40:130)));
F6_A = nanmean(F6_A(time_index-10:time_index+10));


F7_A = squeeze(result(:,1,F7,1)) * 1000000;
F7_A(F7_A==0)=NaN;
F7_A= nanmean(F7_A')' ;
time_index = find(F7_A==nanmax(F7_A(40:130)));
F7_A = nanmean(F7_A(time_index-10:time_index+10));

F8_A = squeeze(result(:,1,F8,1)) * 1000000;
F8_A(F8_A==0)=NaN;
F8_A= nanmean(F8_A')' ;
time_index = find(F8_A==nanmax(F8_A(40:130)));
F8_A = nanmean(F8_A(time_index-10:time_index+10));

F9_A = squeeze(result(:,1,F9,1)) * 1000000;
F9_A(F9_A==0)=NaN;
F9_A= nanmean(F9_A')' ;
time_index = find(F9_A==nanmax(F9_A(40:130)));
F9_A = nanmean(F9_A(time_index-10:time_index+10));

F10_A = squeeze(result(:,1,F10,1)) * 1000000;
F10_A(F10_A==0)=NaN;
F10_A= nanmean(F10_A')' ;
time_index = find(F10_A==nanmax(F10_A(40:130)));
F10_A = nanmean(F10_A(time_index-10:time_index+10));


F1_B = squeeze(result(:,1,F1,2)) * 1000000;
F1_B(F1_B==0)=NaN;
F1_B= nanmean(F1_B')' ;
time_index = find(F1_B==nanmax(F1_B(40:130)));
F1_B = nanmean(F1_B(time_index-10:time_index+10));

F2_B = squeeze(result(:,1,F2,2)) * 1000000;
F2_B(F2_B==0)=NaN;
F2_B= nanmean(F2_B')' ;
time_index = find(F2_B==nanmax(F2_B(40:130)));
F2_B = nanmean(F2_B(time_index-10:time_index+10));

F3_B = squeeze(result(:,1,F3,2)) * 1000000;
F3_B(F3_B==0)=NaN;
F3_B= nanmean(F3_B')' ;
time_index = find(F3_B==nanmax(F3_B(40:130)));
F3_B = nanmean(F3_B(time_index-10:time_index+10));

F4_B = squeeze(result(:,1,F4,2)) * 1000000;
F4_B(F4_B==0)=NaN;
F4_B= nanmean(F4_B')' ;
time_index = find(F4_B==nanmax(F4_B(40:130)));
F4_B = nanmean(F4_B(time_index-10:time_index+10));

F5_B = squeeze(result(:,1,F5,2)) * 1000000;
F5_B(F5_B==0)=NaN;
F5_B= nanmean(F5_B')' ;
time_index = find(F5_B==nanmax(F5_B(40:130)));
F5_B = nanmean(F5_B(time_index-10:time_index+10));

F6_B = squeeze(result(:,1,F6,2)) * 1000000;
F6_B(F6_B==0)=NaN;
F6_B= nanmean(F6_B')' ;
time_index = find(F6_B==nanmax(F6_B(40:130)));
F6_B = nanmean(F6_B(time_index-10:time_index+10));


F7_B = squeeze(result(:,1,F7,2)) * 1000000;
F7_B(F7_B==0)=NaN;
F7_B= nanmean(F7_B')' ;
time_index = find(F7_B==nanmax(F7_B(40:130)));
F7_B = nanmean(F7_B(time_index-10:time_index+10));

F8_B = squeeze(result(:,1,F8,2)) * 1000000;
F8_B(F8_B==0)=NaN;
F8_B= nanmean(F8_B')' ;
time_index = find(F8_B==nanmax(F8_B(40:130)));
F8_B = nanmean(F8_B(time_index-10:time_index+10));

F9_B = squeeze(result(:,1,F9,2)) * 1000000;
F9_B(F9_B==0)=NaN;
F9_B= nanmean(F9_B')' ;
time_index = find(F9_B==nanmax(F9_B(40:130)));
F9_B = nanmean(F9_B(time_index-10:time_index+10));

F10_B = squeeze(result(:,1,F10,2)) * 1000000;
F10_B(F10_B==0)=NaN;
F10_B= nanmean(F10_B')' ;
time_index = find(F10_B==nanmax(F10_B(40:130)));
F10_B = nanmean(F10_B(time_index-10:time_index+10));


oxy_peak_neutral{groupID}{subID} = [F1_A F2_A F3_A F4_A F5_A F6_A F7_A F8_A F9_A F10_A]'; %all peaks 
oxy_peak_incongruent{groupID}{subID} = [F1_B F2_B F3_B F4_B F5_B F6_B F7_B F8_B F9_B F10_B]'; %all peaks 

 
All_oxy{groupID}{subID} = squeeze(result(:,1,:,:));
All_deoxy{groupID}{subID} = squeeze(result(:,2,:,:));


end
end
