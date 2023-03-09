function processData
% This script performs frequency analysis, generates spectra and calculates
% measures of eye movements (data for figure 4)
%
% Companion code for:
%
% N-methyl d-aspartate receptor hypofunction reduces steady state visual
% evoked potentials (2023)

%where are data located and where should results be saved
sourceFolder = strrep(pwd,'code','data\');
dataFolder = [sourceFolder 'combined\'];
targetFolder = [sourceFolder 'processed\'];

%load useInfo
useInfo = load([sourceFolder 'processed\' 'useInfo']);
useInfo = useInfo.useInfo;


%% frequency analysis

%where are files with signal located
dataFiles = dir(dataFolder);
dataFiles = {dataFiles.name};
dataFiles(1:2) = [];

%information about trials and resulting freqeuncy resolution
time = -800:2600;   %trial aligned to stimulus onset

maxSamplePoints = max(2.^nextpow2(sum(time>=501 & time<=1500)));
signalLength = sum(time>=501 & time<=1500);  %use center 1000 ms of stimulus presentation
frequency = maxSamplePoints /(signalLength/1000)*(0:(maxSamplePoints/2))/maxSamplePoints; %frequency vector
frequency(frequency>maxSamplePoints/2) = [];


%load files
catchCntr = 0;
totalElectrodeCntr = 0;
for fileCntr = 1:length(dataFiles)
    
    %load data
    tempFile = load([dataFolder dataFiles{fileCntr}]);
    tempFile = tempFile.data;
    
    uCond = unique(tempFile.lfp.trialInfo.conIdent{1});
    nrElectrodes = size(tempFile.lfp.signal{1},3);
    for electrodeCntr = 1:nrElectrodes
        totalElectrodeCntr = totalElectrodeCntr +1;
        
        for conditionCntr = 1:numel(uCond)

            try
                %select trials
                conIdent = useInfo.trialSelection.saline{totalElectrodeCntr,conditionCntr+1};
                tempSignalSaline = tempFile.lfp.signal{1}(time>=501 & time<=1500,conIdent,electrodeCntr);
                tempSignalSalineMT = tempFile.lfp.signal{1}(:,conIdent,electrodeCntr);
                conIdent = useInfo.trialSelection.ketamine{totalElectrodeCntr,conditionCntr+1};
                tempSignalKetamine = tempFile.lfp.signal{2}(time>=501 & time<=1500,conIdent,electrodeCntr);
                tempSignalKetamineMT = tempFile.lfp.signal{2}(:,conIdent,electrodeCntr);
                
                %resample signal 
                tempSignalSaline = resample(tempSignalSaline,maxSamplePoints,signalLength);
                tempSignalKetamine = resample(tempSignalKetamine,maxSamplePoints,signalLength);

                %calculate evoked power
                fftSalineTemp = fft(mean(tempSignalSaline,2,'omitnan'));  
                fftSalineTemp = abs(fftSalineTemp);                            	%magnitude
                fftSalineTemp = fftSalineTemp./maxSamplePoints;                 %correct spectrum for number of samples
                fftSalineTemp = fftSalineTemp(1:maxSamplePoints/2+1,:);
                fftSalineTemp(2:end-1,:) = 2* fftSalineTemp(2:end-1,:);

                %repeat for ketamine
                fftKetamineTemp = fft(mean(tempSignalKetamine,2,'omitnan'));     
                fftKetamineTemp = abs(fftKetamineTemp);                         %magnitude
                fftKetamineTemp = fftKetamineTemp./maxSamplePoints;           	%correct spectrum for number of samples
                fftKetamineTemp = fftKetamineTemp(1:maxSamplePoints/2+1,:);
                fftKetamineTemp(2:end-1,:) = 2* fftKetamineTemp(2:end-1,:);

                tempEvokedSaline = fftSalineTemp;
                tempEvokedKetamine = fftKetamineTemp;

                dividedSignal.saline.evoked(:,conditionCntr,totalElectrodeCntr) = tempEvokedSaline(frequency<=125);
                dividedSignal.ketamine.evoked(:,conditionCntr,totalElectrodeCntr) = tempEvokedKetamine(frequency<=125);


                %caclulate induced power using the chronux toolbox
                mtParms.trialave   	= 0;
                mtParms.pad       	= -1;   %no padding, we resampled to ^2 instead
                mtParms.Fs        	= maxSamplePoints;
                mtParms.tapers     	= [2 3];
                mtParms.err       	= [1 0.05];
                mtParms.fpass     	= [0 125];
                mtParms.stepsize    = maxSamplePoints;
                mtParms.window      = maxSamplePoints;

                %total power     
                [mtSalineTemp,~,~,~,~] =  mtspecgramc(tempSignalSaline,[mtParms.window/mtParms.Fs mtParms.stepsize/mtParms.Fs],mtParms);
                [mtKetamineTemp,~,~,~,~] =  mtspecgramc(tempSignalKetamine,[mtParms.window/mtParms.Fs mtParms.stepsize/mtParms.Fs],mtParms);
                mtSalineTemp = sqrt(mtSalineTemp);      %magnitude
                mtKetamineTemp = sqrt(mtKetamineTemp);  %magnitude

                tempTotalSaline = mean(mtSalineTemp,2,'omitnan')';
                tempTotalKetamine = mean(mtKetamineTemp,2,'omitnan')';       

                %evoked power
                [mtSalineTemp,~,~,~,~] =  mtspecgramc(mean(tempSignalSaline,2,'omitnan'),[mtParms.window/mtParms.Fs mtParms.stepsize/mtParms.Fs],mtParms);
                [mtKetamineTemp,~,~,~,~] =  mtspecgramc(mean(tempSignalKetamine,2,'omitnan'),[mtParms.window/mtParms.Fs mtParms.stepsize/mtParms.Fs],mtParms);
                mtSalineTemp = sqrt(mtSalineTemp);
                mtKetamineTemp = sqrt(mtKetamineTemp);

                tempEvokedSaline = mtSalineTemp;
                tempEvokedKetamine = mtKetamineTemp;

                %calculate induced power by subtracting 
                tempInducedSaline = tempTotalSaline-tempEvokedSaline;
                tempInducedKetamine = tempTotalKetamine-tempEvokedKetamine;

                dividedSignal.saline.induced(:,conditionCntr,totalElectrodeCntr) = tempInducedSaline;
                dividedSignal.ketamine.induced(:,conditionCntr,totalElectrodeCntr) = tempInducedKetamine;


                %calculate baseline for total power using the chronux toolbox
                mtParms.trialave   	= 0;
                mtParms.pad       	= -1;
                mtParms.Fs        	= 1000;
                mtParms.tapers     	= [2 3];
                mtParms.err       	= [1 0.05];
                mtParms.fpass     	= [0 125];
                mtParms.stepsize    = 50;
                mtParms.window      = 600;


                %total power        
                [mtSalineTemp,mtTime,mtFreq,~,~] =  mtspecgramc(tempSignalSalineMT,[mtParms.window/mtParms.Fs mtParms.stepsize/mtParms.Fs],mtParms);
                [mtKetamineTemp,~,~,~,~] =  mtspecgramc(tempSignalKetamineMT,[mtParms.window/mtParms.Fs mtParms.stepsize/mtParms.Fs],mtParms);
                mtTime = (mtTime-0.801)*1000;
                mtSalineTemp = sqrt(mtSalineTemp);
                mtKetamineTemp = sqrt(mtKetamineTemp);  %use magnitude

                tempTotalSalineBaseline = squeeze(mean(mean(mtSalineTemp(mtTime>=-450 & mtTime<=-350,:,:),1,'omitnan'),3));
                tempTotalKetamineBaseline = squeeze(mean(mean(mtKetamineTemp(mtTime>=-450 & mtTime<=-350,:,:),1,'omitnan'),3));


                dividedSignal.saline.totalBaseline(:,conditionCntr,totalElectrodeCntr) = tempTotalSalineBaseline;
                dividedSignal.ketamine.totalBaseline(:,conditionCntr,totalElectrodeCntr) = tempTotalKetamineBaseline;

            catch
                catchCntr = catchCntr+1;

                nanPlaceHolder = nan(sum(frequency<=125),1);
                dividedSignal.saline.evoked(:,conditionCntr,totalElectrodeCntr) = nanPlaceHolder;
                dividedSignal.ketamine.evoked(:,conditionCntr,totalElectrodeCntr) = nanPlaceHolder;
                dividedSignal.saline.induced(:,conditionCntr,totalElectrodeCntr) = nanPlaceHolder;
                dividedSignal.ketamine.induced(:,conditionCntr,totalElectrodeCntr) = nanPlaceHolder;

                nanPlaceHolder = nan(numel(mtFreq),1);
                dividedSignal.saline.totalBaseline(:,conditionCntr,totalElectrodeCntr) = nanPlaceHolder;
                dividedSignal.ketamine.totalBaseline(:,conditionCntr,totalElectrodeCntr) = nanPlaceHolder;
            end
        end
    end
end

   
%save divided data
save([targetFolder 'dividedSignal.mat'],'dividedSignal','-v7.3');

               

%% permutation testing at each frequency for plots

%which electrodes to use
useSelection = useInfo.doUse;

for conditionCntr = 1:numel(uCond)
    for freqCntr = 1:sum(frequency<=125)

        % evoked power
        permResult = permutationTest2(squeeze(dividedSignal.ketamine.evoked(freqCntr,conditionCntr,useSelection(:,conditionCntr)))*(freqCntr-1),squeeze(dividedSignal.saline.evoked(freqCntr,conditionCntr,useSelection(:,conditionCntr)))*(freqCntr-1),'groupNames',{'ketamine';'saline'},'paired',true);

        splitSpectrum.evoked.condition(conditionCntr).difference.mean(freqCntr) = permResult.difference;
        splitSpectrum.evoked.condition(conditionCntr).difference.ci(freqCntr,:) = permResult.ci;
        splitSpectrum.evoked.condition(conditionCntr).difference.zScore(freqCntr) = permResult.zScore;
        splitSpectrum.evoked.condition(conditionCntr).difference.permPVal(freqCntr) = permResult.permPVal;
        splitSpectrum.evoked.condition(conditionCntr).saline.mean(freqCntr) = permResult.saline.mean;
        splitSpectrum.evoked.condition(conditionCntr).saline.ci(freqCntr,:) = permResult.saline.ci;
        splitSpectrum.evoked.condition(conditionCntr).ketamine.mean(freqCntr) = permResult.ketamine.mean;
        splitSpectrum.evoked.condition(conditionCntr).ketamine.ci(freqCntr,:) = permResult.ketamine.ci;


        % induced power
        permResult = permutationTest2(squeeze(dividedSignal.ketamine.induced(freqCntr,conditionCntr,useSelection(:,conditionCntr)))*(freqCntr-1),squeeze(dividedSignal.saline.induced(freqCntr,conditionCntr,useSelection(:,conditionCntr)))*(freqCntr-1),'groupNames',{'ketamine';'saline'},'paired',true);

        splitSpectrum.induced.condition(conditionCntr).difference.mean(freqCntr) = permResult.difference;
        splitSpectrum.induced.condition(conditionCntr).difference.ci(freqCntr,:) = permResult.ci;
        splitSpectrum.induced.condition(conditionCntr).difference.zScore(freqCntr) = permResult.zScore;
        splitSpectrum.induced.condition(conditionCntr).difference.permPVal(freqCntr) = permResult.permPVal;
        splitSpectrum.induced.condition(conditionCntr).saline.mean(freqCntr) = permResult.saline.mean;
        splitSpectrum.induced.condition(conditionCntr).saline.ci(freqCntr,:) = permResult.saline.ci;
        splitSpectrum.induced.condition(conditionCntr).ketamine.mean(freqCntr) = permResult.ketamine.mean;
        splitSpectrum.induced.condition(conditionCntr).ketamine.ci(freqCntr,:) = permResult.ketamine.ci;

    end
end


%save splitPower 
save([targetFolder 'splitSpectrum.mat'],'splitSpectrum','-v7.3');              

          

%% calculate fixation instability and and fixation inaccuracy
eyeDataPerTrial(length(dataFiles)).fixAccuracy.saline = NaN;

electrodeStartNr = 1;
for fileCntr = 1:length(dataFiles)
    
    %load data
    tempFile = load([dataFolder dataFiles{fileCntr}]);
    tempFile = tempFile.data;
    
   	displacementThreshold = 0.01;    %eyelink 1000 system has an accuracy of 0.01 degrees visual angle
    subject = tempFile.subject;
    
    electrodeStopNr = electrodeStartNr + size(tempFile.lfp.signal{1},3) -1;
    for conditionCntr = 1:numel(uCond)
        
        try
        %which trials fit criteria (i.e., no fixtion breaks, correct condition
        %and initiated within 60 minutes after injection)
        useTrialsSaline = tempFile.eye.trialInfo.useTrials{1} & tempFile.eye.trialInfo.realTime{1}<=60 & tempFile.eye.trialInfo.conIdent{1}==uCond(conditionCntr);
        useTrialsKetamine = tempFile.eye.trialInfo.useTrials{2} & tempFile.eye.trialInfo.realTime{2}<=60 & tempFile.eye.trialInfo.conIdent{2}==uCond(conditionCntr);
        

        if sum(isnan(tempFile.eye.signal{1}(time>500 & time<=1500,find(useTrialsSaline,1),1)))==0
            xCoordinatesSaline = tempFile.eye.signal{1}(time>500 & time<=1500,useTrialsSaline,1);
            yCoordinatesSaline = tempFile.eye.signal{1}(time>500 & time<=1500,useTrialsSaline,2);
            xCoordinatesKetamine = tempFile.eye.signal{2}(time>500 & time<=1500,useTrialsKetamine,1);
            yCoordinatesKetamine = tempFile.eye.signal{2}(time>500 & time<=1500,useTrialsKetamine,2);
        else
           	xCoordinatesSaline = tempFile.eye.signal{1}(time>500 & time<=1500,useTrialsSaline,4);
            yCoordinatesSaline = tempFile.eye.signal{1}(time>500 & time<=1500,useTrialsSaline,5);
            xCoordinatesKetamine = tempFile.eye.signal{2}(time>500 & time<=1500,useTrialsKetamine,4);
            yCoordinatesKetamine = tempFile.eye.signal{2}(time>500 & time<=1500,useTrialsKetamine,5);
        end
        

        distanceFromFixSaline  = mean(sqrt(xCoordinatesSaline.^2 + yCoordinatesSaline.^2));
        distanceFromFixKetamine  = mean(sqrt(xCoordinatesKetamine.^2 + yCoordinatesKetamine.^2));
        
        distanceTraveledSaline  = abs(diff(sqrt(xCoordinatesSaline.^2 + yCoordinatesSaline.^2)));
        distanceTraveledSaline(distanceTraveledSaline<displacementThreshold) = 0;
        distanceTraveledSaline = sum(distanceTraveledSaline);
        
        distanceTraveledKetamine  = abs(diff(sqrt(xCoordinatesKetamine.^2 + yCoordinatesKetamine.^2)));
        distanceTraveledKetamine(distanceTraveledKetamine<displacementThreshold) = 0;
        distanceTraveledKetamine = sum(distanceTraveledKetamine);
        
        eyeDataPerTrial(fileCntr).subject = subject;
        eyeDataPerTrial(fileCntr).fixInaccuracy.condition(conditionCntr).saline = distanceFromFixSaline;
        eyeDataPerTrial(fileCntr).fixInaccuracy.condition(conditionCntr).ketamine = distanceFromFixKetamine;
        eyeDataPerTrial(fileCntr).fixInstability.condition(conditionCntr).saline = distanceTraveledSaline;
        eyeDataPerTrial(fileCntr).fixInstability.condition(conditionCntr).ketamine = distanceTraveledKetamine;
        
        %put data in same format as magnitude measures for evoked, induced
        %and baseline activity
        eyeData.fixInaccuracy.saline(electrodeStartNr:electrodeStopNr,conditionCntr) = mean(eyeDataPerTrial(fileCntr).fixInaccuracy.condition(conditionCntr).saline);
        eyeData.fixInaccuracy.ketamine(electrodeStartNr:electrodeStopNr,conditionCntr) = mean(eyeDataPerTrial(fileCntr).fixInaccuracy.condition(conditionCntr).ketamine);
        eyeData.fixInstability.saline(electrodeStartNr:electrodeStopNr,conditionCntr) = mean(eyeDataPerTrial(fileCntr).fixInstability.condition(conditionCntr).saline);
        eyeData.fixInstability.ketamine(electrodeStartNr:electrodeStopNr,conditionCntr) = mean(eyeDataPerTrial(fileCntr).fixInstability.condition(conditionCntr).ketamine);
        
        catch
            eyeDataPerTrial(fileCntr).subject = subject;
            eyeDataPerTrial(fileCntr).fixInaccuracy.condition(conditionCntr).saline = NaN;
            eyeDataPerTrial(fileCntr).fixInaccuracy.condition(conditionCntr).ketamine = NaN;
            eyeDataPerTrial(fileCntr).fixInstability.condition(conditionCntr).saline = NaN;
            eyeDataPerTrial(fileCntr).fixInstability.condition(conditionCntr).ketamine = NaN;
            
        	eyeData.fixInaccuracy.saline(electrodeStartNr:electrodeStopNr,conditionCntr) = NaN;
            eyeData.fixInaccuracy.ketamine(electrodeStartNr:electrodeStopNr,conditionCntr) = NaN;
            eyeData.fixInstability.saline(electrodeStartNr:electrodeStopNr,conditionCntr) = NaN;
            eyeData.fixInstability.ketamine(electrodeStartNr:electrodeStopNr,conditionCntr) = NaN;
        end
        
    end
    electrodeStartNr = electrodeStopNr+1;
end
        
%eyeData
save([targetFolder 'eyeData.mat'],'eyeData','-v7.3');             
save([targetFolder 'eyeDataPerTrial.mat'],'eyeDataPerTrial','-v7.3');

end
