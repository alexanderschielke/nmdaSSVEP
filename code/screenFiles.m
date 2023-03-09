function screenFiles
% This script performs quality checks to determine whether recordings sites
% should be included and counts number trials that were included or
% excluded. Additionally, signals for figures 2 and 3 are generated
%
% Companion code for:
%
% N-methyl d-aspartate receptor hypofunction reduces steady state visual
% evoked potentials (2023)
% Alexander Schielke & Bart Krekelberg
% Center for Molecular and Behavioral Neuroscience
% Rutgers University - Newark 

%where are data located and where should results be saved
sourceFolder = strrep(pwd,'code','data\combined\');
targetFolder = strrep(pwd,'code','data\processed\');

fileNames = dir(sourceFolder);
fileNames = {fileNames.name};
fileNames(1:2) = [];
time = -800:2600;
%by file
totalElectrodeCntr = 0;
for fileCntr = 1:length(fileNames)

    tempFile = load([sourceFolder fileNames{fileCntr}]);
    tempFile = tempFile.data;
    uCond = unique(tempFile.lfp.trialInfo.conIdent{1});
    fileName = strrep(fileNames{fileCntr},'.mat','');
    fileName = str2double(fileName(5:end));

%where 
    for electrodeCntr = 1:size(tempFile.lfp.signal{1},3)
        totalElectrodeCntr = totalElectrodeCntr+1;

        %track trials excluded because of
        %   -blinks
        %   -missing data
        %   -artifacts
        eyeClosedSaline = tempFile.lfp.trialInfo.outlierInfo{1}.electrode(electrodeCntr).eyesClosed;
        nanTrialsSaline = tempFile.lfp.trialInfo.outlierInfo{1}.electrode(electrodeCntr).nanTrials;
        artifactsSaline = tempFile.lfp.trialInfo.outlierInfo{1}.electrode(electrodeCntr).artifact;
        eyeClosedKetamine = tempFile.lfp.trialInfo.outlierInfo{2}.electrode(electrodeCntr).eyesClosed;
        nanTrialsKetamine = tempFile.lfp.trialInfo.outlierInfo{2}.electrode(electrodeCntr).nanTrials;
        artifactsKetamine = tempFile.lfp.trialInfo.outlierInfo{2}.electrode(electrodeCntr).artifact;

        uEyeClosedSaline = eyeClosedSaline;
        uNanTrialsSaline = nanTrialsSaline;
        uNanTrialsSaline(uEyeClosedSaline) = 0;
        uArtifactsSaline = artifactsSaline;
        uArtifactsSaline(uEyeClosedSaline | uNanTrialsSaline) = 0;
        uEyeClosedKetamine = eyeClosedKetamine;
        uNanTrialsKetamine = nanTrialsKetamine;
        uNanTrialsKetamine(uEyeClosedKetamine) = 0;
        uArtifactsKetamine = artifactsKetamine;
        uArtifactsKetamine(uEyeClosedKetamine | uNanTrialsKetamine) = 0;

    	%trials that do not contains artifacts, nans or blinks 
        allIncludedTrialsSaline = ~((tempFile.lfp.trialInfo.outlierInfo{1}.electrode(electrodeCntr).eyesClosed + ...
                                tempFile.lfp.trialInfo.outlierInfo{1}.electrode(electrodeCntr).nanTrials + ...
                                tempFile.lfp.trialInfo.outlierInfo{1}.electrode(electrodeCntr).artifact)>0);

        allIncludedTrialsKetamine = ~((tempFile.lfp.trialInfo.outlierInfo{2}.electrode(electrodeCntr).eyesClosed + ...
                                tempFile.lfp.trialInfo.outlierInfo{2}.electrode(electrodeCntr).nanTrials + ...
                                tempFile.lfp.trialInfo.outlierInfo{2}.electrode(electrodeCntr).artifact)>0);

        %what are the conditions
        conditionsSaline = tempFile.lfp.trialInfo.conIdent{1};
        conditionsKetamine = tempFile.lfp.trialInfo.conIdent{2};

        %time after injection that a trial was completed
        sessionTimeSaline = tempFile.lfp.trialInfo.realTime{1};
        sessionTimeKetamine = tempFile.lfp.trialInfo.realTime{2};

        subselectionSaline = sessionTimeSaline<=60;
        subselectionKetamine = sessionTimeKetamine<=60;
        allIncludedTrialsSaline = allIncludedTrialsSaline & subselectionSaline;
        allIncludedTrialsKetamine= allIncludedTrialsKetamine & subselectionKetamine;
        useInfo.subject(totalElectrodeCntr,:) = tempFile.subject;
        useInfo.session(totalElectrodeCntr) = fileName;
        useInfo.electrode(totalElectrodeCntr) = electrodeCntr;

        useInfo.rejectedTrials.saline.nrFixationBreakTrials(totalElectrodeCntr,1) = sum(uEyeClosedSaline(subselectionSaline));
        useInfo.rejectedTrials.saline.nrNaNTrials(totalElectrodeCntr,1) = sum(uNanTrialsSaline(subselectionSaline));
        useInfo.rejectedTrials.saline.nrArtifactTrials(totalElectrodeCntr,1) = sum(uArtifactsSaline(subselectionSaline));
        useInfo.rejectedTrials.ketamine.nrFixationBreakTrials(totalElectrodeCntr,1) = sum(uEyeClosedKetamine(subselectionKetamine));
        useInfo.rejectedTrials.ketamine.nrNaNTrials(totalElectrodeCntr,1) = sum(uNanTrialsKetamine(subselectionKetamine));
        useInfo.rejectedTrials.ketamine.nrArtifactTrials(totalElectrodeCntr,1) = sum(uArtifactsKetamine(subselectionKetamine));
        

        for conditionCntr = 1:numel(uCond)

            useInfo.rejectedTrials.saline.nrFixationBreakTrials(totalElectrodeCntr,conditionCntr+1) = sum(uEyeClosedSaline(subselectionSaline & conditionsSaline==uCond(conditionCntr)));
            useInfo.rejectedTrials.saline.nrNaNTrials(totalElectrodeCntr,conditionCntr+1) = sum(uNanTrialsSaline(subselectionSaline & conditionsSaline==uCond(conditionCntr)));
            useInfo.rejectedTrials.saline.nrArtifactTrials(totalElectrodeCntr,conditionCntr+1) = sum(uArtifactsSaline(subselectionSaline & conditionsSaline==uCond(conditionCntr)));
            useInfo.rejectedTrials.ketamine.nrFixationBreakTrials(totalElectrodeCntr,conditionCntr+1) = sum(uEyeClosedKetamine(subselectionKetamine & conditionsKetamine==uCond(conditionCntr)));
            useInfo.rejectedTrials.ketamine.nrNaNTrials(totalElectrodeCntr,conditionCntr+1) = sum(uNanTrialsKetamine(subselectionKetamine & conditionsKetamine==uCond(conditionCntr)));
            useInfo.rejectedTrials.ketamine.nrArtifactTrials(totalElectrodeCntr,conditionCntr+1) = sum(uArtifactsKetamine(subselectionKetamine & conditionsKetamine==uCond(conditionCntr)));


            %include only trials of recordings that have at least 20 trials for that condition
            useInfo.trialSelection.saline{totalElectrodeCntr,conditionCntr+1} = allIncludedTrialsSaline & conditionsSaline==uCond(conditionCntr);
            useInfo.trialSelection.saline{totalElectrodeCntr,conditionCntr+1} = (useInfo.trialSelection.saline{totalElectrodeCntr,conditionCntr+1}) * sum(useInfo.trialSelection.saline{totalElectrodeCntr,conditionCntr+1})>=20;
            useInfo.trialSelection.ketamine{totalElectrodeCntr,conditionCntr+1} = allIncludedTrialsKetamine & conditionsKetamine==uCond(conditionCntr);
            useInfo.trialSelection.ketamine{totalElectrodeCntr,conditionCntr+1} = useInfo.trialSelection.ketamine{totalElectrodeCntr,conditionCntr+1} * sum(useInfo.trialSelection.ketamine{totalElectrodeCntr,conditionCntr+1})>=20;

            useInfo.rejectedTrials.saline.sufficientTrials(totalElectrodeCntr,conditionCntr) = sum(useInfo.trialSelection.saline{totalElectrodeCntr,conditionCntr+1})>=20;
            useInfo.rejectedTrials.ketamine.sufficientTrials(totalElectrodeCntr,conditionCntr) = sum(useInfo.trialSelection.ketamine{totalElectrodeCntr,conditionCntr+1})>=20;


            %calculate SNR by using only first 60 minutes
            tempMeanSal = mean(tempFile.lfp.signal{1}(:,allIncludedTrialsSaline & conditionsSaline==uCond(conditionCntr),electrodeCntr),2,'omitnan');
            baselineStd = std(abs(tempMeanSal(time>=-500 & time<0)));
            onsetMean = mean(abs(tempMeanSal(time>=51 & time<=250)));
            onsetSNR = (onsetMean)/baselineStd;
            useInfo.snr.saline(totalElectrodeCntr,conditionCntr) = onsetSNR;
            averageSignal.saline(:,conditionCntr,totalElectrodeCntr) = tempMeanSal;

            tempMeanKet = mean(tempFile.lfp.signal{2}(:,allIncludedTrialsKetamine & conditionsKetamine==uCond(conditionCntr),electrodeCntr),2,'omitnan');
            baselineStd = std(abs(tempMeanKet(time>=-500 & time<0)));
            onsetMean = mean(abs(tempMeanKet(time>=51 & time<=250)));
            onsetSNR = (onsetMean)/baselineStd;
            useInfo.snr.ketamine(totalElectrodeCntr,conditionCntr) = onsetSNR;
            averageSignal.ketamine(:,conditionCntr,totalElectrodeCntr) = tempMeanKet;

            %use only first 60 minutes to calculate zScored signal
            combinedSignal = mean([tempFile.lfp.signal{1}(:,allIncludedTrialsSaline & conditionsSaline==uCond(conditionCntr),electrodeCntr) tempFile.lfp.signal{2}(:,allIncludedTrialsKetamine & conditionsKetamine==uCond(conditionCntr),electrodeCntr)],2);
            combinedSignalStd = std(combinedSignal(time>=-500 & time<0));
            signalSalZScored = tempMeanSal/combinedSignalStd;
            signalKetZScored = tempMeanKet/combinedSignalStd;

            zScoredSignal.saline(:,conditionCntr,totalElectrodeCntr) = signalSalZScored;
            zScoredSignal.ketamine(:,conditionCntr,totalElectrodeCntr) = signalKetZScored;

        end
    end
end



%% book keeping: create information about rejected electrodes and trials
%how many electrodes were rejected becaue of snr (we only use the 0 Hz
%condition to determine SNR
rejectInfo.snr.saline = useInfo.snr.saline(:,1)<2.5 & useInfo.snr.ketamine(:,1)>=2.5;   %only saline recording is bad
rejectInfo.snr.ketamine = useInfo.snr.ketamine(:,1)<2.5 & useInfo.snr.saline(:,1)>=2.5; %only ketamine recording is bad
rejectInfo.snr.both = useInfo.snr.ketamine(:,1)<2.5 & useInfo.snr.saline(:,1)<2.5; %both, saline and ketamine recording are bad

%how many electrodes that were not excluded so far were rejected because
%either the saline or ketamine recording had too few trials left
rejectInfo.insufficientTrials.saline = (useInfo.rejectedTrials.ketamine.sufficientTrials-useInfo.rejectedTrials.saline.sufficientTrials)==1;
rejectInfo.insufficientTrials.saline(useInfo.snr.saline(:,1)<2.5 | useInfo.snr.ketamine(:,1)<2.5,:) = 0;    %do not count electrodes that have already been excluded
rejectInfo.insufficientTrials.ketamine = (useInfo.rejectedTrials.ketamine.sufficientTrials-useInfo.rejectedTrials.saline.sufficientTrials)==-1;
rejectInfo.insufficientTrials.ketamine(useInfo.snr.saline(:,1)<2.5 | useInfo.snr.ketamine(:,1)<2.5,:) = 0;
rejectInfo.insufficientTrials.both = (useInfo.rejectedTrials.saline.sufficientTrials+ useInfo.rejectedTrials.ketamine.sufficientTrials)==0;
rejectInfo.insufficientTrials.both(useInfo.snr.saline(:,1)<2.5 & useInfo.snr.ketamine(:,1)<2.5,:) = 0;


%combining information about electrodes that had sufficients trials and a
%good SNR provides the  inclusion critereon we will use for all further
%analyses for this study
useInfo.doUse = useInfo.snr.saline(:,1)>=2.5 & useInfo.snr.ketamine(:,1)>=2.5 & ...     %electrodes with an snr of at least 2.5
                useInfo.rejectedTrials.saline.sufficientTrials & useInfo.rejectedTrials.ketamine.sufficientTrials;%electrodes with at least 20 trials are included


%of the included electrodes, what is the number trials that were rejected
    % -because the monkey lost fixation (e.g., due to blink) during a trial
    for conditionCntr = 2:size(useInfo.rejectedTrials.saline.nrFixationBreakTrials,2)
        rejectedTrials.saline.perCondition.blinks.mean(conditionCntr-1) = mean(useInfo.rejectedTrials.saline.nrFixationBreakTrials(useInfo.doUse(:,conditionCntr-1),conditionCntr));
        rejectedTrials.saline.perCondition.blinks.std(conditionCntr-1) = std(useInfo.rejectedTrials.saline.nrFixationBreakTrials(useInfo.doUse(:,conditionCntr-1),conditionCntr));
        rejectedTrials.saline.perCondition.blinks.range(conditionCntr-1,:) = [min(useInfo.rejectedTrials.saline.nrFixationBreakTrials(useInfo.doUse(:,conditionCntr-1),conditionCntr)) max(useInfo.rejectedTrials.saline.nrFixationBreakTrials(useInfo.doUse(:,conditionCntr-1),conditionCntr))];
        rejectedTrials.ketamine.perCondition.blinks.mean(conditionCntr-1) = mean(useInfo.rejectedTrials.ketamine.nrFixationBreakTrials(useInfo.doUse(:,conditionCntr-1),conditionCntr));
        rejectedTrials.ketamine.perCondition.blinks.std(conditionCntr-1) = std(useInfo.rejectedTrials.ketamine.nrFixationBreakTrials(useInfo.doUse(:,conditionCntr-1),conditionCntr));
        rejectedTrials.ketamine.perCondition.blinks.range(conditionCntr-1,:) = [min(useInfo.rejectedTrials.ketamine.nrFixationBreakTrials(useInfo.doUse(:,conditionCntr-1),conditionCntr)) max(useInfo.rejectedTrials.ketamine.nrFixationBreakTrials(useInfo.doUse(:,conditionCntr-1),conditionCntr))];
    end
    rejectedTrials.saline.acrossConditions.blinks.mean = mean(rejectedTrials.saline.perCondition.blinks.mean);
    rejectedTrials.saline.acrossConditions.blinks.std = std(rejectedTrials.saline.perCondition.blinks.mean);
    rejectedTrials.saline.acrossConditions.blinks.range = [min(rejectedTrials.saline.perCondition.blinks.range(:)) max(rejectedTrials.saline.perCondition.blinks.range(:))];
   	rejectedTrials.ketamine.acrossConditions.blinks.mean = mean(rejectedTrials.ketamine.perCondition.blinks.mean);
    rejectedTrials.ketamine.acrossConditions.blinks.std = std(rejectedTrials.ketamine.perCondition.blinks.mean);
    rejectedTrials.ketamine.acrossConditions.blinks.range = [min(rejectedTrials.ketamine.perCondition.blinks.range(:)) max(rejectedTrials.ketamine.perCondition.blinks.range(:))];
    
    % -because a trial had missing data n the recording of local field
    % potentials (this has apparently not happened)
  	for conditionCntr = 2:size(useInfo.rejectedTrials.saline.nrFixationBreakTrials,2)
        rejectedTrials.saline.perCondition.nans.mean(conditionCntr-1) = mean(useInfo.rejectedTrials.saline.nrNaNTrials(useInfo.doUse(:,conditionCntr-1),conditionCntr));
        rejectedTrials.saline.perCondition.nans.std(conditionCntr-1) = std(useInfo.rejectedTrials.saline.nrNaNTrials(useInfo.doUse(:,conditionCntr-1),conditionCntr));
        rejectedTrials.saline.perCondition.nans.range(conditionCntr-1,:) = [min(useInfo.rejectedTrials.saline.nrNaNTrials(useInfo.doUse(:,conditionCntr-1),conditionCntr)) max(useInfo.rejectedTrials.saline.nrNaNTrials(useInfo.doUse(:,conditionCntr-1),conditionCntr))];
        rejectedTrials.ketamine.perCondition.nans.mean(conditionCntr-1) = mean(useInfo.rejectedTrials.ketamine.nrNaNTrials(useInfo.doUse(:,conditionCntr-1),conditionCntr));
        rejectedTrials.ketamine.perCondition.nans.std(conditionCntr-1) = std(useInfo.rejectedTrials.ketamine.nrNaNTrials(useInfo.doUse(:,conditionCntr-1),conditionCntr));
        rejectedTrials.ketamine.perCondition.nans.range(conditionCntr-1,:) = [min(useInfo.rejectedTrials.ketamine.nrNaNTrials(useInfo.doUse(:,conditionCntr-1),conditionCntr)) max(useInfo.rejectedTrials.ketamine.nrNaNTrials(useInfo.doUse(:,conditionCntr-1),conditionCntr))];
    end
    rejectedTrials.saline.acrossConditions.nans.mean = mean(rejectedTrials.saline.perCondition.nans.mean);
    rejectedTrials.saline.acrossConditions.nans.std = std(rejectedTrials.saline.perCondition.nans.mean);
    rejectedTrials.saline.acrossConditions.nans.range = [min(rejectedTrials.saline.perCondition.nans.range(:)) max(rejectedTrials.saline.perCondition.nans.range(:))];
   	rejectedTrials.ketamine.acrossConditions.nans.mean = mean(rejectedTrials.ketamine.perCondition.nans.mean);
    rejectedTrials.ketamine.acrossConditions.nans.std = std(rejectedTrials.ketamine.perCondition.nans.mean);
    rejectedTrials.ketamine.acrossConditions.nans.range = [min(rejectedTrials.ketamine.perCondition.nans.range(:)) max(rejectedTrials.ketamine.perCondition.nans.range(:))];

    % -because a trial had artifacts in the lfp recording (e.g., movement) 
  	for conditionCntr = 2:size(useInfo.rejectedTrials.saline.nrFixationBreakTrials,2)
        rejectedTrials.saline.perCondition.artifacts.mean(conditionCntr-1) = mean(useInfo.rejectedTrials.saline.nrArtifactTrials(useInfo.doUse(:,conditionCntr-1),conditionCntr));
        rejectedTrials.saline.perCondition.artifacts.std(conditionCntr-1) = std(useInfo.rejectedTrials.saline.nrArtifactTrials(useInfo.doUse(:,conditionCntr-1),conditionCntr));
        rejectedTrials.saline.perCondition.artifacts.range(conditionCntr-1,:) = [min(useInfo.rejectedTrials.saline.nrArtifactTrials(useInfo.doUse(:,conditionCntr-1),conditionCntr)) max(useInfo.rejectedTrials.saline.nrArtifactTrials(useInfo.doUse(:,conditionCntr-1),conditionCntr))];
        rejectedTrials.ketamine.perCondition.artifacts.mean(conditionCntr-1) = mean(useInfo.rejectedTrials.ketamine.nrArtifactTrials(useInfo.doUse(:,conditionCntr-1),conditionCntr));
        rejectedTrials.ketamine.perCondition.artifacts.std(conditionCntr-1) = std(useInfo.rejectedTrials.ketamine.nrArtifactTrials(useInfo.doUse(:,conditionCntr-1),conditionCntr));
        rejectedTrials.ketamine.perCondition.artifacts.range(conditionCntr-1,:) = [min(useInfo.rejectedTrials.ketamine.nrArtifactTrials(useInfo.doUse(:,conditionCntr-1),conditionCntr)) max(useInfo.rejectedTrials.ketamine.nrArtifactTrials(useInfo.doUse(:,conditionCntr-1),conditionCntr))];
    end
    rejectedTrials.saline.acrossConditions.artifacts.mean = mean(rejectedTrials.saline.perCondition.artifacts.mean);
    rejectedTrials.saline.acrossConditions.artifacts.std = std(rejectedTrials.saline.perCondition.artifacts.mean);
    rejectedTrials.saline.acrossConditions.artifacts.range = [min(rejectedTrials.saline.perCondition.artifacts.range(:)) max(rejectedTrials.saline.perCondition.artifacts.range(:))];
   	rejectedTrials.ketamine.acrossConditions.artifacts.mean = mean(rejectedTrials.ketamine.perCondition.artifacts.mean);
    rejectedTrials.ketamine.acrossConditions.artifacts.std = std(rejectedTrials.ketamine.perCondition.artifacts.mean);
    rejectedTrials.ketamine.acrossConditions.artifacts.range = [min(rejectedTrials.ketamine.perCondition.artifacts.range(:)) max(rejectedTrials.ketamine.perCondition.artifacts.range(:))];
    
%how many trials did each electrode we included supply
for conditionCntr =  2:size(useInfo.trialSelection.saline,2)
    for electrodeCntr = 1:size(useInfo.trialSelection.saline,1)
        nrTrialsIncludedSaline(electrodeCntr,conditionCntr-1) = sum(useInfo.trialSelection.saline{electrodeCntr,conditionCntr});
        nrTrialsIncludedKetamine(electrodeCntr,conditionCntr-1) = sum(useInfo.trialSelection.ketamine{electrodeCntr,conditionCntr});
    end
end
   
nrTrialsIncludedSaline(~useInfo.doUse) = NaN;
nrTrialsIncludedKetamine(~useInfo.doUse) = NaN;
nrTrialsIncludedCombined = [nrTrialsIncludedSaline; nrTrialsIncludedKetamine];

includedTrials.saline.perCondition.mean = mean(nrTrialsIncludedSaline,'omitnan');
includedTrials.ketamine.perCondition.mean = mean(nrTrialsIncludedKetamine,'omitnan');
includedTrials.saline.perCondition.std = std(nrTrialsIncludedSaline,'omitnan');
includedTrials.ketamine.perCondition.std = std(nrTrialsIncludedKetamine,'omitnan');
includedTrials.saline.perCondition.range = [min(nrTrialsIncludedSaline)', max(nrTrialsIncludedSaline)'];
includedTrials.ketamine.perCondition.range = [min(nrTrialsIncludedKetamine)', max(nrTrialsIncludedKetamine)'];
includedTrials.combined.perCondition.mean = mean(nrTrialsIncludedCombined,'omitnan');
includedTrials.combined.perCondition.std = std(nrTrialsIncludedCombined,'omitnan');
includedTrials.combined.perCondition.range = [min(nrTrialsIncludedCombined)', max(nrTrialsIncludedCombined)'];



%save the information we want to report in a readable format
    %electrodes removed for poor snr
    electrodesRemovedForLowSNR.saline = sum(rejectInfo.snr.saline); %because only saline snr was below 2.5
    electrodesRemovedForLowSNR.ketamine = sum(rejectInfo.snr.ketamine); %because only saline snr was below 2.5
    electrodesRemovedForLowSNR.both = sum(rejectInfo.snr.both); %because only saline snr was below 2.5
    electrodesRemovedForLowSNR.total = electrodesRemovedForLowSNR.saline + electrodesRemovedForLowSNR.ketamine + electrodesRemovedForLowSNR.both ; %because only saline snr was below 2.5
    
    %electrodes removed for too low nr trials
    electrodesRemovedForTrialCount.saline = sum(rejectInfo.insufficientTrials.saline);
    electrodesRemovedForTrialCount.ketamine = sum(rejectInfo.insufficientTrials.ketamine);
    electrodesRemovedForTrialCount.both = sum(rejectInfo.insufficientTrials.both);
    electrodesRemovedForTrialCount.total = electrodesRemovedForTrialCount.saline+electrodesRemovedForTrialCount.ketamine+electrodesRemovedForTrialCount.both;
    electrodesRemovedForTrialCount.unique = sum((sum(rejectInfo.insufficientTrials.saline,2) + sum(rejectInfo.insufficientTrials.ketamine,2))>0);
    
    %how many electrodes were ultimately included
    electrodesIncluded.totalIncluded = sum(useInfo.doUse);
    electrodesIncluded.totalExcluded = sum(~useInfo.doUse);
    electrodesIncluded.total = electrodesIncluded.totalIncluded+electrodesIncluded.totalExcluded ;
    
%how many trials did we reject (on average) and for what reason
    conditionNames = {'0Hz'; '5Hz'; '10Hz'; '20Hz'; '40Hz'; 'all'};
    indexCntr = 0;
    for conditionCntr = 1:6
        indexCntr = indexCntr+1;
        if conditionCntr < 6
        reportRejectedTrials.blinks_saline_mean(indexCntr,:) = rejectedTrials.saline.perCondition.blinks.mean(conditionCntr);
        reportRejectedTrials.blinks_saline_std(indexCntr,:) = rejectedTrials.saline.perCondition.blinks.std(conditionCntr);
        reportRejectedTrials.blinks_saline_range(indexCntr,:) = rejectedTrials.saline.perCondition.blinks.range(conditionCntr,:);
        reportRejectedTrials.blinks_ketamine_mean(indexCntr,:) = rejectedTrials.ketamine.perCondition.blinks.mean(conditionCntr);
        reportRejectedTrials.blinks_ketamine_std(indexCntr,:) = rejectedTrials.ketamine.perCondition.blinks.std(conditionCntr);
        reportRejectedTrials.blinks_ketamine_range(indexCntr,:) = rejectedTrials.ketamine.perCondition.blinks.range(conditionCntr,:);
        
        reportRejectedTrials.nans_saline_mean(indexCntr,:) = rejectedTrials.saline.perCondition.nans.mean(conditionCntr);
        reportRejectedTrials.nans_saline_std(indexCntr,:) = rejectedTrials.saline.perCondition.nans.std(conditionCntr);
        reportRejectedTrials.nans_saline_range(indexCntr,:) = rejectedTrials.saline.perCondition.nans.range(conditionCntr,:);
        reportRejectedTrials.nans_ketamine_mean(indexCntr,:) = rejectedTrials.ketamine.perCondition.nans.mean(conditionCntr);
        reportRejectedTrials.nans_ketamine_std(indexCntr,:) = rejectedTrials.ketamine.perCondition.nans.std(conditionCntr);
        reportRejectedTrials.nans_ketamine_range(indexCntr,:) = rejectedTrials.ketamine.perCondition.nans.range(conditionCntr,:);
        
       	reportRejectedTrials.artifacts_saline_mean(indexCntr,:) = rejectedTrials.saline.perCondition.artifacts.mean(conditionCntr);
        reportRejectedTrials.artifacts_saline_std(indexCntr,:) = rejectedTrials.saline.perCondition.artifacts.std(conditionCntr);
        reportRejectedTrials.artifacts_saline_range(indexCntr,:) = rejectedTrials.saline.perCondition.artifacts.range(conditionCntr,:);
        reportRejectedTrials.artifacts_ketamine_mean(indexCntr,:) = rejectedTrials.ketamine.perCondition.artifacts.mean(conditionCntr);
        reportRejectedTrials.artifacts_ketamine_std(indexCntr,:) = rejectedTrials.ketamine.perCondition.artifacts.std(conditionCntr);
        reportRejectedTrials.artifacts_ketamine_range(indexCntr,:) = rejectedTrials.ketamine.perCondition.artifacts.range(conditionCntr,:);
        else
            reportRejectedTrials.blinks_saline_mean(indexCntr,:) = rejectedTrials.saline.acrossConditions.blinks.mean;
            reportRejectedTrials.blinks_saline_std(indexCntr,:) = rejectedTrials.saline.acrossConditions.blinks.std;
            reportRejectedTrials.blinks_saline_range(indexCntr,:) = rejectedTrials.saline.acrossConditions.blinks.range;
            reportRejectedTrials.blinks_ketamine_mean(indexCntr,:) = rejectedTrials.ketamine.acrossConditions.blinks.mean;
            reportRejectedTrials.blinks_ketamine_std(indexCntr,:) = rejectedTrials.ketamine.acrossConditions.blinks.std;
            reportRejectedTrials.blinks_ketamine_range(indexCntr,:) = rejectedTrials.ketamine.acrossConditions.blinks.range;    
            
        	reportRejectedTrials.nans_saline_mean(indexCntr,:) = rejectedTrials.saline.acrossConditions.nans.mean;
            reportRejectedTrials.nans_saline_std(indexCntr,:) = rejectedTrials.saline.acrossConditions.nans.std;
            reportRejectedTrials.nans_saline_range(indexCntr,:) = rejectedTrials.saline.acrossConditions.nans.range;
            reportRejectedTrials.nans_ketamine_mean(indexCntr,:) = rejectedTrials.ketamine.acrossConditions.nans.mean;
            reportRejectedTrials.nans_ketamine_std(indexCntr,:) = rejectedTrials.ketamine.acrossConditions.nans.std;
            reportRejectedTrials.nans_ketamine_range(indexCntr,:) = rejectedTrials.ketamine.acrossConditions.nans.range;    
            
           	reportRejectedTrials.artifacts_saline_mean(indexCntr,:) = rejectedTrials.saline.acrossConditions.artifacts.mean;
            reportRejectedTrials.artifacts_saline_std(indexCntr,:) = rejectedTrials.saline.acrossConditions.artifacts.std;
            reportRejectedTrials.artifacts_saline_range(indexCntr,:) = rejectedTrials.saline.acrossConditions.artifacts.range;
            reportRejectedTrials.artifacts_ketamine_mean(indexCntr,:) = rejectedTrials.ketamine.acrossConditions.artifacts.mean;
            reportRejectedTrials.artifacts_ketamine_std(indexCntr,:) = rejectedTrials.ketamine.acrossConditions.artifacts.std;
            reportRejectedTrials.artifacts_ketamine_range(indexCntr,:) = rejectedTrials.ketamine.acrossConditions.artifacts.range; 
        end
    end
   
%how many trials were ultimately included
    for conditionCntr = 1:5
        reportIncludedTrials.saline_mean(conditionCntr,:) = includedTrials.saline.perCondition.mean(conditionCntr);
        reportIncludedTrials.saline_std(conditionCntr,:) = includedTrials.saline.perCondition.std(conditionCntr);
        reportIncludedTrials.saline_range(conditionCntr,:) = includedTrials.saline.perCondition.range(conditionCntr,:);
       	reportIncludedTrials.ketamine_mean(conditionCntr,:) = includedTrials.ketamine.perCondition.mean(conditionCntr);
        reportIncludedTrials.ketamine_std(conditionCntr,:) = includedTrials.ketamine.perCondition.std(conditionCntr);
        reportIncludedTrials.ketamine_range(conditionCntr,:) = includedTrials.ketamine.perCondition.range(conditionCntr,:);
      	reportIncludedTrials.combined_mean(conditionCntr,:) = includedTrials.combined.perCondition.mean(conditionCntr);
        reportIncludedTrials.combined_std(conditionCntr,:) = includedTrials.combined.perCondition.std(conditionCntr);
        reportIncludedTrials.combined_range(conditionCntr,:) = includedTrials.combined.perCondition.range(conditionCntr,:);
    end
    
    

%save information about data quality and zscored signal
save([targetFolder 'useInfo'], 'useInfo','-v7.3');
save([targetFolder 'averageSignal'], 'averageSignal','-v7.3');
save([targetFolder 'zScoredSignal'], 'zScoredSignal','-v7.3');

end