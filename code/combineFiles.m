function combineFiles
% This script matches lfp recordings and eyelink edf files for each
% recording day and recording pair (respective saline and ketamine pair)
%
% Companion code for:
%
% N-methyl d-aspartate receptor hypofunction reduces steady state visual
% evoked potentials (2023)
% Alexander Schielke & Bart Krekelberg
% Center for Molecular and Behavioral Neuroscience
% Rutgers University - Newark 

%where are data located
sourceFolder = strrep(pwd,'code','data\');
lfpFolder = [sourceFolder 'ssLFP\'];
eyeFolder = [sourceFolder 'ssEye\'];
lfpFiles = dir(lfpFolder);
lfpFiles = {lfpFiles.name};
lfpFiles(1:2) = [];
eyeFiles = dir(eyeFolder);
eyeFiles = {eyeFiles.name};
eyeFiles(1:2) = [];


%identify all files with preprocessed local field potentials
for fileCntr = 1:length(lfpFiles)
    
    tempLfpFile = load([lfpFolder lfpFiles{fileCntr}]);
    tempLfpFile = tempLfpFile.ssLFP;
    
  	%check that the files belong to each other
   	[~, lfpFileNameSalineTemp, ~] = fileparts(tempLfpFile.partnerFile.fullSalinePath);
    [~, lfpFileNameKetamineTemp, ~] = fileparts(tempLfpFile.partnerFile.fullKetaminePath);
    
    lfpFileNames.saline(fileCntr,:) = lfpFileNameSalineTemp;
    lfpFileNames.ketamine(fileCntr,:) = lfpFileNameKetamineTemp;
end


%identify all files with preprocessed eye and pupil information
for fileCntr = 1:length(eyeFiles)
    
    tempEyeFile = load([eyeFolder eyeFiles{fileCntr}]);
    tempEyeFile = tempEyeFile.ssEye;
    
  	%check that the files belong to each other
    [~, eyeFileNameSalineTemp, ~] = fileparts(tempEyeFile.partnerFile.fullSalinePath);
    [~, eyeFileNameKetamineTemp, ~] = fileparts(tempEyeFile.partnerFile.fullKetaminePath);
    edfFileNameSalineTemp = tempEyeFile.partnerFile.edfPairFileName{1};
    [~, edfFileNameKetamineTemp, ~] = fileparts(tempEyeFile.partnerFile.edfPairFileName{2});
    
    eyeFileNames.saline(fileCntr,:) = eyeFileNameSalineTemp;
    eyeFileNames.ketamine(fileCntr,:) = eyeFileNameKetamineTemp;
    edfFileNames.saline(fileCntr,:) = edfFileNameSalineTemp;
    edfFileNames.ketamine(fileCntr,:) = edfFileNameKetamineTemp; 
end

%match files
eyeFileIdx = nan(size(lfpFileNames.saline,1),1);
edfFileMatch = nan(size(lfpFileNames.saline,1),2);
for fileCntr = 1:size(lfpFileNames.saline,1)
    try
        eyeFileIdx(fileCntr) = find(sum((lfpFileNames.saline(fileCntr,:))==(eyeFileNames.saline),2)==size(eyeFileNames.saline,2));
        edfFileMatch(fileCntr,1) = strcmpi(lfpFileNames.saline(fileCntr,:),edfFileNames.saline(eyeFileIdx(fileCntr),:));
        edfFileMatch(fileCntr,2) = strcmpi(lfpFileNames.ketamine(fileCntr,:),edfFileNames.ketamine(eyeFileIdx(fileCntr),:));
    catch
        eyeFileIdx(fileCntr) = NaN;
        edfFileMatch(fileCntr,1) = 0;
        edfFileMatch(fileCntr,2) = 0;
    end
end


%combineFiles
for fileCntr = 1:length(lfpFiles)
    clearvars data
    
    tempLfpFile = load([lfpFolder lfpFiles{fileCntr}]);
    tempLfpFile = tempLfpFile.ssLFP;
    
    %create files
    data.subject = tempLfpFile.partnerFile.salineFile(1:2);
    data.date = tempLfpFile.partnerFile.partnerDayDir;
    data.startDelay = tempLfpFile.partnerFile.startDelay;
    
    data.lfp.signal = tempLfpFile.signal;
    data.lfp.trialInfo =  tempLfpFile.trialInfo;
    
    %is there a matching edf file
    if ~isnan(eyeFileIdx(fileCntr)) && sum(edfFileMatch(fileCntr,:))==2
        tempEyeFile = load([eyeFolder eyeFiles{eyeFileIdx(fileCntr)}]);
        tempEyeFile = tempEyeFile.ssEye;
        
        if isequal(size(tempEyeFile.signal{1},2),size(tempLfpFile.signal{1},2)) && isequal(size(tempEyeFile.signal{2},2),size(tempLfpFile.signal{2},2))
            data.eye.signal = tempEyeFile.signal;
            data.eye.trialInfo = tempEyeFile.trialInfo;
        else
            error('something is wrong')
        end
    else
        data.eye.signal = {NaN;NaN};
        data.eye.trialInfo = [];
    end
        
    save([sourceFolder 'combined\file' num2str(fileCntr)], 'data','-v7.3');
end

end