function varargout = doStats(inputTable,modelNotation,modelName)
%use linear mixed effects model 

   
    lmeModel = fitlme(inputTable,modelNotation);

    residualsModel = residuals(lmeModel);
    outliers =isoutlier(residualsModel);
    
    outlierIdx = find(outliers);
    salineOutliers = outlierIdx(outlierIdx<=numel(outliers)/2);
    ketamineOutliers = outlierIdx(outlierIdx>numel(outliers)/2);
    
    outliers = unique([salineOutliers' (ketamineOutliers-numel(outliers)/2)' ketamineOutliers' (salineOutliers+numel(outliers)/2)']);
    
    inputTable(outliers,:) = [];

    lmeModel = fitlme(inputTable,modelNotation);
    [robustnessTest, ~, ~] = lm.bootstrap(lmeModel,'subjectVariable','observationNr','mode','FITPLUSNOISE','nrHeteroBins',3,'graph',false,'nrMonteCarlo',1000);
    robustnessTest = robustnessTest.significanceMatch;

    effectModel = lmeModel.Coefficients;
    anovaModel = lmeModel.anova;
    
   	columNames = fieldnames(effectModel);
    columNames([1 end-2:end]) = [];
    columNames = [columNames; {'CI';'robust'}];
    if isempty(modelName)
        modelName = 'lmeModel';
        modelNameAnova = 'anovaModel';
    else
        modelNameAnova = ['anova_' modelName];
    end
    termNames = effectModel.Name;
    entryCntr = 0;
    for termCntr = 1:length(termNames)
        entryCntr = entryCntr+1;
        lmeTable.analysisName{entryCntr,1} = modelName;
        lmeTable.term{entryCntr,1} = termNames{termCntr};
        for columCntr = 1:length(columNames)
            if strcmpi(columNames{columCntr},'ci')
                lmeTable.(columNames{columCntr})(entryCntr,:) = [effectModel.Lower(termCntr) effectModel.Upper(termCntr)];
            elseif strcmpi(columNames{columCntr},'robust')
            lmeTable.(columNames{columCntr})(entryCntr,1) = robustnessTest(termCntr);
            else
                lmeTable.(columNames{columCntr})(entryCntr,1) = effectModel.(columNames{columCntr})(termCntr);
            end
        end
    end
    
    varargout{1} = struct2table(lmeTable);
    
    columNames = fieldnames(anovaModel);
    columNames([1 end]) = [];
    termNames = anovaModel.Term;
    entryCntr = 0;
    for termCntr = 1:length(termNames)
        entryCntr = entryCntr+1;
        anovaTable.analysisName{entryCntr,1} = modelNameAnova;
        anovaTable.term{entryCntr,1} = termNames{termCntr};
        for columCntr = 1:length(columNames)
        	anovaTable.(columNames{columCntr})(entryCntr,1) = anovaModel.(columNames{columCntr})(termCntr);
        end
    end
    varargout{2} = struct2table(anovaTable);
    
    interceptText = lm.disp(lmeModel,{'(Intercept)'});
    interceptText = interceptText(2,3:end-1);
    interceptText = strrep(interceptText,newline,' ');
    modelText = lm.disp(lmeModel);
    modelText = modelText(2:end,3:end-1);
    for lineCntr = 1:size(modelText,1)
        modelText(lineCntr,:) = strrep(modelText(lineCntr,:),newline,' ');
    end
    

    [howLong, whichOne] =max([size(interceptText,2) size(modelText,2)]);
    if whichOne==1
        modelText(:,size(modelText,2)+1:howLong) = repmat(' ',1,numel(size(modelText,2)+1:howLong));
    else
        interceptText(:,size(interceptText,2)+1:howLong) = repmat(' ',1,numel(size(interceptText,2)+1:howLong));
    end
    
    varargout{3} = cat(1,interceptText,modelText);
    
    if nargout>3
        conditionNames = ['5Hz  -- '; '10Hz -- '; '20Hz -- ';'40Hz -- '];
        tempContrastText = cell(numel(unique(lmeModel.Variables.condition)),1);
        for conditionCntr = 1:numel(unique(lmeModel.Variables.condition))
            [p,stat,df,delta,CI,str,~,~] = lm.posthoc(lmeModel,{'condition',num2str(conditionCntr+1),'drug','2'},{'condition',num2str(conditionCntr+1),'drug','1'});
            tempContrastTable.condition{conditionCntr,1} = conditionNames(conditionCntr,:);
            tempContrastTable.fValue(conditionCntr,1) = stat;
            tempContrastTable.df(conditionCntr,1) = df;
            tempContrastTable.delta(conditionCntr,1) = delta;
            tempContrastTable.CI(conditionCntr,:) = CI;
            tempContrastTable.pValue(conditionCntr,:) = p;

            tempContrastText{conditionCntr} = str;
        end
        varargout{4} = struct2table(tempContrastTable);
    end
        
    if nargout>4
        howLong =  max(max(cell2mat(cellfun(@size,tempContrastText,'UniformOutput', 0)),[],2));
        contrastText = repmat(' ',length(tempContrastText),howLong);
        for contrastCntr = 1:length(tempContrastText)
            contrastText(contrastCntr,1:size(tempContrastText{contrastCntr},2)) = tempContrastText{contrastCntr};
        end
        varargout{5} = [conditionNames contrastText];
    end
    
end