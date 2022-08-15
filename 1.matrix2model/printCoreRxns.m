function[coreRxns] = printCoreRxns(expMat)
%%
tStart = tic;

%% add path to dependencies:
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/')

%% Load data
    modelRef = getDistributedModel('Recon3DModel_301.mat');
    data = importdata(expMat);

%% mapExpressionToReactions (Transcriptome + HPA) (abs/norm_t1/norm_t2)

    HPA_Level_NA = strfind(data.HPA_Level,'NA');
    HPA_Level_YES = find(cellfun(@isempty,HPA_Level_NA)); %YES implies the protein expression can be low/med/high
    gene = (data.Entrez_ID(HPA_Level_YES)); gene = num2cell(gene); gene = cellfun(@num2str,gene,'uni',0);
    expression.gene = strcat(gene,'.1');
    expression.value = max(data.matrix(HPA_Level_YES,:),[],2);
    [expressionRxns, ~] = mapExpressionToReactions(modelRef, expression);
    
%% Print Core Rxns
    coreRxns = modelRef.rxns(expressionRxns~=-1);
    
%% 
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));

%%
end