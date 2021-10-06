function[iMAT_model_T, iMAT_model_TP, GIMME_model_T, GIMME_model_TP, MBA_model_T, MBA_model_TP] = mat2mod(filename);
tic;
%%
changeCobraSolver('gurobi','all'); clear ans;
%% load reference model
model = getDistributedModel('Recon3DModel_301.mat'); %model
model = changeObjective(model, 'biomass_maintenance', 1);
%% load expression data
data = importdata(filename);

%% NORMALIZE (LOG-TRANSFORM & QUANTILENORM) DATA MATRIX
data.matrix = quantilenorm(log10(data.matrix+1));

%% (Model_Transcriptome)
%% mapExpressionToReactions (Transcriptome)
gene = (data.Entrez_ID); gene = num2cell(gene); gene = cellfun(@num2str,gene,'uni',0);
expression.gene = strcat(gene,'.1');
expression.value = max(data.matrix,[],2); %for the row maximums && max(A(:)); %for the matrix maximum
[expressionRxns, parsedGPR] = mapExpressionToReactions(model, expression);
%% createTissueSpecificModel_iMAT (Transcriptome)
threshold_lb=prctile(expression.value,25); threshold_ub=prctile(expression.value,50); solver='iMAT';
options=struct('expressionRxns',expressionRxns,'threshold_lb',threshold_lb,'threshold_ub',threshold_ub,'solver',solver);
iMAT_model_T = createTissueSpecificModel(model, options);
%% createTissueSpecificModel_GIMME: (Transcriptome)
threshold = prctile(expression.value,50); solver = 'GIMME';
options_GIMME=struct('expressionRxns',expressionRxns,'threshold',threshold,'solver',solver);
options = options_GIMME;
GIMME_model_T = createTissueSpecificModel(model, options);
%% createTissueSpecificModel_MBA: (Transcriptome)
high_set=model.rxns(expressionRxns>=prctile(expression.value,50));
medium_set=model.rxns(expressionRxns<prctile(expression.value,50) & expressionRxns>prctile(expression.value,25));
solver = 'MBA';
options=struct('high_set',{high_set},'medium_set',{medium_set},'solver',{solver});
changeCobraSolver('gurobi','all');
MBA_model_T = createTissueSpecificModel(model, options);
%% createTissueSpecificModel_FastCore: (Transcriptome)
% core = find(ismember(model.rxns,model.rxns(expressionRxns>=prctile(expression.value,50)))==1);
% solver = 'fastCore';
% epsilon = 1e-6;
% options=struct('core',{core},'solver',{solver},'epsilon',{epsilon});
% FastCore_model_T = createTissueSpecificModel(model, options);

%% (Model_Transcriptome_HPA)
%% mapExpressionToReactions (Transcriptome + HPA)
HPA_Level_NA = strfind(data.HPA_Level,'NA');
HPA_Level_YES = find(cellfun(@isempty,HPA_Level_NA)); %YES implies the protein expression can be low/med/high
gene = (data.Entrez_ID(HPA_Level_YES)); gene = num2cell(gene); gene = cellfun(@num2str,gene,'uni',0);
expression.gene = strcat(gene,'.1');
expression.value = max(data.matrix(HPA_Level_YES,:),[],2);
[expressionRxns, parsedGPR] = mapExpressionToReactions(model, expression);
%% createTissueSpecificModel_iMAT (Transcriptome + HPA)
threshold_lb=prctile(expression.value,25); threshold_ub=prctile(expression.value,50); solver='iMAT';
options=struct('expressionRxns',expressionRxns,'threshold_lb',threshold_lb,'threshold_ub',threshold_ub,'solver',solver);
iMAT_model_TP = createTissueSpecificModel(model, options);
%% createTissueSpecificModel_GIMME: (Transcriptome + HPA)
threshold = prctile(expression.value,50); solver = 'GIMME';
options_GIMME=struct('expressionRxns',expressionRxns,'threshold',threshold,'solver',solver);
options = options_GIMME;
GIMME_model_TP = createTissueSpecificModel(model, options);
%% createTissueSpecificModel_MBA: (Transcriptome + HPA)
high_set=model.rxns(expressionRxns>=prctile(expression.value,50));
medium_set=model.rxns(expressionRxns<prctile(expression.value,50) & expressionRxns>prctile(expression.value,25));
solver = 'MBA';
options=struct('high_set',{high_set},'medium_set',{medium_set},'solver',{solver});
changeCobraSolver('gurobi','all');
MBA_model_TP = createTissueSpecificModel(model, options);
%% createTissueSpecificModel_FastCore: (Transcriptome + HPA)
% core = find(ismember(model.rxns,model.rxns(expressionRxns>=prctile(expression.value,50)))==1);
% solver = 'fastCore';
% epsilon = 1e-6;
% options=struct('core',{core},'solver',{solver},'epsilon',{epsilon});
% FastCore_model_TP = createTissueSpecificModel(model, options);

%%
toc;
end