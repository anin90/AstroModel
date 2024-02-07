function[overlapLit, overlapModRef_Exp2Rxn, overlapMapExp2Rxn, model_out_bm, simpleCheck, fluxInconsistentRxns, model_imBalancedChargeRxns] = modelStats(model);
tic;
changeCobraSolver('gurobi','all'); clear ans;
%% load ref, lit models
model_ref = getDistributedModel('Recon3DModel_301.mat'); %model
astroModel = xls2model('AstroModel_New.xlsx'); 
%% mapExpressionToReactions (Transcriptome + HPA)
data = importdata('GSE73721_HMA_CTX_v9.mat');
%% NORMALIZE (LOG-TRANSFORM & QUANTILENORM) DATA MATRIX
data.matrix = quantilenorm(log10(data.matrix+1));
%%
HPA_Level_NA = strfind(data.HPA_Level,'NA');
HPA_Level_YES = find(cellfun(@isempty,HPA_Level_NA)); %YES implies the protein expression can be low/med/high
gene = (data.Entrez_ID(HPA_Level_YES)); gene = num2cell(gene); gene = cellfun(@num2str,gene,'uni',0);
expression.gene = strcat(gene,'.1');
expression.value = max(data.matrix(HPA_Level_YES,:),[],2);
[expressionRxns, parsedGPR] = mapExpressionToReactions(model_ref, expression);
%% stats
overlapLit = length(intersect(lower(model.rxns),lower(astroModel.rxns)));
overlapModRef_Exp2Rxn = lower(model_ref.rxns(expressionRxns~=-1));
overlapMapExp2Rxn = length(intersect(lower(model.rxns),lower(model_ref.rxns(expressionRxns~=-1))));
model_out_bm = find(strcmp(model.rxns,'biomass_maintenance'));
simpleCheck = verifyModel(model,'simpleCheck',true);
fluxConsistency = verifyModel(model,'fluxConsistency',true); fluxInconsistentRxns = length(model.rxns(fluxConsistency.fluxConsistency.consistentReactionBool~=1));
chargeBalance = verifyModel(model,'chargeBalance',true); model_imBalancedChargeRxns = model.rxns(chargeBalance.chargeBalance.imBalanceCharge~=0);

%%
toc;
end