function[overlapLit, overlapModRef_Exp2Rxn, overlapMapExp2Rxn, model_out_bm, simpleCheck, fluxInconsistentRxns, model_imBalancedChargeRxns] = modelStats(model);
tic;
changeCobraSolver('gurobi','all'); clear ans;
%% load ref, lit models
model_ref = getDistributedModel('Recon3DModel_301.mat'); %model
astroModel = xls2model('AstroModel_New.xlsx'); 
%% mapExpressionToReactions (Transcriptome)
data = importdata('GSE73721_HMA_CTX_v10.mat');
gene = (data.Entrez_ID); gene = num2cell(gene); gene = cellfun(@num2str,gene,'uni',0);
expression.gene = strcat(gene,'.1');
expression.value = max(data.matrix,[],2); %for the row maximums && max(A(:)); %for the matrix maximum
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