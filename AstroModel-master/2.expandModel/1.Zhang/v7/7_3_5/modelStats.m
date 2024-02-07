function[matSize, IntRxns, ExRxns, BiDirExc, OnlyReleaseExc, OnlyUptakeExc] = stats(model)
%%
tic;
changeCobraSolver('gurobi','all');
model = changeObjective(model, 'biomass_maintenance', 0);
%% size of the stoichiometric matrix
matSize = size(model.S);
%%
%% Internal reaction
IntRxns = unique(model.rxns(findExcRxns(model)~=1));
%% Exchange reactions
modelexchanges1 = strmatch('Ex_',model.rxns);
modelexchanges4 = strmatch('EX_',model.rxns);
modelexchanges2 = strmatch('DM_',model.rxns);
modelexchanges3 = strmatch('sink_',model.rxns);
selExc = (find( full((sum(abs(model.S)==1,1) ==1) & (sum(model.S~=0) == 1))))';
modelexchanges = unique([modelexchanges1;modelexchanges2;modelexchanges3;modelexchanges4;selExc]);
ExRxns = unique(model.rxns(modelexchanges));
%% Blood exchanges
BloodRxns = findRxnFromCompartment(model,'[e]');
BloodRxns = BloodRxns(:,[1]);
BloodExcRxns = unique(intersect(ExRxns,BloodRxns));
%% Synapse exchanges
SynapseRxns = findRxnFromCompartment(model,'[s]');
SynapseRxns = SynapseRxns(:,[1]);
SynapseExcRxns = unique(intersect(ExRxns,SynapseRxns));
%% Both uptake and release Exchanges
BiDirExc = model.rxns(model.lb<0 & model.ub>0);
BiDirExc = unique(intersect(ExRxns,BiDirExc));
OpenUptakesBlood = unique(intersect(BloodExcRxns,BiDirExc));
OpenUptakesSynapse = unique(intersect(SynapseExcRxns,BiDirExc));
OpenUptakesSink = unique(intersect(model.rxns(modelexchanges3),BiDirExc));

%% Only released exchanges
OnlyReleaseExc = model.rxns(model.lb==0 & model.ub>0);
OnlyReleaseExc = unique(intersect(ExRxns,OnlyReleaseExc));

%% Only uptake exchanges
OnlyUptakeExc = model.rxns(model.lb<0 & model.ub<=0);
OnlyUptakeExc = unique(intersect(ExRxns,OnlyUptakeExc));

%% inactive exchanges
closedRxns = model.rxns(model.lb==0 & model.ub==0);

%% active exchanges after FVA:
% ans = findBlockedReaction(model);
% ans_bl = ans';
% ans_bl(strmatch('EX_',ans_bl));
% ans_bl(strmatch('DM_',ans_bl));
% ans_bl(strmatch('sink_',ans_bl));

%%
toc;
end