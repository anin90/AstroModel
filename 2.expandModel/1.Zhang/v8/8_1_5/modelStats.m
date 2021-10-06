function[total_rxns, Total_ExRxns, EX_Rxns, DM_Rxns, Sink_Rxns, blockedRxns, blockedEX, blockedDM, blockedSink] = modelStats(model)
%%
tic;
changeCobraSolver('gurobi','all');
model = changeObjective(model, 'biomass_maintenance', 0);
%% number of reactions within system
total_rxns = length(model.rxns);
%% number of EX_, DM_ & sink_ rxns;
modelexchanges1 = strmatch('Ex_',model.rxns);
modelexchanges4 = strmatch('EX_',model.rxns);
modelexchanges2 = strmatch('DM_',model.rxns);
modelexchanges3 = strmatch('sink_',model.rxns);
selExc = (find( full((sum(abs(model.S)==1,1) ==1) & (sum(model.S~=0) == 1))))';
modelexchanges = unique([modelexchanges1;modelexchanges2;modelexchanges3;modelexchanges4;selExc]);
Total_ExRxns = unique(model.rxns(modelexchanges));
%% input: EX, DM, sink rxns
EX_Rxns = Total_ExRxns(strmatch('EX_',Total_ExRxns));
DM_Rxns = Total_ExRxns(strmatch('DM_',Total_ExRxns));
Sink_Rxns = Total_ExRxns(strmatch('sink_',Total_ExRxns));
%% output: find all blocked rxns
blockedRxns = findBlockedReaction(model)';
blockedEX = blockedRxns(strmatch('EX_',blockedRxns));
blockedDM = blockedRxns(strmatch('DM_',blockedRxns));
blockedSink = blockedRxns(strmatch('sink_',blockedRxns));
end
