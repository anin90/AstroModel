function[T,rxnAbsent] = Test4CuratedRxns(model, filename);
tic;

%% LOAD WT MODEL
changeCobraSolver('gurobi','all'); clear ans;
model = changeObjective(model, 'biomass_maintenance', 0);

%% IMPORT RXNS
rxnList = importdata(filename);

%% INDEX of RXNS
[tf,loc] = ismember(model.rxns, rxnList); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);

%% FIND RXNS PRESENT & ABSENT IN MODEL
rxnPresent = model.rxns(idx);
rxnAbsent = setdiff(rxnList, rxnPresent);

%% FIND MIN/MAX of RXNS
[min, max] = fluxVariability(model);
min_rxnPresent = min(idx);
max_rxnPresent = max(idx);

%% PRINT RESULTS TO TABLE
VarNames = {'rxnPresent', 'min_rxnPresent', 'max_rxnPresent'};
T = table(rxnPresent, min_rxnPresent, max_rxnPresent,'VariableNames',VarNames);

%%
toc;

end
