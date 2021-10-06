function[min_idx, max_idx, idx] = Test4CuratedRxns(model, filename);
tic;

%% LOAD WT MODEL
changeCobraSolver('gurobi','all'); clear ans;
model = changeObjective(model, 'biomass_maintenance', 0);

%% IMPORT DATA
rxnList = importdata(filename);

%% 
[idx] = find(ismember(model.rxns,rxnList));
[minFlux, maxFlux] = fluxVariability(model);

min_idx = minFlux(idx);
max_idx = maxFlux(idx);

%%
toc;

end
