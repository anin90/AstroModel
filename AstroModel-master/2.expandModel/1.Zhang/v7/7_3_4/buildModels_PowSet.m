function[PowRxnList,model_KO] = multiModels2(model,filename);
tic;

%% Load WT model
changeCobraSolver('gurobi','all'); clear ans;
model = changeObjective(model, 'biomass_maintenance', 0);

%% perform FVA for WT model:
% [minFlux_WT, maxFlux_WT] = fluxVariability(model);

%% Load rxns for single_KO
delimiterIn = '\t';
rxnList = importdata(filename,delimiterIn);
PowRxnList = PowerSet(rxnList)';
%% build single_KO models and perform FVA
model_KO=cell(rxnList);
parfor ind = 1:length(rxnList');
    model_single_KO = changeRxnBounds(model,{rxnList(ind)},0,'b');
    model_KO{ind} = model_single_KO;
    % printFluxBounds(allModels{1},{'INSTt4'});
end

%%
toc;
end