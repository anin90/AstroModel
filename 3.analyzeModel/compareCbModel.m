function[FluxDiff] = compareCbModel(model_a, model_b)
tic;

changeCobraSolver('gurobi','all'); clear ans;

%% FVA model_a & model_b
model_a = changeObjective(model_a, 'biomass_maintenance', 0);
[minFlux_a, maxFlux_a] = fluxVariability(model_a);

model_b = changeObjective(model_b, 'biomass_maintenance', 0);
[minFlux_b, maxFlux_b] = fluxVariability(model_b);

%% EXTRACT RXNS COMMON TO model_a & model_b
VarNames_a = {'Rxns_model_a', 'minFlux_a', 'maxFlux_a'};
VarNames_b = {'Rxns_model_b', 'minFlux_b', 'maxFlux_b'};
flux_a = table(model_a.rxns, minFlux_a, maxFlux_a, 'VariableNames',VarNames_a);
flux_b = table(model_b.rxns, minFlux_b, maxFlux_b, 'VariableNames',VarNames_b);
rxnIntersect = intersect(flux_a.Rxns_model_a, flux_b.Rxns_model_b);
fluxIntersect_a = flux_a(ismember(flux_a.Rxns_model_a, rxnIntersect)==1,:);
fluxIntersect_b = flux_b(ismember(flux_b.Rxns_model_b, rxnIntersect)==1,:);
fluxIntersect_a_sorted = sortrows(fluxIntersect_a); fluxIntersect_b_sorted = sortrows(fluxIntersect_b);
fluxIntersect_a_vs_b_logical = strcmp(fluxIntersect_a_sorted.Rxns_model_a, fluxIntersect_b_sorted.Rxns_model_b);
fluxIntersect_a_vs_b_logical_all = all(fluxIntersect_a_vs_b_logical);

%% CALCULATE FLUX SPAN RATIO (FSr) for RXNS COMMON to model_a & model_b
FS_a = fluxIntersect_a_sorted.maxFlux_a-fluxIntersect_a_sorted.minFlux_a;
FS_b = bsxfun(@minus, fluxIntersect_b_sorted.maxFlux_b, fluxIntersect_b_sorted.minFlux_b);
FSr = bsxfun(@rdivide, abs(FS_a), abs(FS_b));

%% REACTIONS WITH DIFFERENTIAL FSr
FSr_significant_L = rxnIntersect(FSr>1.5);
FSr_significant_H = rxnIntersect(FSr<0.8);

%% SUBSYSTEMS DISRUPTED
% FSr_significant_L
[tf,loc] = ismember(model_a.rxns, FSr_significant_L); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
subSystem_L = model_a.subSystems(idx);

% FSr_significant_H
[tf,loc] = ismember(model_a.rxns, FSr_significant_H); [~,p] = sort(loc(tf)); idx = find(tf); idx = idx(p);
subSystem_H = model_a.subSystems(idx);

%% PRINT RESULTS TO TABLE
VarNames_FS = {'Rxn', 'minFlux_a', 'maxFlux_a', 'FluxSpan_a', 'minFlux_b', 'maxFlux_b', 'FluxSpan_b', 'Flux_Span_Ratio'};
VarNames_L = {'FSr_significant_L', 'subSystem_L'};
VarNames_H = {'FSr_significant_H', 'subSystem_H'};

fluxSpanTable = table(rxnIntersect, fluxIntersect_a_sorted.minFlux_a, fluxIntersect_a_sorted.maxFlux_a, FS_a, fluxIntersect_b_sorted.minFlux_b, fluxIntersect_b_sorted.maxFlux_b, FS_b, FSr, 'VariableNames',VarNames_FS);
fluxGoneLow = table(FSr_significant_L, subSystem_L, 'VariableNames',VarNames_L);
fluxGoneHigh = table(FSr_significant_H, subSystem_H, 'VariableNames',VarNames_H);

FluxDiff = struct('fluxSpanTable', {fluxSpanTable}, 'fluxGoneLow', {fluxGoneLow}, 'fluxGoneHigh', {fluxGoneHigh}, 'fluxIntersect_a_vs_b_logical_all', {fluxIntersect_a_vs_b_logical_all});

% Rand.fluxSpan(ismember(Rand.fluxSpan.rxnIntersect, {'MI1PP', 'INSTt4', 'PGS', 'PDE1', 'BPNT', 'PGMT', 'ADNCYC'})==1,:)
% Rand.fluxSpan(ismember(Rand.fluxSpan.rxnIntersect, intersect(findRxnsFromMets(iAstro_Primary_T,{'pyr[c]'}, 'consumersOnly', true), findRxnsFromMets(iAstro_Primary_TP,{'pyr[c]'}, 'consumersOnly', true)))==1,:)

%%
toc;
end