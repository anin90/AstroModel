%below code displays minFlux & maxFlux thru 'rxns' for FVA (no objective specified):
tic;
model = MBA_model;
model = changeObjective(model, 'biomass_maintenance', 0);
[minFluxModel, maxFluxModel] = fluxVariability(model,0); % '0' for no objective specifying
rxns = {'biomass_maintenance','EX_glc_D[e]','EX_o2[e]','O2t','DM_atp_c_','CAATPS', 'EX_gln_L[s]', 'EX_glu_L[s]'};
for ind = 1:length(rxns)
	rnxIndex = find(strcmp(model.rxns, rxns(ind)));
	disp(minFluxModel(rnxIndex));
	disp(maxFluxModel(rnxIndex));
end
ExcRxns = model.rxns(findExcRxns(model));
ExcMax = intersect(model.rxns(maxFluxModel~=0),ExcRxns);
ExcMin = intersect(model.rxns(minFluxModel~=0),ExcRxns);
ExcUnion = union(ExcMax,ExcMin);
ExcIntersect = intersect(ExcMax,ExcMin);

clear model
clear minFluxModel
clear maxFluxModel
clear maxFlux
clear rnxIndex
clear ind
clear rxns
%clear ExcRxns
%clear ExcMax
%clear ExcMin
%clear ExcUnion
%clear ExcIntersect
toc;