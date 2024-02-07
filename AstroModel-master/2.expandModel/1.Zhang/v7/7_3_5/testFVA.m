%below code displays minFlux & maxFlux thru 'rxns' for FVA (no objective specified):
tic;
%%
model = MBA_model_7_3_4_ASM;
model = changeObjective(model, 'biomass_maintenance', 1);
[minFluxModel, maxFluxModel] = fluxVariability(model,0); % '0' for no objective specifying
rxns = {'biomass_maintenance','EX_ca2[e]','EX_cl[e]','EX_fol[e]','EX_ile_L[e]','EX_mg2[e]','EX_na1[e]','EX_no3[e]','EX_pnto_R[e]','EX_prgstrn[e]','EX_selni[e]','EX_so4[e]','EX_zn2[e]'};
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

%%
clearvars -except MBA_model_7_3 MBA_model_7_3_4_ASM fluxInconsistentRxns
%%
toc;