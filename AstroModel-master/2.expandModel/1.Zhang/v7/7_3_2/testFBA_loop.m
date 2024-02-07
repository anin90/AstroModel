%below code displays minFlux & maxFlux thru 'rxns' for FVA optimized for biomass_maintenance:
tic;
model = MBA_model;
%changeCobraSolver('gurobi','all');
rxns = {'biomass_maintenance','EX_glc_D[e]','EX_o2[e]','O2t','DM_atp_c_','CAATPS', 'EX_gln_L[s]', 'EX_glu_L[s]'};
for ind = 1:length(rxns)
	modelNew = changeObjective(model, rxns(ind));
	modelNewMin = optimizeCbModel(modelNew,'min');
	modelNewMax = optimizeCbModel(modelNew,'max');
	disp(modelNewMin.obj)
	disp(modelNewMax.obj)
end

clear model
clear rxns
clear ind
clear modelNew
clear modelNewMin
clear modelNewMax
toc;
