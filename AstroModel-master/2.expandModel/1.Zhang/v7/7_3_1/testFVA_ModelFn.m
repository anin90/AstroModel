%below code displays minFlux & maxFlux thru 'rxns' for FVA (no objective specified):

tic;
model = MBA_model;
model = changeObjective(model, 'biomass_maintenance', 0);
%[minFluxModel, maxFluxModel] = fluxVariability(model,0); % '0' for no objective specifying
rxns = {'ASPTA','ALATA_L','GLUDxi','HMR_9802','ACONT','GLUDC','GTHS','r0354'};
%rxns = {'biomass_maintenance','sink_oaa[c]','sink_pyr[c]','DM_atp_c_','sink_nh4[c]','DM_atp[c]','sink_4abut[c]','sink_gthrd[c]','DM_atp[c]','sink_succoa[m]','sink_accoa[c]','sink_accoa[c]','sink_r5p[c]','sink_pmtcoa[c]','sink_amp[c]','sink_imp[c]','sink_gmp[c]','sink_txa2[r]','sink_nrpphr[c]','sink_adrnl[c]','sink_melatn[c]','sink_hista[c]','sink_nh4[c]','sink_acac[m]','sink_pyr[m]','sink_gln-L[c]','DM_atp[c]','sink_coa[m]','sink_accoa[m]','sink_dlnlcgcoa[c]','sink_ach[c]','sink_creat[c]','sink_pcreat[c]','sink_oaa[m]','ABTArm','GABA aminotransferase','METAT','sink_crtn[c]','sink_prostgh2[c]','sink_prostgd2[r]','sink_prostge2[r]','sink_prostgi2[r]','sink_leuktrE4[c]','sink_C06314[c]','sink_3mox4hoxm[c]'
%};
for ind = 1:length(rxns)
	rnxIndex = find(strcmp(model.rxns, rxns(ind)));
	%disp(minFluxModel(rnxIndex));
	%disp(maxFluxModel(rnxIndex));
    %find(strcmp(model.rxns,rxns(ind)))
    disp(printRxnFormula(model,rxns(ind)));
    printFluxBounds(model,rxns(ind));
end

clear model
clear minFluxModel
clear maxFluxModel
clear maxFlux
clear rnxIndex
clear ind
clear rxns
toc;