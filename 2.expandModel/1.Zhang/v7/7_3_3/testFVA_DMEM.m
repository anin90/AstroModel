%below code displays minFlux & maxFlux thru 'rxns' for FVA (no objective specified):
tic;
model = MBA_model_7_3_2;
%% DMEM_AC_START
%% define closed model
modelClosed = model;
modelexchanges1 = strmatch('Ex_',modelClosed.rxns);
modelexchanges4 = strmatch('EX_',modelClosed.rxns);
modelexchanges2 = strmatch('DM_',modelClosed.rxns);
modelexchanges3 = strmatch('sink_',modelClosed.rxns);
selExc = (find( full((sum(abs(modelClosed.S)==1,1) ==1) & (sum(modelClosed.S~=0) == 1))))';

%set bounds of all exchanges to: lb=0 & ub=1000
modelexchanges = unique([modelexchanges1;modelexchanges2;modelexchanges3;modelexchanges4;selExc]);
modelClosed.lb(find(ismember(modelClosed.rxns,modelClosed.rxns(modelexchanges))))=0;
modelClosed.ub(find(ismember(modelClosed.rxns,modelClosed.rxns(modelexchanges))))=1000;
model = modelClosed;
%% add DMEM missing exchanges into the model
%% CAN'T BE ADDED! B'cuz mets absent in cytosol!
%model=addExchangeRxn(model,{'thm[e]'},-1000,1000); %NA_Ex;NA_t;NA_[c]
%model=addExchangeRxn(model,{'cortsn[e]'},-1000,1000); %NA_Ex;NA_t;NA_[c]
%model=addExchangeRxn(model,{'tststerone[e]'},-1000,1000); %NA_Ex;NA_t;NA_[c]
%model=addExchangeRxn(model,{'fol[e]'},-1000,1000); %NA_Ex;NA_t;NA_[c]
%model=addExchangeRxn(model,{'bilirub[e]'},-1000,1000); %NA_Ex;NA_t;NA_[c]
%model=addExchangeRxn(model,{'prgstrn[e]'},-1000,1000); %NA_Ex;NA_t;NA_[c]
%% add DMEM missing exchanges & transport into the model
model=addExchangeRxn(model,{'pydxn[e]'},-1000,1000); %NA_Ex;NA_t
model = addReaction(model,'PYDXNtr','reactionFormula','pydxn[e] <=> pydxn[c] ', 'subSystem','Transport, extracellular', 'geneRule','');

model=addExchangeRxn(model,{'prostge2[e]'},-1000,1000); %NA_Ex;NA_t
model = addReaction(model,'PROSTGE2t','reactionFormula','hco3[c] + prostge2[e] <=> hco3[e] + prostge2[c] ', 'subSystem','Transport, extracellular', 'geneRule','28231.1 or 6578.1 or 28232.1 or 6579.1');

model=addExchangeRxn(model,{'val_L[e]'},-1000,1000); %NA_Ex;NA_t
model = addReaction(model,'VALt4','reactionFormula','na1[e] + val_L[e] -> na1[c] + val_L[c] ', 'subSystem','Transport, extracellular', 'geneRule','11254.1 or (340024.1 and 57393.1) or (340024.1 and 59272.1)');

model=addExchangeRxn(model,{'lys_L[e]'},-1000,1000); %NA_Ex;NA_t
model = addReaction(model,'LYSt4','reactionFormula','lys_L[e] + na1[e] -> lys_L[c] + na1[c] ', 'subSystem','Transport, extracellular', 'geneRule','11254.1 or 6584.1');

model=addExchangeRxn(model,{'arg_L[e]'},-1000,1000); %NA_Ex;NA_t
model = addReaction(model,'ARGt4','reactionFormula','arg_L[e] + na1[e] -> arg_L[c] + na1[c] ', 'subSystem','Transport, extracellular', 'geneRule','11254.1');

model=addExchangeRxn(model,{'chsterol[e]'},-1000,1000); %NA_Ex;NA_t
model = addReaction(model,'r1050','reactionFormula','chsterol[e] <=> chsterol[c] ', 'subSystem','Transport, extracellular', 'geneRule','');

model=addExchangeRxn(model,{'inost[e]'},-1000,1000); %NA_Ex;NA_t
model = addReaction(model,'INSTt4','reactionFormula','inost[e] + na1[e] <=> inost[c] + na1[c] ', 'subSystem','Transport, extracellular', 'geneRule','6526.1');

model=addExchangeRxn(model,{'Nacasp[e]'},-1000,1000); %NA_Ex;NA_t
model = addReaction(model,'NACASPt','reactionFormula','3.0 na1[c] + Nacasp[c] <=> 3.0 na1[e] + Nacasp[e] ', 'subSystem','Transport, extracellular', 'geneRule','');

%model=addExchangeRxn(model,{'retinol_cis_11[e]'},-1000,1000); %NA_Ex;NA_t
%model = addReaction(model,'','reactionFormula','', 'subSystem','Transport, extracellular', 'geneRule','');

model=addExchangeRxn(model,{'chol[e]'},-1000,1000); %NA_Ex;NA_t
model = addReaction(model,'CHOLt4','reactionFormula','chol[e] + na1[e] <=> chol[c] + na1[c] ', 'subSystem','Transport, extracellular', 'geneRule','60482.1');
%% add DMEM missing exchanges into the model
model=addExchangeRxn(model,{'na1[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'cl[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'hco3[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'avite1[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'ser_L[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'prostge1[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'gln_L[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'thr_L[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'leu_L[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'asn_L[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'glu_L[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'his_L[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'tyr_L[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'ribflv[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'pnto_R[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'so4[e]'},-1000,1000); %NA_Ex
model=addExchangeRxn(model,{'ca2[e]'},-1000,1000); %NA_Ex
%% add other DMEM exchanges into the model (these might be already present in model_v7_3_2)
model=addExchangeRxn(model,{'glc_D[e]'},-1000,1000);
model=addExchangeRxn(model,{'gly[e]'},-1000,1000);
model=addExchangeRxn(model,{'creat[e]'},-1000,1000);
model=addExchangeRxn(model,{'ala_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'pro_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'cys_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'ile_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'asp_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'met_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'phe_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'pi[e]'},-1000,1000);
model=addExchangeRxn(model,{'trp_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'k[e]'},-1000,1000);
model=addExchangeRxn(model,{'urea[e]'},-1000,1000);

%% set constant bounds for synaptic exchanges
model=changeRxnBounds(model, 'EX_h[s]',-1000,1000);
model=changeRxnBounds(model, 'EX_glu_L[s]',-1000,0);    %only uptake
model=changeRxnBounds(model, 'EX_k[s]',-1000,1000);
model=changeRxnBounds(model, 'EX_gln_L[s]',0,1000);     %only release
model=changeRxnBounds(model, 'EX_lac_L[s]',0,1000);     %only release
model=changeRxnBounds(model, 'EX_ala_L[s]',0,1000);     %only release
model=changeRxnBounds(model, 'EX_cit[s]',0,1000);       %only release
model=changeRxnBounds(model, 'EX_h2o[s]',-1000,0);      %only uptake
model=changeRxnBounds(model, 'EX_nh4[s]',0,1000);       %only release

%% set DMEM bounds (Prerna_Bhalla & Swagatika_Sahoo, 3-Oct-2019)
model=changeRxnBounds(model, 'EX_na1[e]',-0.0107908,'l');
model=changeRxnBounds(model, 'EX_cl[e]',-0.003405409,'l');
model=changeRxnBounds(model, 'EX_hco3[e]',-0.001835165,'l');
model=changeRxnBounds(model, 'EX_glc_D[e]',-0.001331019,'l');
model=changeRxnBounds(model, 'EX_gly[e]',-0.001071271,'l');
model=changeRxnBounds(model, 'EX_pydxn[e]',-0.000985144,'l');
model=changeRxnBounds(model, 'EX_creat[e]',-0.000985028,'l');
model=changeRxnBounds(model, 'EX_avite1[e]',-0.000967395,'l');
model=changeRxnBounds(model, 'EX_ala_L[e]',-0.000888583,'l');
model=changeRxnBounds(model, 'EX_ser_L[e]',-0.000769956,'l');
model=changeRxnBounds(model, 'EX_prostge1[e]',-0.00070929,'l');
model=changeRxnBounds(model,'EX_prostge2[e]'-0.00070929,'l');
model=changeRxnBounds(model,'EX_val_L[e]',-0.000709227,'l');
model=changeRxnBounds(model,'EX_gln_L[e]',-0.000708203,'l');
model=changeRxnBounds(model,'EX_thr_L[e]',-0.00069818,'l');
model=changeRxnBounds(model,'EX_pro_L[e]',-0.000687625,'l');
model=changeRxnBounds(model,'EX_cys_L[e]',-0.000675083,'l');
model=changeRxnBounds(model,'EX_ile_L[e]',-0.000636882,'l');
model=changeRxnBounds(model,'EX_leu_L[e]',-0.000636882,'l');
model=changeRxnBounds(model,'EX_asn_L[e]',-0.000599212,'l');
model=changeRxnBounds(model,'EX_asp_L[e]',-0.000594779,'l');
model=changeRxnBounds(model,'EX_lys_L[e]',-0.000583155,'l');
model=changeRxnBounds(model,'EX_met_L[e]',-0.000538946,'l');
model=changeRxnBounds(model,'EX_glu_L[e]',-0.000538075,'l');
model=changeRxnBounds(model,'EX_his_L[e]',-0.000521523,'l');
model=changeRxnBounds(model,'EX_phe_L[e]',-0.000495896,'l');
model=changeRxnBounds(model,'EX_arg_L[e]',-0.000474548,'l');
model=changeRxnBounds(model,'EX_tyr_L[e]',-0.000460846,'l');
model=changeRxnBounds(model,'EX_ribflv[e]',-0.000442834,'l');
model=changeRxnBounds(model,'EX_pi[e]',-0.00043873,'l');
model=changeRxnBounds(model,'EX_trp_L[e]',-0.000390908,'l');
model=changeRxnBounds(model,'EX_pnto_R[e]',-0.000349726,'l');
model=changeRxnBounds(model,'EX_so4[e]',-0.00033809,'l');
model=changeRxnBounds(model,'EX_chsterol[e]',-0.00033395,'l');
model=changeRxnBounds(model,'EX_k[e]',-0.000299951,'l');
model=changeRxnBounds(model,'EX_inost[e]',-0.000166522,'l');
model=changeRxnBounds(model,'EX_ca2[e]',-0.000141384,'l');
model=changeRxnBounds(model,'EX_Nacasp[e]',-0.000136472,'l');
model=changeRxnBounds(model,'EX_chol[e]',-0.000119369,'l');
model=changeRxnBounds(model,'EX_urea[e]',-0.000111,'l');

%model=changeRxnBounds(model,'EX_thm[e]',-0.000628089,'l');
%model=changeRxnBounds(model,'EX_cortsn[e]',-0.000577981,'l');
%model=changeRxnBounds(model,'EX_tststerone[e]',-0.000577861,'l');
%model=changeRxnBounds(model,'EX_fol[e]',-0.000377586,'l');
%model=changeRxnBounds(model,'EX_bilirub[e]',-0.000285065,'l');
%model=changeRxnBounds(model,'EX_retinol_cis_11[e]',-0.000130912,'l');
%model=changeRxnBounds(model,'EX_prgstrn[e]',-0.000106002,'l');
%% Other requirements for model to run in DMEM:
% model=changeRxnBounds(model, 'EX_h2o[e]',-1000,1000);
% model=changeRxnBounds(model, 'EX_atp[e]',-1000,1000);
% model=changeRxnBounds(model, 'EX_pchol_hs[e]',-1000,1000);
% model=changeRxnBounds(model, 'EX_adp[e]',-1000,1000);

% model=changeRxnBounds(model, 'EX_met_L[e]',-1000,1000);
% model=changeRxnBounds(model, 'EX_cys_L[e]',-1000,1000);
% model=changeRxnBounds(model, 'EX_phe_L[e]',-1000,1000);
% model=changeRxnBounds(model, 'EX_gly[e]',-1000,1000);
% model=changeRxnBounds(model, 'EX_pi[e]',-1000,1000);
% model=changeRxnBounds(model, 'EX_ala_L[e]',-1000,1000);
% model=changeRxnBounds(model, 'EX_ile_L[e]',-1000,1000);
% model=changeRxnBounds(model, 'EX_pro_L[e]',-1000,1000);
% model=changeRxnBounds(model, 'EX_trp_L[e]',-1000,1000);
% model=changeRxnBounds(model, 'EX_chsterol[e]',-1000,1000);
% model=changeRxnBounds(model, 'EX_asp_L[e]',-1000,1000);
% model=changeRxnBounds(model, 'EX_adp[e]',-1000,1000);
% model=changeRxnBounds(model, 'EX_adp[e]',-1000,1000);

%% DMEM_AC_END
%% 
model = model;
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
%%
clearvars -except AstroModel_Lewis_2010 MBA_model_7_3_1 MBA_model_7_3_2 Imbalanced_NonExcRxns
%%
toc;