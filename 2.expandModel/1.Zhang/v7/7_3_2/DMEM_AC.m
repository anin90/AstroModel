tic;

model = MBA_model;
changeCobraSolver('gurobi','all');

%add DMEM missing exchanges into the model
model = addExchangeRxn(model, {'na1[e]'}, -1000, 1000);

%add DMEM missing transport Rxns
model = addReaction(model,'NAt','reactionFormula','na1[e] <=> na1[c] ', 'subSystem','Transport, extracellular', 'geneRule','6526.1 or 6523.1 or 6528.1');
model = addReaction(model,'CLCFTRte','reactionFormula','cl[c] <=> cl[e] ', 'subSystem','Transport, extracellular', 'geneRule','1080.1');
model = addReaction(model,'PYDXNtr','reactionFormula','pydxn[e] <=> pydxn[c] ', 'subSystem','Transport, extracellular', 'geneRule','');
model = addReaction(model,'AVITE1t','reactionFormula','avite1[e] -> avite1[c] ', 'subSystem','Transport, extracellular', 'geneRule','949.1 and 29881.1');

%fluxConsistency = verifyModel(model,'fluxConsistency',true); 
%fluxInconsistentRxns = model.rxns(fluxConsistency.fluxConsistency.consistentReactionBool~=1);

%define closed model
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

%%DMEM composition (Prerna_Bhalla & Swagatika_Sahoo, 3-Oct-2019)
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
model=changeRxnBounds(model,'EX_thm[e]',-0.000628089,'l');
model=changeRxnBounds(model,'EX_asn_L[e]',-0.000599212,'l');
model=changeRxnBounds(model,'EX_asp_L[e]',-0.000594779,'l');
model=changeRxnBounds(model,'EX_lys_L[e]',-0.000583155,'l');
model=changeRxnBounds(model,'EX_cortsn[e]',-0.000577981,'l');
model=changeRxnBounds(model,'EX_tststerone[e]',-0.000577861,'l');
model=changeRxnBounds(model,'EX_met_L[e]',-0.000538946,'l');
model=changeRxnBounds(model,'EX_glu_L[e]',-0.000538075,'l');
model=changeRxnBounds(model,'EX_his_L[e]',-0.000521523,'l');
model=changeRxnBounds(model,'EX_phe_L[e]',-0.000495896,'l');
model=changeRxnBounds(model,'EX_arg_L[e]',-0.000474548,'l');
model=changeRxnBounds(model,'EX_tyr_L[e]',-0.000460846,'l');
model=changeRxnBounds(model,'EX_ribflv[e]',-0.000442834,'l');
model=changeRxnBounds(model,'EX_pi[e]',-0.00043873,'l');
model=changeRxnBounds(model,'EX_trp_L[e]',-0.000390908,'l');
model=changeRxnBounds(model,'EX_fol[e]',-0.000377586,'l');
model=changeRxnBounds(model,'EX_pnto_R[e]',-0.000349726,'l');
model=changeRxnBounds(model,'EX_so4[e]',-0.00033809,'l');
model=changeRxnBounds(model,'EX_chsterol[e]',-0.00033395,'l');
model=changeRxnBounds(model,'EX_k[e]',-0.000299951,'l');
model=changeRxnBounds(model,'EX_bilirub[e]',-0.000285065,'l');
model=changeRxnBounds(model,'EX_inost[e]',-0.000166522,'l');
model=changeRxnBounds(model,'EX_ca2[e]',-0.000141384,'l');
model=changeRxnBounds(model,'EX_Nacasp[e]',-0.000136472,'l');
model=changeRxnBounds(model,'EX_retinol_cis_11[e]',-0.000130912,'l');
model=changeRxnBounds(model,'EX_chol[e]',-0.000119369,'l');
model=changeRxnBounds(model,'EX_urea[e]',-0.000111,'l');
model=changeRxnBounds(model,'EX_prgstrn[e]',-0.000106002,'l');

DMEM_model = model;

%optimizeCbModel(MBA_model,'max');
optimizeCbModel(DMEM_model,'max')

clearvars -except AstroModel_Lewis_2010 MBA_model MBA_model_7_3_1 TestSolutionName_Brain TestSolution_Brain Imbalanced_NonExcRxns DMEM_model

toc;
