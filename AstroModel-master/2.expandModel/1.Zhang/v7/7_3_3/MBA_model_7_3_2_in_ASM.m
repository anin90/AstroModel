tic;
model = MBA_model_7_3_2;
%% ASM_AC_START
%% define closed model (as per Sanity check!)
modelClosed = model;
modelexchanges1 = strmatch('Ex_',modelClosed.rxns);
modelexchanges4 = strmatch('EX_',modelClosed.rxns);
% modelexchanges2 = strmatch('DM_',modelClosed.rxns);
% modelexchanges3 = strmatch('sink_',modelClosed.rxns);
% selExc = (find( full((sum(abs(modelClosed.S)==1,1) ==1) & (sum(modelClosed.S~=0) == 1))))';

% set bounds of all exchanges to: lb=0 & ub=1000
modelexchanges = unique([modelexchanges1; modelexchanges4]);
modelClosed.lb(find(ismember(modelClosed.rxns,modelClosed.rxns(modelexchanges))))=0;
% modelClosed.ub(find(ismember(modelClosed.rxns,modelClosed.rxns(modelexchanges))))=1000;
model = modelClosed;

%% close EX_rxns_AC:
% modelClosed = model;
% modelEX_Rxns = strmatch('EX_',modelClosed.rxns);
% modelClosed.lb(find(ismember(modelClosed.rxns,modelClosed.rxns(modelEX_Rxns))))=0;
% modelClosed.ub(find(ismember(modelClosed.rxns,modelClosed.rxns(modelEX_Rxns))))=1000;
% model = modelClosed; 

%% addExchangeRxns as per media composition..
% add ASM_media_EX_unconstrained:
model=addExchangeRxn(model,{'ala_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'gln_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'arg_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'asn_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'asp_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'ca2[e]'},-1000,1000);
model=addExchangeRxn(model,{'chol[e]'},-1000,1000);
model=addExchangeRxn(model,{'cl[e]'},-1000,1000);
model=addExchangeRxn(model,{'cys_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'fe3[e]'},-1000,1000);
model=addExchangeRxn(model,{'fol[e]'},-1000,1000);
model=addExchangeRxn(model,{'glc_D[e]'},-1000,1000);
model=addExchangeRxn(model,{'glu_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'gly[e]'},-1000,1000);
model=addExchangeRxn(model,{'hco3[e]'},-1000,1000);
model=addExchangeRxn(model,{'his_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'ile_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'inost[e]'},-1000,1000);
model=addExchangeRxn(model,{'k[e]'},-1000,1000);
model=addExchangeRxn(model,{'leu_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'lys_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'met_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'mg2[e]'},-1000,1000);
model=addExchangeRxn(model,{'na1[e]'},-1000,1000);
model=addExchangeRxn(model,{'ncam[e]'},-1000,1000);
model=addExchangeRxn(model,{'no3[e]'},-1000,1000);
model=addExchangeRxn(model,{'phe_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'pi[e]'},-1000,1000);
model=addExchangeRxn(model,{'pnto_R[e]'},-1000,1000);
model=addExchangeRxn(model,{'prgstrn[e]'},-1000,1000);
model=addExchangeRxn(model,{'pro_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'ptrc[e]'},-1000,1000);
model=addExchangeRxn(model,{'pydx[e]'},-1000,1000);
model=addExchangeRxn(model,{'pyr[e]'},-1000,1000);
model=addExchangeRxn(model,{'ribflv[e]'},-1000,1000);
model=addExchangeRxn(model,{'selni[e]'},-1000,1000);
model=addExchangeRxn(model,{'ser_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'so4[e]'},-1000,1000);
model=addExchangeRxn(model,{'thm[e]'},-1000,1000);
model=addExchangeRxn(model,{'thr_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'trp_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'tyr_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'val_L[e]'},-1000,1000);
model=addExchangeRxn(model,{'zn2[e]'},-1000,1000);

% now update ASM_bounds, as per TFS:
model=changeRxnBounds(model,'EX_ala_L[e]',-200,'l');
model=changeRxnBounds(model,'EX_gln_L[e]',-200,'l');
model=changeRxnBounds(model,'EX_arg_L[e]',-0.39810428,'l');
model=changeRxnBounds(model,'EX_asn_L[e]',-0.0055333334,'l');
model=changeRxnBounds(model,'EX_asp_L[e]',-10,'l');
model=changeRxnBounds(model,'EX_ca2[e]',-1.8018018,'l');
model=changeRxnBounds(model,'EX_chol[e]',-0.028571429,'l');
model=changeRxnBounds(model,'EX_cl[e]',-51.724136,'l');
model=changeRxnBounds(model,'EX_cys_L[e]',-0.2603306,'l');
model=changeRxnBounds(model,'EX_fe3[e]',-0.00024752476,'l');
model=changeRxnBounds(model,'EX_fol[e]',-0.009070295,'l');
model=changeRxnBounds(model,'EX_glc_D[e]',-25,'l');
model=changeRxnBounds(model,'EX_glu_L[e]',-10,'l');
model=changeRxnBounds(model,'EX_gly[e]',-10,'l');
model=changeRxnBounds(model,'EX_hco3[e]',-26.190475,'l');
model=changeRxnBounds(model,'EX_his_L[e]',-0.2,'l');
model=changeRxnBounds(model,'EX_ile_L[e]',-0.8015267,'l');
model=changeRxnBounds(model,'EX_inost[e]',-0.04,'l');
model=changeRxnBounds(model,'EX_k[e]',-5.3333335,'l');
model=changeRxnBounds(model,'EX_leu_L[e]',-0.8015267,'l');
model=changeRxnBounds(model,'EX_lys_L[e]',-0.7978142,'l');
model=changeRxnBounds(model,'EX_met_L[e]',-0.20134228,'l');
model=changeRxnBounds(model,'EX_mg2[e]',-0.8136842,'l');
model=changeRxnBounds(model,'EX_na1[e]',-51.724136,'l');
model=changeRxnBounds(model,'EX_ncam[e]',-0.032786883,'l');
model=changeRxnBounds(model,'EX_no3[e]',-0.00024752476,'l');
model=changeRxnBounds(model,'EX_phe_L[e]',-0.4,'l');
model=changeRxnBounds(model,'EX_pi[e]',-0.9057971,'l');
model=changeRxnBounds(model,'EX_pnto_R[e]',-0.008385744,'l');
model=changeRxnBounds(model,'EX_prgstrn[e]',-0.0020033708,'l');
model=changeRxnBounds(model,'EX_pro_L[e]',-10,'l');
model=changeRxnBounds(model,'EX_ptrc[e]',-10.006211,'l');
model=changeRxnBounds(model,'EX_pydx[e]',-0.019607844,'l');
model=changeRxnBounds(model,'EX_pyr[e]',-0.22727273,'l');
model=changeRxnBounds(model,'EX_ribflv[e]',-0.0010638298,'l');
model=changeRxnBounds(model,'EX_selni[e]',-0.0030057803,'l');
model=changeRxnBounds(model,'EX_ser_L[e]',-10,'l');
model=changeRxnBounds(model,'EX_so4[e]',-0.0006736111,'l');
model=changeRxnBounds(model,'EX_thm[e]',-0.011869436,'l');
model=changeRxnBounds(model,'EX_thr_L[e]',-0.79831934,'l');
model=changeRxnBounds(model,'EX_trp_L[e]',-0.078431375,'l');
model=changeRxnBounds(model,'EX_tyr_L[e]',-0.39779004,'l');
model=changeRxnBounds(model,'EX_val_L[e]',-0.8034188,'l');
model=changeRxnBounds(model,'EX_zn2[e]',-0.0006736111,'l');

% now update ASM_bounds on overlapping synaptic exchanges, as per TFS:
% model=changeRxnBounds(model,'EX_glu_L[s]',-10,'l');
% model=changeRxnBounds(model,'EX_gln_L[s]',-200,'l');
% model=changeRxnBounds(model,'EX_k[s]',-5.3333335,'l');

%% check transport reactions for ASM exhcanges..
Cyto_trspRxns.gln_L=findTrspRxnFromMet(model, {'gln_L[e]'}, 'c');
Cyto_trspRxns.arg_L=findTrspRxnFromMet(model, {'arg_L[e]'}, 'c');
Cyto_trspRxns.asn_L=findTrspRxnFromMet(model, {'asn_L[e]'}, 'c');
Cyto_trspRxns.asp_L=findTrspRxnFromMet(model, {'asp_L[e]'}, 'c');
Cyto_trspRxns.ca2=findTrspRxnFromMet(model, {'ca2[e]'}, 'c');
Cyto_trspRxns.chol=findTrspRxnFromMet(model, {'chol[e]'}, 'c');
Cyto_trspRxns.cl=findTrspRxnFromMet(model, {'cl[e]'}, 'c');
Cyto_trspRxns.cys_L=findTrspRxnFromMet(model, {'cys_L[e]'}, 'c');
Cyto_trspRxns.fe3=findTrspRxnFromMet(model, {'fe3[e]'}, 'c');
Cyto_trspRxns.fol=findTrspRxnFromMet(model, {'fol[e]'}, 'c');
Cyto_trspRxns.glc_D=findTrspRxnFromMet(model, {'glc_D[e]'}, 'c');
Cyto_trspRxns.glu_L=findTrspRxnFromMet(model, {'glu_L[e]'}, 'c');
Cyto_trspRxns.gly=findTrspRxnFromMet(model, {'gly[e]'}, 'c');
Cyto_trspRxns.hco3=findTrspRxnFromMet(model, {'hco3[e]'}, 'c');
Cyto_trspRxns.his_L=findTrspRxnFromMet(model, {'his_L[e]'}, 'c');
Cyto_trspRxns.ile_L=findTrspRxnFromMet(model, {'ile_L[e]'}, 'c');
Cyto_trspRxns.inost=findTrspRxnFromMet(model, {'inost[e]'}, 'c');
Cyto_trspRxns.k=findTrspRxnFromMet(model, {'k[e]'}, 'c');
Cyto_trspRxns.leu_L=findTrspRxnFromMet(model, {'leu_L[e]'}, 'c');
Cyto_trspRxns.lys_L=findTrspRxnFromMet(model, {'lys_L[e]'}, 'c');
Cyto_trspRxns.met_L=findTrspRxnFromMet(model, {'met_L[e]'}, 'c');
Cyto_trspRxns.mg2=findTrspRxnFromMet(model, {'mg2[e]'}, 'c');
Cyto_trspRxns.na1=findTrspRxnFromMet(model, {'na1[e]'}, 'c');
Cyto_trspRxns.ncam=findTrspRxnFromMet(model, {'ncam[e]'}, 'c');
Cyto_trspRxns.no3=findTrspRxnFromMet(model, {'no3[e]'}, 'c');
Cyto_trspRxns.phe_L=findTrspRxnFromMet(model, {'phe_L[e]'}, 'c');
Cyto_trspRxns.pi=findTrspRxnFromMet(model, {'pi[e]'}, 'c');
Cyto_trspRxns.pnto_R=findTrspRxnFromMet(model, {'pnto_R[e]'}, 'c');
Cyto_trspRxns.prgstrn=findTrspRxnFromMet(model, {'prgstrn[e]'}, 'c');
Cyto_trspRxns.pro_L=findTrspRxnFromMet(model, {'pro_L[e]'}, 'c');
Cyto_trspRxns.ptrc=findTrspRxnFromMet(model, {'ptrc[e]'}, 'c');
Cyto_trspRxns.pydx=findTrspRxnFromMet(model, {'pydx[e]'}, 'c');
Cyto_trspRxns.pyr=findTrspRxnFromMet(model, {'pyr[e]'}, 'c');
Cyto_trspRxns.ribflv=findTrspRxnFromMet(model, {'ribflv[e]'}, 'c');
Cyto_trspRxns.selni=findTrspRxnFromMet(model, {'selni[e]'}, 'c');
Cyto_trspRxns.ser_L=findTrspRxnFromMet(model, {'ser_L[e]'}, 'c');
Cyto_trspRxns.so4=findTrspRxnFromMet(model, {'so4[e]'}, 'c');
Cyto_trspRxns.thm=findTrspRxnFromMet(model, {'thm[e]'}, 'c');
Cyto_trspRxns.thr_L=findTrspRxnFromMet(model, {'thr_L[e]'}, 'c');
Cyto_trspRxns.trp_L=findTrspRxnFromMet(model, {'trp_L[e]'}, 'c');
Cyto_trspRxns.tyr_L=findTrspRxnFromMet(model, {'tyr_L[e]'}, 'c');
Cyto_trspRxns.val_L=findTrspRxnFromMet(model, {'val_L[e]'}, 'c');
Cyto_trspRxns.zn2=findTrspRxnFromMet(model, {'zn2[e]'}, 'c');

trspRxns.metBool = structfun(@isempty, Cyto_trspRxns); % will return '1' if empty
trspRxns.met = fieldnames(Cyto_trspRxns);

%% add missing TransportRxns to bring media metabolites into cytosol..
model = addReaction(model,'ARGt4','reactionFormula','arg_L[e] + na1[e] -> arg_L[c] + na1[c] ', 'subSystem','Transport, extracellular', 'geneRule','11254.1');
model = addReaction(model,'CHOLt4','reactionFormula','chol[e] + na1[e] <=> chol[c] + na1[c] ', 'subSystem','Transport, extracellular', 'geneRule','60482.1');
model = addReaction(model,'r0963','reactionFormula','fol[e] -> fol[c] ', 'subSystem','Transport, extracellular', 'geneRule','2348.1 or 2350.1 or 2352.1');
model = addReaction(model,'INSTt4','reactionFormula','inost[e] + na1[e] <=> inost[c] + na1[c] ', 'subSystem','Transport, extracellular', 'geneRule','6526.1');
model = addReaction(model,'LYSt4','reactionFormula','lys_L[e] + na1[e] -> lys_L[c] + na1[c] ', 'subSystem','Transport, extracellular', 'geneRule','11254.1 or 6584.1');
model = addReaction(model,'HMR_9586','reactionFormula','mg2[e] -> mg2[c] ', 'subSystem','Transport, extracellular', 'geneRule','254428.1 or 84102.1');
model = addReaction(model,'NCAMUP','reactionFormula','ncam[e] -> ncam[c] ', 'subSystem','Transport, extracellular', 'geneRule','');
model = addReaction(model,'PRGSTRNt','reactionFormula','prgstrn[e] <=> prgstrn[c] ', 'subSystem','Transport, extracellular', 'geneRule','');
model = addReaction(model,'PTCRTD','reactionFormula','ptrc[e] <=> ptrc[c] ', 'subSystem','Transport, extracellular', 'geneRule','');
model = addReaction(model,'PYDXtr','reactionFormula','pydx[e] <=> pydx[c] ', 'subSystem','Transport, extracellular', 'geneRule','');
model = addReaction(model,'THMtrbc','reactionFormula','thm[e] <=> thm[c] ', 'subSystem','Transport, extracellular', 'geneRule','');
model = addReaction(model,'VALt4','reactionFormula','na1[e] + val_L[e] -> na1[c] + val_L[c] ', 'subSystem','Transport, extracellular', 'geneRule','11254.1 or (340024.1 and 57393.1) or (340024.1 and 59272.1)');
model = addReaction(model,'r2073','reactionFormula','h[e] + zn2[e] -> h[c] + zn2[c] ', 'subSystem','Transport, extracellular', 'geneRule','6556.1');

%% link transported metabolites to internal reactions, if not already linked..
model = addReaction(model,'FOLR2','reactionFormula','fol[c] + nadph[c] -> dhf[c] + nadp[c] ', 'subSystem','Transport, extracellular', 'geneRule','1719.1');
model = addReaction(model,'NMNS','reactionFormula','h[c] + ncam[c] + prpp[c] -> nmn[c] + ppi[c] ', 'subSystem','Transport, extracellular', 'geneRule','10135.1');
model = addReaction(model,'TMDPK','reactionFormula','atp[c] + thm[c] -> amp[c] + h[c] + thmpp[c] ', 'subSystem','Transport, extracellular', 'geneRule','27010.1');

%% addition of misc. reactions for biomass_maintenance..

% model = addSinkReactions(model,{'h2o[c]'});
% model = addSinkReactions(model,{'h[c]'});
% model = addSinkReactions(model,{'atp[c]'});
% model = addSinkReactions(model,{'adp[c]'});
% model = addSinkReactions(model,{'pi[c]'});
% model = addSinkReactions(model,{'glu_L[c]'});
% model = addSinkReactions(model,{'asp_L[c]'});
% model = addSinkReactions(model,{'gtp[c]'});
% model = addSinkReactions(model,{'ala_L[c]'});
% model = addSinkReactions(model,{'asn_L[c]'});
% model = addSinkReactions(model,{'cys_L[c]'});
% model = addSinkReactions(model,{'gln_L[c]'});
% model = addSinkReactions(model,{'gly[c]'});
% model = addSinkReactions(model,{'ser_L[c]'});
% model = addSinkReactions(model,{'thr_L[c]'});
% model = addSinkReactions(model,{'lys_L[c]'});
% model = addSinkReactions(model,{'arg_L[c]'});
% model = addSinkReactions(model,{'met_L[c]'});
% model = addSinkReactions(model,{'pail_hs[c]'});
% model = addSinkReactions(model,{'pchol_hs[c]'});
% model = addSinkReactions(model,{'chsterol[c]'});
% model = addSinkReactions(model,{'utp[c]'});
% model = addSinkReactions(model,{'g6p[c]'});
% model = addSinkReactions(model,{'his_L[c]'});
% model = addSinkReactions(model,{'tyr_L[c]'});
% model = addSinkReactions(model,{'ile_L[c]'});
% model = addSinkReactions(model,{'leu_L[c]'});
% model = addSinkReactions(model,{'trp_L[c]'});
% model = addSinkReactions(model,{'phe_L[c]'});
% model = addSinkReactions(model,{'pro_L[c]'});
% model = addSinkReactions(model,{'ps_hs[c]'});
% model = addSinkReactions(model,{'val_L[c]'});

% model = addSinkReactions(model,{'ctp[c]'}); % essential
% model = addSinkReactions(model,{'pe_hs[c]'}); % essential
% model = addSinkReactions(model,{'pglyc_hs[c]'}); % essential
% model = addSinkReactions(model,{'clpn_hs[c]'}); % essential
% model = addSinkReactions(model,{'sphmyln_hs[c]'}); % essential

%% misc
model=changeRxnBounds(model,'EX_o2[e]',-1000,'l');
model=addExchangeRxn(model,{'cytd[e]'},-1000,1000);
model=addExchangeRxn(model,{'pe_hs[e]'},-1000,1000);
model = addReaction(model,'PEt','reactionFormula','pe_hs[e] <=> pe_hs[c] ', 'subSystem','Transport, extracellular', 'geneRule','');
model=changeRxnBounds(model,'EX_mag_hs[e]',-1000,'l');
model = addReaction(model,'CLS_hs','reactionFormula','cdpdag_hs[c] + pglyc_hs[c] -> clpn_hs[c] + cmp[c] + h[c] ', 'subSystem','Glycerophospholipid metabolism', 'geneRule','54675.1');
model = addReaction(model,'DSAT','reactionFormula','Rtotalcoa[c] + sphgn[c] -> coa[c] + dhcrm_hs[c] + h[c] ', 'subSystem','Sphingolipid metabolism', 'geneRule','');

%% ASM_AC_END
%%
MBA_model_7_3_2_ASM = model;

%fluxConsistency = verifyModel(MBA_model_7_3_3,'fluxConsistency',true); MBA_model_7_3_3.rxns(fluxConsistency.fluxConsistency.consistentReactionBool~=1);

%%
clearvars -except MBA_model_7_3_1 MBA_model_7_3_2 MBA_model_7_3_2_ASM
%%
toc;
















