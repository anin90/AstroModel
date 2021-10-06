function [TestSolution,TestSolutionName,TestedRxns,PercTestedRxns] = test4HumanFctExt(model,test,optionSinks)
% test for the ~50 human brain functions
%
% USAGE:
%     [TestSolution,TestSolutionName,TestedRxns,PercTestedRxns] = test4HumanFctExt(model,test,optionSinks)
%
% INPUT:
%    model:             model structure (Recon1, with desired in silico condition)
%    test:              possible statements: Recon1, IECori, IEC, all (default)
%                       (choose IECori if you intend to test the IEC model OR a model that
%                       contains lumen ('u') as compartment otw choose IEC);
%                       all check for Recon1 and IEC
%    option:            if true = set sink reactions to 0 (default, leave unchanged).
%                       Note that all lb's of exchanges and demands will be set to 0
%
% OUTPUT:
%    TestSolution:      array containing the optimal value for the different tests
%    TestSolutionName:  array containing the names  for the different tests
%
% .. Authors:
%       - Ines Thiele, 09/05/09
%       - MKA, 03/04/12 some of the reaction names have changed in newer versions of Recon1
%         and Recon2. Comment setup til line 146 if using an old version.
%       - MKA, 24/05/12 finds correct EX_reactions and changes these to zero
%       - IT, 07/20/12 added tests for sIEC model
%       - AH, 07/12/17 minor changes to constraints that were resulting in infeasible models
%       - Anirudh S. Chellappa, 03/03/20
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
% model=changeRxnBounds(model,'EX_o2[e]',-1000,'l');
model=addExchangeRxn(model,{'cytd[e]'},-1000,1000);
model=addExchangeRxn(model,{'pe_hs[e]'},-1000,1000);
model = addReaction(model,'PEt','reactionFormula','pe_hs[e] <=> pe_hs[c] ', 'subSystem','Transport, extracellular', 'geneRule','');
model=changeRxnBounds(model,'EX_mag_hs[e]',-1000,'l');
model = addReaction(model,'CLS_hs','reactionFormula','cdpdag_hs[c] + pglyc_hs[c] -> clpn_hs[c] + cmp[c] + h[c] ', 'subSystem','Glycerophospholipid metabolism', 'geneRule','54675.1');
model = addReaction(model,'DSAT','reactionFormula','Rtotalcoa[c] + sphgn[c] -> coa[c] + dhcrm_hs[c] + h[c] ', 'subSystem','Sphingolipid metabolism', 'geneRule','');

%% ASM_AC_END
%%
if nargin<2
    test = 'all';
end
if nargin<3
    optionSinks = 0; % do not close
end

if optionSinks
    % close sink reactions
    % model.lb(strmatch('sink_',model.rxns))=0;
end

TestSolution = [];
%%
% S ='';
%  S = warning('QUERY','VERBOSE');
%% Setup

% for organ atlas derived from Harvey only
if strcmp(test,'Harvey')
    model.rxns = regexprep(model.rxns,'\[bc\]','\(e\)');
end

diary('Test4Functions_diary.txt');
TestedRxns =[];
tol = 1e-6;
% fixes the model met names in cases the compartments are not given with ()
model.mets = regexprep(model.mets,'[','(');
model.mets = regexprep(model.mets,']',')');
model.mets = regexprep(model.mets,'_','-');
model.mets = regexprep(model.mets,'-FSLASH-','/');

model.rxns = regexprep(model.rxns,'\[','\(');
model.rxns = regexprep(model.rxns,'\]','\)');

% replace reaction names
new = {'DM_atp_c_'
    'EX_gln_L(e)'
    'EX_glu_L(e)'
    'EX_lac_L(e)'
    'EX_pro_L(e)'
    'EX_cys_L(e)'
    'EX_lys_L(e)'
    'EX_arg_L(e)'
    'EX_his_L(e)'
    'EX_glc_D(e)'
    'CYOR_u10m'
    'NADH2_u10m'
    'EX_4hpro(e)'
    % due to innermitochondrial membrane representation in recon3
    'ASPGLUmi'
    'ATPS4mi'
    'CYOR_u10mi'
    'Htmi'
    'NADH2_u10mi'
    'CYOOm3i'
    'CYOOm2i'};

original = {'DM_atp(c)'
    'EX_gln-L(e)'
    'EX_glu-L(e)'
    'EX_lac-L(e)'
    'EX_pro-L(e)'
    'EX_cys-L(e)'
    'EX_lys-L(e)'
    'EX_arg-L(e)'
    'EX_his-L(e)'
    'EX_glc(e)'
    'CYOR-u10m'
    'NADH2-u10m'
    'EX_4HPRO'
    'ASPGLUm'
    'ATPS4m'
    'CYOR-u10m'
    'Htm'
    'NADH2-u10m'
    'CYOOm3'
    'CYOOm2'};

for i=1:length(new)
    A = find(ismember(model.rxns,new(i,1)));
    model.rxns(A,1)= original(i,1);
end

%replace metabolite names

new_mets = {'Ser-Gly-Ala-X-Gly(r)'
    'Ser-Thr(g)'
    'Ser-Thr(l)'
    'ksii-core2(g)'
    'ksii-core4(g)'
    'ksii-core2(l)'
    'ksii-core4(l)'
    'cspg-a(l)'
    'cspg-b(l)'
    'cspg-c(l)'
    'cspg-d(l)'
    'cspg-e(l)'
    'cspg-a(g)'
    'cspg-b(g)'
    'cspg-c(g)'
    'cspg-d(g)'
    'cspg-e(g)'
    'galgluside-hs(g)'
    'gluside-hs(g)'
    'galgalgalthcrm-hs(g)'
    'acgagbside-hs(g)'
    'acnacngalgbside-hs(g)'
    'gd1b2-hs(g)'
    'gd1c-hs(g)'
    'gq1balpha-hs(g)'
    'dag-hs(c)'
    'pe-hs(c)'
    'tag-hs(c)'
    'cs-pre(g)'
    'crmp-hs(c)'
    'sphmyln-hs(c)'
    'pail-hs(c)'
    'pail45p-hs(c)'
    'pail4p-hs(c)'
    'dolichol-L(c)'
    'dolmanp-L(r)'
    'dolichol-U(c)'
    'dolmanp-U(r)'
    'dolichol-L(r)'
    'dolichol-U(r)'
    'gpi-prot-hs(r)'
    'g3m8mpdol-L(r)'
    'g3m8mpdol-U(r)'
    'gp1c-hs(g)'
    'dsTn-antigen(g)'
    'sTn-antigen(g)'
    'Tn-antigen(g)'
    };
original_mets = {'Ser-Gly/Ala-X-Gly(r)'
    'Ser/Thr(g)'
    'Ser/Thr(l)'
    'ksii_core2(g)'
    'ksii_core4(g)'
    'ksii_core2(l)'
    'ksii_core4(l)'
    'cspg_a(l)'
    'cspg_b(l)'
    'cspg_c(l)'
    'cspg_d(l)'
    'cspg_e(l)'
    'cspg_a(g)'
    'cspg_b(g)'
    'cspg_c(g)'
    'cspg_d(g)'
    'cspg_e(g)'
    'galgluside_hs(g)'
    'gluside_hs(g)'
    'galgalgalthcrm_hs(g)'
    'acgagbside_hs(g)'
    'acnacngalgbside_hs(g)'
    'gd1b2_hs(g)'
    'gd1c_hs(g)'
    'gq1balpha_hs(g)'
    'dag_hs(c)'
    'pe_hs(c)'
    'tag_hs(c)'
    'cs_pre(g)'
    'crmp_hs(c)'
    'sphmyln_hs(c)'
    'pail_hs(c)'
    'pail45p_hs(c)'
    'pail4p_hs(c)'
    'dolichol_L(c)'
    'dolmanp_L(r)'
    'dolichol_U(c)'
    'dolmanp_U(r)'
    'dolichol_L(r)'
    'dolichol_U(r)'
    'gpi_prot_hs(r)'
    'g3m8mpdol_L(r)'
    'g3m8mpdol_U(r)'
    'gp1c_hs(g)'
    'dsTn_antigen(g)'
    'sTn_antigen(g)'
    'Tn_antigen(g)'
    };
for i=1:length(new_mets)
    met = new_mets(i,1);
    A = find(ismember(model.mets,met));
    model.mets(A,1)= original_mets(i,1);
end

for i=1:length(new_mets)
    M = regexprep(new_mets(i,1),'(','[');
    M = regexprep(M,')',']');
    met = new_mets(i,1);
    A = find(ismember(M,met));
    model.mets(A,1)= original_mets(i,1);
end

% close sink reactions
% model.lb(strmatch('DM_',model.rxns))=0;

%aerobic
model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=0;

model.c(find(model.c)) = 0;
modelOri = model;
k = 1;
RPMI_composition={'EX_ala_L(e)','EX_arg-L(e)','EX_asn_L(e)','EX_asp_L(e)','EX_cys-L(e)','EX_gln-L(e)','EX_glu-L(e)','EX_gly(e)','EX_his-L(e)','EX_ile_L(e)','EX_leu_L(e)','EX_lys-L(e)','EX_met_L(e)','EX_phe_L(e)','EX_4HPRO','EX_pro-L(e)','EX_ser_L(e)','EX_thr_L(e)','EX_trp_L(e)','EX_tyr_L(e)','EX_val_L(e)','EX_ascb_L(e)','EX_btn(e)','EX_chol(e)','EX_pnto_R(e)','EX_fol(e)','EX_ncam(e)','EX_pydxn(e)','EX_ribflv(e)','EX_thm(e)','EX_cbl1(e)','EX_inost(e)','EX_ca2(e)','EX_fe3(e)','EX_k(e)','EX_hco3(e)','EX_na1(e)','EX_pi(e)','EX_glc(e)','EX_hxan(e)','EX_lnlc(e)','EX_lipoate(e)','EX_ptrc(e)','EX_pyr(e)','EX_thymd(e)','EX_etha(e)','EX_gthrd(e)'};

if strcmp(test,'Recon1') || strcmp(test,'all') || strcmp(test,'Harvey')
    
	%% Generic
    %% "Human Recon 1 human  biomass"
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'biomass_maintenance'))=1;
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Human Recon 1 test human biomass (maintenance)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Mets for below rxns are present in media:
    %% BrainTasks_Old
    %% asp-L -> oaa / transamination pair
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'asp-L(c)','oaa(c)'},[-10 -10; 0 1000]);
    model.c(ismember(model.rxns,'sink_oaa(c)'))=1;
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'asp-L -> oaa / transamination pair';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% ala-L(c) -> pyr(c) / alanine transaminase
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'ala-L(c)','pyr(c)'},[-200 -200; 0 1000]);
    model.c(ismember(model.rxns,'sink_pyr(c)'))=1;
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ala-L(c) -> pyr(c) / alanine transaminase';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA    

	%% ATP max, aerobic, glu-L/ glutamate dehydrogenase -> energy by TCA
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glu-L(e)'))=-10;model.ub(ismember(model.rxns,'EX_glu-L(e)'))=-10;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    model.c(ismember(model.rxns,'DM_atp(c)'))=1;
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ATP max, aerobic, glu-L/ glutamate dehydrogenase -> energy by TCA';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% gln-L(c) -> nh4(c)/ glutaminase
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'gln-L(c)','nh4(c)'},[-200 -200; 0 1000]);
    model.lb(ismember(model.rxns,'EX_nh4(e)'))=0; model.ub(ismember(model.rxns,'EX_nh4(e)'))=1000;
    model.c(ismember(model.rxns,'sink_nh4(c)'))=1;
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'gln-L(c) -> nh4(c)/ glutaminase';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% glu-L(c) -> 4abut(c)/ GABA synthesis from glutamate
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'glu-L(c)','4abut(c)'},[-10 -10; 0 1000]);
    model.c(ismember(model.rxns,'sink_4abut(c)'))=1;
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero')
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glu-L(c) -> 4abut(c)/ GABA synthesis from glutamate';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA    

    %% cys_L + glu_L + gly -> ghtrd/ glutathione synthesis
	model = modelOri;
	model.c(find(model.c)) = 0;
	[model] = addSinkReactions(model,{'cys-L(c)','glu-L(c)','gly(c)','gthrd(c)'},[-0.2603306 -0.2603306; -10 -10; -10 -10; 0 1000]);
	model.c(ismember(model.rxns,'sink_gthrd(c)'))=1;
	if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'cys_L + glu_L + gly -> ghtrd/ glutathione synthesis';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA   
    
	%% ATP max, aerobic, glc/ glycolysis
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glc(e)'))=-25;model.ub(ismember(model.rxns,'EX_glc(e)'))=-25;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    model.c(ismember(model.rxns,'DM_atp(c)'))=1;
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ATP max, aerobic, glc/ glycolysis';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% val-L(c) -> succoa(m)/ valine degradation
	model = modelOri;
	model.c(find(model.c)) = 0;
	[model] = addSinkReactions(model,{'val-L(c)','succoa(m)'},[-0.8034188 -0.8034188; 0 1000]);
	model.c(ismember(model.rxns,'sink_succoa(m)'))=1;
	if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'val-L(c) -> succoa(m)/ valine degradation';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA  
    
    %% leu-L(c) -> accoa(c)/ leucine degradation
	model = modelOri;
	model.c(find(model.c)) = 0;
	[model] = addSinkReactions(model,{'leu-L(c)','accoa(c)'},[-0.8015267 -0.8015267; 0 1000]);
	model.c(ismember(model.rxns,'sink_accoa(c)'))=1;
	if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'leu-L(c) -> accoa(c)/ leucine degradation';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA  
    
    %% ile-L(c) -> accoa(c)/ isoleucine degradation
	model = modelOri;
	model.c(find(model.c)) = 0;
	[model] = addSinkReactions(model,{'ile-L(c)','accoa(c)'},[-0.8015267 -0.8015267; 0 1000]);
	model.c(ismember(model.rxns,'sink_accoa(c)'))=1;
	if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ile-L(c) -> accoa(c)/ isoleucine degradation';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA  
    
    %% tyr_L(c) -> nrpphr(c)/ norepinephrine synthesis
	model = modelOri;
	model.c(find(model.c)) = 0;
	[model] = addSinkReactions(model,{'tyr-L(c)','nrpphr(c)'},[-0.39779004 -0.39779004; 0 1000]);
	model.c(ismember(model.rxns,'sink_nrpphr(c)'))=1;
	if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'tyr_L(c) -> nrpphr(c)/ norepinephrine synthesis';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% tyr_L(c) -> adrnl(c)/ ephinephrine synthesis
	model = modelOri;
	model.c(find(model.c)) = 0;
	[model] = addSinkReactions(model,{'tyr-L(c)','adrnl(c)'},[-0.39779004 -0.39779004; 0 1000]);
	model.c(ismember(model.rxns,'sink_adrnl(c)'))=1;
	if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'tyr_L(c) -> adrnl(c)/ ephinephrine synthesis';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA  
    
    %% trp_L(c) -> melatn(c)
	model = modelOri;
	model.c(find(model.c)) = 0;
	[model] = addSinkReactions(model,{'trp-L(c)','melatn(c)'},[-0.078431375 -0.078431375; 0 1000]);
	model.c(ismember(model.rxns,'sink_melatn(c)'))=1;
	if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'trp_L(c) -> melatn(c)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% his_L(c) -> hista(c)
	model = modelOri;
	model.c(find(model.c)) = 0;
	[model] = addSinkReactions(model,{'his-L(c)','hista(c)'},[-0.2 -0.2; 0 1000]);
	model.c(ismember(model.rxns,'sink_hista(c)'))=1;
	if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'his_L(c) -> hista(c)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% gly -> co2 + nh4/ glycine cleavage system
	model = modelOri;
	model.c(find(model.c)) = 0;
	[model] = addSinkReactions(model,{'gly(c)','co2(c)','nh4(c)'},[-10 -10; 0.1 1000; 0.1 1000]);
	model.c(ismember(model.rxns,'sink_nh4(c)'))=1;
	if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'gly -> co2 + nh4/ glycine cleavage system';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

	%% BrainTasks_New
    %% glu_L(c) -> gln_L(c)/ glutamine synthase
	model = modelOri;
	model.c(find(model.c)) = 0;
	[model] = addSinkReactions(model,{'glu-L(c)','gln-L(c)'},[-10 -10; 0 1000]);
	model.c(ismember(model.rxns,'sink_gln-L(c)'))=1;
	if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glu_L(c) -> gln_L(c)/ glutamine synthase';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% ATP max, aerobic, pyruvate/ pyruvate dehydrogenase-->TCA->energy
    model = modelOri;
    model.c(find(model.c)) = 0;
    
    model.lb(ismember(model.rxns,'EX_pyr(e)'))=-0.22727273;model.ub(ismember(model.rxns,'EX_pyr(e)'))=-0.22727273;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    
    model.c(ismember(model.rxns,'DM_atp(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ATP max, aerobic, pyruvate/ pyruvate dehydrogenase-->TCA->energy';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% Pyruvate -> Lactate, aerobic
    model = modelOri;
    model.c(find(model.c)) = 0;
    
    model.lb(ismember(model.rxns,'EX_pyr(e)'))=-0.22727273;model.ub(ismember(model.rxns,'EX_pyr(e)'))=-0.22727273;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addSinkReactions(model,{'pyr(m)','lac-L(m)'},[-0.22727273 -0.22727273; 0 1000]);
	model.c(ismember(model.rxns,'sink_lac-L(m)'))=1;
   
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Pyruvate -> Lactate, aerobic';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
        
    %% cys_L(m) -> coa(m)/ CoA synthesis from cysteine
	model = modelOri;
	model.c(find(model.c)) = 0;
	[model] = addSinkReactions(model,{'cys-L(m)','coa(m)'},[-0.2603306 -0.2603306; 0 1000]);
	model.c(ismember(model.rxns,'sink_coa(m)'))=1;
	if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'cys_L(m) -> coa(m)/ CoA synthesis from cysteine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
	
    %% chol[c] -> ach[c]/ acetyl-choline synthesis in brain
	model = modelOri;
	model.c(find(model.c)) = 0;
	[model] = addSinkReactions(model,{'chol(c)','ach(c)'},[-0.028571429 -0.028571429; 0 1000]);
	model.c(ismember(model.rxns,'sink_ach(c)'))=1;
	if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'chol[c] -> ach[c]/ acetyl-choline synthesis in brain';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% pyr[m] -> oaa[m]/ pyruvate carboxylase
	model = modelOri;
	model.c(find(model.c)) = 0;
	[model] = addSinkReactions(model,{'pyr(m)','oaa(m)'},[-0.22727273 -0.22727273; 0 1000]);
	model.c(ismember(model.rxns,'sink_oaa(m)'))=1;
	if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'pyr[m] -> oaa[m]/ pyruvate carboxylase';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% GABA aminotransferase
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_glu-L(e)'))=-10;model.ub(ismember(model.rxns,'EX_glu-L(e)'))=-10;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    model.c(ismember(model.rxns,'ABTArm'))=1;
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'GABA aminotransferase';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% methionine adenosyltransferase
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_met-L(e)'))=-0.20134228;model.ub(ismember(model.rxns,'EX_met-L(e)'))=-0.20134228;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    model.c(ismember(model.rxns,'METAT'))=1;
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'methionine adenosyltransferase';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% creatine synthesis
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_arg-L(e)'))=-0.39810428;model.ub(ismember(model.rxns,'EX_arg-L(e)'))=-0.39810428;
    model.lb(ismember(model.rxns,'EX_gly(e)'))=-10;model.ub(ismember(model.rxns,'EX_gly(e)'))=-10;
    [model] = addSinkReactions(model,{'crtn(c)'},[0 1000]);
    model.c(ismember(model.rxns,'sink_crtn(c)'))=1;
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'creatine synthesis';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% astrocyte function tests:    
    %% alanine demand synapse, aerobic
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'ala-L(s)'});
    model.c(ismember(model.rxns,'DM_ala-L(s)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'alanine demand synapse, aerobic';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    %% lactate demand synapse, aerobic
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'lac-L(s)'});
    model.c(ismember(model.rxns,'DM_lac-L(s)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'lactate demand synapse, aerobic';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    %% citrate demand synapse, aerobic
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'cit(s)'});
    model.c(ismember(model.rxns,'DM_cit(s)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'citrate demand synapse, aerobic';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
     %% glutamine demand synapse, aerobic
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'gln-L(s)'});
    model.c(ismember(model.rxns,'DM_gln-L(s)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glutamine demand synapse, aerobic';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
     %% glutamate synaptic uptake, aerobic
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    
    [model] = addSinkReactions(model,{'glu-L(s)','glu-L(c)'},[-10 -1; 0 1000]);
    model.c(ismember(model.rxns,'sink_glu-L(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glutamate synaptic uptake, aerobic';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% some general tests:
    %% co2 extracellular release, aerobic
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    
    [model] = addDemandReaction(model, {'co2(e)'});
    model.c(ismember(model.rxns,'DM_co2(e)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'co2 extracellular release, aerobic';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
   
    %% co2 intracellular, aerobic
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    
    [model] = addDemandReaction(model, {'co2(c)'});
    model.c(ismember(model.rxns,'DM_co2(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'co2 intracellular, aerobic';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% h2o2 intracellular, aerobic
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    
    [model] = addDemandReaction(model, {'h2o2(c)'});
    model.c(ismember(model.rxns,'DM_h2o2(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'h2o2 intracellular, aerobic';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% h202 intracellular, aerobic, KO Dopa[c] -> DOPAL[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    model.lb(ismember(model.rxns,'42A12BOOX'))=0;model.ub(ismember(model.rxns,'42A12BOOX'))=0;
    
    [model] = addDemandReaction(model, {'h2o2(c)'});
    model.c(ismember(model.rxns,'DM_h2o2(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'h202 intracellular, aerobic, KO Dopa[c] -> DOPAL[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% h202 intracellular, aerobic, Max Dopa[c] -> DOPAL[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    model.lb(ismember(model.rxns,'42A12BOOX'))=5.4804;model.ub(ismember(model.rxns,'42A12BOOX'))=1000;
    
    [model] = addDemandReaction(model, {'h2o2(c)'});
    model.c(ismember(model.rxns,'DM_h2o2(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'h202 intracellular, aerobic, Max Dopa[c] -> DOPAL[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% DOPAL intracellular, aerobic
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    model.lb(ismember(model.rxns,'42A12BOOX'))=0;model.ub(ismember(model.rxns,'42A12BOOX'))=1000;
    
    [model] = addDemandReaction(model, {'34dhpac(c)'});
    model.c(ismember(model.rxns,'DM_34dhpac(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'DOPAL intracellular, aerobic';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% DOPAL intracellular, aerobic, KO Dopa[c] -> DOPAL[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    model.lb(ismember(model.rxns,'42A12BOOX'))=0;model.ub(ismember(model.rxns,'42A12BOOX'))=0;
    
    [model] = addDemandReaction(model, {'34dhpac(c)'});
    model.c(ismember(model.rxns,'DM_34dhpac(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'DOPAL intracellular, aerobic, KO Dopa[c] -> DOPAL[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
       
    %% DOPAL intracellular, aerobic, Max Dopa[c] -> DOPAL[c]
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    model.lb(ismember(model.rxns,'42A12BOOX'))=5.4804;model.ub(ismember(model.rxns,'42A12BOOX'))=1000;
    
    [model] = addDemandReaction(model, {'34dhpac(c)'});
    model.c(ismember(model.rxns,'DM_34dhpac(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'DOPAL intracellular, aerobic, Max Dopa[c] -> DOPAL[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

	%% Mets for below rxns NOT present in media:
    %% BrainTasks_Old
    %% ATP max, aerobic, citrate/TCA cycle
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_cit(e)'))=0;model.ub(ismember(model.rxns,'EX_cit(e)'))=0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    model.c(ismember(model.rxns,'DM_atp(c)'))=1;
  
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ATP max, aerobic, citrate/TCA cycle';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
     %% g6p(c) -> r5p(c)/ HMP shunt
	model = modelOri;
	model.c(find(model.c)) = 0;
	[model] = addSinkReactions(model,{'g6p(c)','r5p(c)'},[0 0; 0 1000]);
	model.c(ismember(model.rxns,'sink_r5p(c)'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'g6p(c) -> r5p(c)/ HMP shunt';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA  
    
    %% accoa(c) -> pmtcoa(c)/ fatty acid synthesis
	model = modelOri;
	model.c(find(model.c)) = 0;
	[model] = addSinkReactions(model,{'accoa(c)','pmtcoa(c)'},[0 0; 0 1000]);
	model.c(ismember(model.rxns,'sink_pmtcoa(c)'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'accoa(c) -> pmtcoa(c)/ fatty acid synthesis';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% adenine -> amp/ salvage of adenine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_adn(e)'))=0;model.ub(ismember(model.rxns,'EX_adn(e)'))=0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    model.c(ismember(model.rxns,'sink_amp(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'adenine -> amp/ salvage of adenine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% hypoxanthine -> imp/ salvage of hypoxanthine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_hxan(e)'))=0;model.ub(ismember(model.rxns,'EX_hxan(e)'))=0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    model.c(ismember(model.rxns,'sink_imp(c)'))=1;
   
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'hypoxanthine -> imp/ salvage of hypoxanthine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% guanine -> gmp/ salvage of guanine
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_gua(e)'))=0;model.ub(ismember(model.rxns,'EX_gua(e)'))=0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    model.c(ismember(model.rxns,'sink_gmp(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'guanine -> gmp/ salvage of guanine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% arachd(r) -> txa2(r)
	model = modelOri;
	model.c(find(model.c)) = 0;
	[model] = addSinkReactions(model,{'arachd(r)','txa2(r)'},[0 0; 0 1000]);
	model.c(ismember(model.rxns,'sink_txa2(r)'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'arachd(r) -> txa2(r)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA  
    
    %% BrainTasks_New
    %% bhb(m) -> acac(m)/ ketone body utilization
	model = modelOri;
	model.c(find(model.c)) = 0;
	[model] = addSinkReactions(model,{'bhb(m)','acac(m)'},[0 0; 0 1000]);
	model.c(ismember(model.rxns,'sink_acac(m)'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'bhb(m) -> acac(m)/ ketone body utilization';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% mal_L(m) -> pyr(m)/ malic enzyme
	model = modelOri;
	model.c(find(model.c)) = 0;
	[model] = addSinkReactions(model,{'mal-L(m)','pyr(m)'},[0 0; 0 1000]);
	model.c(ismember(model.rxns,'sink_pyr(m)'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'mal_L(m) -> pyr(m)/ malic enzyme';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    
    %% occoa[m] -> accoa[m]/ octanoate oxidation
	model = modelOri;
	model.c(find(model.c)) = 0;
	[model] = addSinkReactions(model,{'occoa(m)','accoa(m)'},[0 0; 0 1000]);
	model.c(ismember(model.rxns,'sink_accoa(m)'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'occoa[m] -> accoa[m]/ octanoate oxidation';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% lnlncgcoa[c] -> dlnlcgcoa[c]/ fatty acid elongation
	model = modelOri;
	model.c(find(model.c)) = 0;
	[model] = addSinkReactions(model,{'lnlncgcoa(c)','dlnlcgcoa(c)'},[0 0; 0 1000]);
	model.c(ismember(model.rxns,'sink_dlnlcgcoa(c)'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'lnlncgcoa[c] -> dlnlcgcoa[c]/ fatty acid elongation';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    

    %% phosphocreatine -> creatine/ cytosolic creatine kinase
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_pcreat(e)'))=0;model.ub(ismember(model.rxns,'EX_pcreat(e)'))=0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    model.c(ismember(model.rxns,'sink_creat(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'phosphocreatine -> creatine/ cytosolic creatine kinase';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% creatine -> phosphocreatine/mitochondrial creatine kinase
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_creat(e)'))=0;model.ub(ismember(model.rxns,'EX_creat(e)'))=0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    model.c(ismember(model.rxns,'sink_pcreat(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'creatine -> phosphocreatine/mitochondrial creatine kinase';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    
    %% arachd(c) -> prostgh2(c)/ prostaglandin synthesis
	model = modelOri;
	model.c(find(model.c)) = 0;
	[model] = addSinkReactions(model,{'arachd(c)','prostgh2(c)'},[0 0; 0 1000]);
	model.c(ismember(model.rxns,'sink_prostgh2(c)'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'arachd(c) -> prostgh2(c)/ prostaglandin synthesis';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% arachd(c) -> prostgd2(r)/ prostaglandin synthesis
	model = modelOri;
	model.c(find(model.c)) = 0;
	[model] = addSinkReactions(model,{'arachd(c)','prostgd2(r)'},[0 0; 0 1000]);
	model.c(ismember(model.rxns,'sink_prostgd2(r)'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'arachd(c) -> prostgd2(r)/ prostaglandin synthesis';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% arachd(c) -> prostge2(r)/ prostaglandin synthesis
	model = modelOri;
	model.c(find(model.c)) = 0;
	[model] = addSinkReactions(model,{'arachd(c)','prostge2(r)'},[0 0; 0 1000]);
	model.c(ismember(model.rxns,'sink_prostge2(r)'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'arachd(c) -> prostge2(r)/ prostaglandin synthesis';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% arachd(c) -> prostgi2(r)/ prostaglandin synthesis
	model = modelOri;
	model.c(find(model.c)) = 0;
	[model] = addSinkReactions(model,{'arachd(c)','prostgi2(r)'},[0 0; 0 1000]);
	model.c(ismember(model.rxns,'sink_prostgi2(r)'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'arachd(c) -> prostgi2(r)/ prostaglandin synthesis';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% arachd(c) -> leuktrE4(c)/ leukotriene synthesis
	model = modelOri;
	model.c(find(model.c)) = 0;
	[model] = addSinkReactions(model,{'arachd(c)','leuktrE4(c)'},[0 0; 0 1000]);
	model.c(ismember(model.rxns,'sink_leuktrE4(c)'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'arachd(c) -> leuktrE4(c)/ leukotriene synthesis';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% arachd(c) -> C06314(c)/ lipoxin synthesis
	model = modelOri;
	model.c(find(model.c)) = 0;
	[model] = addSinkReactions(model,{'arachd(c)','C06314(c)'},[0 0; 0 1000]);
	model.c(ismember(model.rxns,'sink_C06314(c)'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'arachd(c) -> C06314(c)/ lipoxin synthesis';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% nrpphr(c) -> 3mox4hoxm(c)/ degradation of norepinephrine
	model = modelOri;
	model.c(find(model.c)) = 0;
	[model] = addSinkReactions(model,{'nrpphr(c)','3mox4hoxm(c)'},[0 0; 0 1000]);
	model.c(ismember(model.rxns,'sink_3mox4hoxm(c)'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'nrpphr(c) -> 3mox4hoxm(c)/ degradation of norepinephrine';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
 
          
end

TestSolutionName(:,2) = num2cell(TestSolution);
TestedRxns = unique(TestedRxns);
TestedRxns = intersect(modelOri.rxns,TestedRxns); % only those reactions that are also in modelOri not those that have been added to the network
PercTestedRxns = length(TestedRxns)*100/length(modelOri.rxns);


