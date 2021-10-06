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
% %% define closed model (as per Sanity check!)
% modelClosed = model;
% modelexchanges1 = strmatch('Ex_',modelClosed.rxns);
% modelexchanges4 = strmatch('EX_',modelClosed.rxns);
% % modelexchanges2 = strmatch('DM_',modelClosed.rxns);
% % modelexchanges3 = strmatch('sink_',modelClosed.rxns);
% % selExc = (find( full((sum(abs(modelClosed.S)==1,1) ==1) & (sum(modelClosed.S~=0) == 1))))';
% 
% % set bounds of all exchanges to: lb=0 & ub=1000
% modelexchanges = unique([modelexchanges1; modelexchanges4;]);
% modelClosed.lb(find(ismember(modelClosed.rxns,modelClosed.rxns(modelexchanges))))=0;
% % modelClosed.ub(find(ismember(modelClosed.rxns,modelClosed.rxns(modelexchanges))))=1000;
% model = modelClosed;

% ASM_AC_END
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
%model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=0;

model.c(find(model.c)) = 0;
modelOri = model;
k = 1;
%RPMI_composition={'EX_ala_L(e)','EX_arg-L(e)','EX_asn_L(e)','EX_asp_L(e)','EX_cys-L(e)','EX_gln-L(e)','EX_glu-L(e)','EX_gly(e)','EX_his-L(e)','EX_ile_L(e)','EX_leu_L(e)','EX_lys-L(e)','EX_met_L(e)','EX_phe_L(e)','EX_4HPRO','EX_pro-L(e)','EX_ser_L(e)','EX_thr_L(e)','EX_trp_L(e)','EX_tyr_L(e)','EX_val_L(e)','EX_ascb_L(e)','EX_btn(e)','EX_chol(e)','EX_pnto_R(e)','EX_fol(e)','EX_ncam(e)','EX_pydxn(e)','EX_ribflv(e)','EX_thm(e)','EX_cbl1(e)','EX_inost(e)','EX_ca2(e)','EX_fe3(e)','EX_k(e)','EX_hco3(e)','EX_na1(e)','EX_pi(e)','EX_glc(e)','EX_hxan(e)','EX_lnlc(e)','EX_lipoate(e)','EX_ptrc(e)','EX_pyr(e)','EX_thymd(e)','EX_etha(e)','EX_gthrd(e)'};

 if strcmp(test,'Recon1') || strcmp(test,'all') || strcmp(test,'Harvey')
%     modelOri = modelOri;
%     mediumCompounds = {'EX_ala_L(e)', 'EX_gln-L(e)', 'EX_arg-L(e)', 'EX_asn_L(e)', 'EX_asp_L(e)', 'EX_ca2(e)', 'EX_chol(e)', 'EX_cl(e)', 'EX_cys-L(e)', 'EX_fe3(e)', 'EX_fol(e)', 'EX_glc(e)', 'EX_glu-L(e)', 'EX_gly(e)', 'EX_hco3(e)', 'EX_his-L(e)', 'EX_ile-L(e)', 'EX_inost(e)', 'EX_k(e)', 'EX_leu_L(e)', 'EX_lys-L(e)', 'EX_met_L(e)', 'EX_mg2(e)', 'EX_na1(e)', 'EX_ncam(e)', 'EX_no3(e)', 'EX_phe_L(e)', 'EX_pi(e)', 'EX_pnto-R(e)', 'EX_prgstrn(e)', 'EX_pro-L(e)', 'EX_ptrc(e)', 'EX_pydx(e)', 'EX_pyr(e)', 'EX_ribflv(e)', 'EX_selni(e)', 'EX_ser-L(e)', 'EX_so4(e)', 'EX_thm(e)', 'EX_thr-L(e)', 'EX_trp_L(e)', 'EX_tyr-L(e)', 'EX_val_L(e)', 'EX_zn2(e)','sink_chsterol(c)','sink_clpn-hs(c)','sink_ctp(c)','sink_gtp(c)','sink_ile-L(c)','sink_pail-hs(c)','sink_pe-hs(c)','sink_pglyc-hs(c)','sink_ps-hs(c)','sink_sphmyln-hs(c)','sink_utp(c)'}; 
%  
%     I = strmatch('EX_', modelOri.rxns);
%     I = strmatch('sink_', modelOri.rxns);
%     I = strmatch('DM_', modelOri.rxns);
%     
%     for i=1:length(I);
%         Ex= I(i);
%         modelOri.lb(Ex,1) = 0;
%         if modelOri.ub(Ex,1) < 0;
%             modelOri.ub(Ex,1)=1
%         end
%         modelOri.ub(Ex,1) = 1000;% uncomment to run for tcell models
%     end
%      modelOri.lb(find(ismember(modelOri.rxns,mediumCompounds)))=-10;
   
	%% ASTROCYTE_FUNCTION_TESTS
    %% "Human Recon3D biomass_maintenance"
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'biomass_maintenance'))=1;
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Human Recon3D biomass_maintenance';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% "Ca ATPase"
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'CAATPS'))=1;
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Ca ATPase';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Aspartate Transaminase / asp_L -> oaa
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'ASPTA'))=1;
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Aspartate Transaminase / asp_L -> oaa';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Aspartate Transaminase / asp_L -> oaa (mito)
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'ASPTAm'))=1;
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Aspartate Transaminase / asp_L -> oaa (mito)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
        
    %% L-Alanine Transaminase / ala_L -> pyr
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'ALATA_L'))=1;
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'L-Alanine Transaminase / ala_L -> pyr';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA    

	%% Glutamate Dehydrogenase (NAD), Mitochondrial / glu_L -> nh4
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'GLUDxm'))=1;
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Glutamate Dehydrogenase (NAD), Mitochondrial / glu_L -> nh4';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% Glutaminase / gln_L -> nh4
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'HMR_9802'))=1;
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Glutaminase / gln_L -> nh4';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Glutaminase / gln_L -> nh4 (mito)
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'GLUNm'))=1;
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Glutaminase / gln_L -> nh4 (mito)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
        
    %% Glutamate Decarboxylase / glu_L -> 4abut
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'GLUDC'))=1;
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero')
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Glutamate Decarboxylase / glu_L -> 4abut';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA    

    %% Glutathione Synthetase / cys_L + glu_L + gly -> ghtrd
	model = modelOri;
	model.c(find(model.c)) = 0;
	model.c(ismember(model.rxns,'GTHS'))=1;
	if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Glutathione Synthetase / cys_L + glu_L + gly -> ghtrd';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA   
        
    %% Histidine Decarboxylase / his_L -> hista
	model = modelOri;
	model.c(find(model.c)) = 0;
	model.c(ismember(model.rxns,'HISDC'))=1;
	if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Histidine Decarboxylase / his_L -> hista';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% Glutamine Synthetase / glu_L -> gln_L
	model = modelOri;
	model.c(find(model.c)) = 0;
	model.c(ismember(model.rxns,'GLNS'))=1;
	if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Glutamine Synthetase / glu_L -> gln_L';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% Pyruvate Dehydrogenase / pyr -> accoa
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'PDHm'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Pyruvate Dehydrogenase / pyr -> accoa';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% L-Lactate Dehydrogenase / pyr -> lac_L
    model = modelOri;
    model.c(find(model.c)) = 0;
	model.c(ismember(model.rxns,'LDH_L'))=1;
   
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'L-Lactate Dehydrogenase / pyr -> lac_L';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% L-Lactate Dehydrogenase / pyr -> lac_L (mito)
    model = modelOri;
    model.c(find(model.c)) = 0;
	model.c(ismember(model.rxns,'LDH_Lm'))=1;
   
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'L-Lactate Dehydrogenase / pyr -> lac_L (mito)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% Pyruvate Carboxylase / pyr -> oaa
	model = modelOri;
	model.c(find(model.c)) = 0;
	model.c(ismember(model.rxns,'PCm'))=1;
	if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Pyruvate Carboxylase / pyr -> oaa';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
        	        	        
    %% Acetyl-choline synthesis in brain / chol -> ach
	model = modelOri;
	model.c(find(model.c)) = 0;
	[model] = addSinkReactions(model,{'chol(c)','ach(c)'},[-10 -0.01; 0 1000]);
	model.c(ismember(model.rxns,'sink_ach(c)'))=1;
	if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Acetyl-choline synthesis in brain / chol -> ach';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
        
    %% 4-Aminobutyrate Transaminase, Reversible, Mitochondrial / 4abut -> glu_L
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'ABTArm'))=1;
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = '4-Aminobutyrate Transaminase, Reversible, Mitochondrial / 4abut -> glu_L';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% Methionine Adenosyltransferase / met_L -> amet
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'METAT'))=1;
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Methionine Adenosyltransferase / met_L -> amet';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% ASTROCYTE-SYNAPSE EXCHANGES (START)    
    %% Exchange of L-Alanine, Synapse
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'EX_ala_L(s)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Exchange of L-Alanine, Synapse';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% Exchange of L-Lactate, Synapse
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'EX_lac_L(s)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Exchange of L-Lactate, Synapse';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% Exchange of Citrate, Synapse
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'EX_cit(s)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Exchange of Citrate, Synapse';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% Exchange of L-Glutamine, Synapse
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'EX_gln_L(s)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Exchange of L-Glutamine, Synapse';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% Exchange of L-Glutamate, Synapse
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'EX_glu_L(s)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Exchange of L-Glutamate, Synapse';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% Exchange of Glutathionie, Synapse
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'EX_gthrd(s)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Exchange of Glutathionie, Synapse';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% Exchange of Potassium, Synapse
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'EX_k(s)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Exchange of Potassium, Synapse';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Exchange of Cholesterol, Synapse
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'EX_chsterol(s)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Exchange of Cholesterol, Synapse';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Exchange of Hydrogen sulfide, Synapse
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'EX_HC00250(s)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Exchange of Hydrogen sulfide, Synapse';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Exchange of Ascorbic acid, Synapse
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'EX_ascb_L(s)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Exchange of Ascorbic acid, Synapse';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Exchange of ATP, Synapse
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'EX_atp(s)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Exchange of ATP, Synapse';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% ASTROCYTE-SYNAPSE EXCHANGES (END)    
    %% Adenosine Kinase / adn -> amp
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'ADNK1'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Adenosine Kinase / adn -> amp';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Adenosine Kinase / adn -> amp (mito)
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'ADNK1m'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Adenosine Kinase / adn -> amp (mito)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
        
    %% Hypoxanthine Phosphoribosyltransferase (Hypoxanthine) / hxan -> imp
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'HXPRT'))=1;
   
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Hypoxanthine Phosphoribosyltransferase (Hypoxanthine) / hxan -> imp';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% Guanine Phosphoribosyltransferase / gua -> gmp
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'GUAPRT'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Guanine Phosphoribosyltransferase / gua -> gmp';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
        
    %% (R)-3-Hydroxybutanoate:NAD+ Oxidoreductase / bhb -> acac
	model = modelOri;
	model.c(find(model.c)) = 0;
	model.c(ismember(model.rxns,'BDHm'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = '(R)-3-Hydroxybutanoate:NAD+ Oxidoreductase / bhb -> acac';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% Malic Enzyme (NAD), Mitochondrial / mal_L -> pyr (mito)
	model = modelOri;
	model.c(find(model.c)) = 0;
	model.c(ismember(model.rxns,'ME1m'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Malic Enzyme (NAD), Mitochondrial / mal_L -> pyr (mito)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Malic Enzyme (NADP) / mal_L -> pyr
	model = modelOri;
	model.c(find(model.c)) = 0;
	model.c(ismember(model.rxns,'ME2'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Malic Enzyme (NADP) / mal_L -> pyr';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
            
    %% ATP Creatine Kinase / creat -> pcreat (mito)
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'CK'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ATP Creatine Kinase / creat -> pcreat (mito)';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% ATP Creatine Kinase, Cytosolic / creat -> pcreat
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.c(ismember(model.rxns,'CKc'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ATP Creatine Kinase, Cytosolic / creat -> pcreat';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
           
    %% Prostaglandin G/H Synthase / arachd -> prostgh2
	model = modelOri;
	model.c(find(model.c)) = 0;
	model.c(ismember(model.rxns,'PGS'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Prostaglandin G/H Synthase / arachd -> prostgh2';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% 3-Sulfino-L-Alanine Carboxy-Lyase / 3sala -> hyptaur
	model = modelOri;
	model.c(find(model.c)) = 0;
	model.c(ismember(model.rxns,'3SALACBOXL'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = '3-Sulfino-L-Alanine Carboxy-Lyase / 3sala -> hyptaur';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Hypotaurine Oxidase / hyptaur -> taur
	model = modelOri;
	model.c(find(model.c)) = 0;
	model.c(ismember(model.rxns,'HYPTROX'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Hypotaurine Oxidase / hyptaur -> taur';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% Exchange of Taurine / taur[e] <=> 
	model = modelOri;
	model.c(find(model.c)) = 0;
	model.c(ismember(model.rxns,'EX_taur(e)'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Exchange of Taurine / taur[e] <=> ';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Serine racemase / ser_L <=> ser_D
	model = modelOri;
	model.c(find(model.c)) = 0;
	model.c(ismember(model.rxns,'SRR'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Serine racemase / ser_L <=> ser_D';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% PI-Cycle Reactions (N=14)        
    %% Myo-Inositol-1-Phosphate Synthase / g6p -> mi1p_D
	model = modelOri;
	model.c(find(model.c)) = 0;
	model.c(ismember(model.rxns,'MI1PS'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Myo-Inositol-1-Phosphate Synthase / g6p -> mi1p_D';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% Myo-Inositol 1-Phosphatase / mi1p_D -> inost
	model = modelOri;
	model.c(find(model.c)) = 0;
	model.c(ismember(model.rxns,'MI1PP'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Myo-Inositol 1-Phosphatase / mi1p_D -> inost';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% Transport of Inositol via Sodium Symport / inost[e] <=> inost[c]
	model = modelOri;
	model.c(find(model.c)) = 0;
	model.c(ismember(model.rxns,'INSTt4'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Transport of Inositol via Sodium Symport / inost[e] <=> inost[c]';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% Phosphatidylinositol Synthase / inost <=> pail_hs
	model = modelOri;
	model.c(find(model.c)) = 0;
	model.c(ismember(model.rxns,'CDIPTr'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Phosphatidylinositol Synthase / inost <=> pail_hs';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% Phosphatidylinositol 4-Kinase / pail_hs -> pail4p_hs
	model = modelOri;
	model.c(find(model.c)) = 0;
	model.c(ismember(model.rxns,'PIK4'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Phosphatidylinositol 4-Kinase / pail_hs -> pail4p_hs';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% Phosphatidylinositol 4-Phosphate 5-Kinase / pail4p_hs -> pail45p_hs
	model = modelOri;
	model.c(find(model.c)) = 0;
	model.c(ismember(model.rxns,'PI4P5K'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Phosphatidylinositol 4-Phosphate 5-Kinase / pail4p_hs -> pail45p_hs';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% Phosphatidylinositol 4, 5-Bisphosphate Phospholipase C / pail45p_hs -> dag_hs + mi145p
	model = modelOri;
	model.c(find(model.c)) = 0;
	model.c(ismember(model.rxns,'PI45PLC'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Phosphatidylinositol 4, 5-Bisphosphate Phospholipase C / pail45p_hs -> dag_hs + mi145p';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% Inositol-1, 4, 5-Trisphosphate 5-Phosphatase / mi145p -> mi14p
	model = modelOri;
	model.c(find(model.c)) = 0;
	model.c(ismember(model.rxns,'MI145PP'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Inositol-1, 4, 5-Trisphosphate 5-Phosphatase / mi145p -> mi14p';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% Inositol-1, 4-Bisphosphate 4-Phosphatase / mi14p -> mi1p_D
	model = modelOri;
	model.c(find(model.c)) = 0;
	model.c(ismember(model.rxns,'MI14P4P'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Inositol-1, 4-Bisphosphate 4-Phosphatase / mi14p -> mi1p_D';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% Diacylglycerol Phosphate Kinase / dag_hs <=> pa_hs
	model = modelOri;
	model.c(find(model.c)) = 0;
	model.c(ismember(model.rxns,'DAGK_hs'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Diacylglycerol Phosphate Kinase / dag_hs <=> pa_hs';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% Phosphatidate Cytidylyltransferase / pa_hs -> cdpdag_hs
	model = modelOri;
	model.c(find(model.c)) = 0;
	model.c(ismember(model.rxns,'CDS'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Phosphatidate Cytidylyltransferase / pa_hs -> cdpdag_hs';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% Inositol-1, 4, 5-Trisphosphate 3-Kinase / mi145p[c] -> mi1345p
	model = modelOri;
	model.c(find(model.c)) = 0;
	model.c(ismember(model.rxns,'MI145PK'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Inositol-1, 4, 5-Trisphosphate 3-Kinase / mi145p[c] -> mi1345p';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% Inositol-1, 3, 4, 5-Triphosphate 6-Kinase / mi1345p -> mi13456p
	model = modelOri;
	model.c(find(model.c)) = 0;
	model.c(ismember(model.rxns,'MI1345PKc'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Inositol-1, 3, 4, 5-Triphosphate 6-Kinase / mi1345p -> mi13456p';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
    %% Inositol-1, 3, 4, 5, 6-Pentakisphosphate 2-Kinase / mi13456p -> minohp
	model = modelOri;
	model.c(find(model.c)) = 0;
	model.c(ismember(model.rxns,'MI13456PK'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'min','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Inositol-1, 3, 4, 5, 6-Pentakisphosphate 2-Kinase / mi13456p -> minohp';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

end

TestSolutionName(:,2) = num2cell(TestSolution);
TestedRxns = unique(TestedRxns);
TestedRxns = intersect(modelOri.rxns,TestedRxns); % only those reactions that are also in modelOri not those that have been added to the network
PercTestedRxns = length(TestedRxns)*100/length(modelOri.rxns);
