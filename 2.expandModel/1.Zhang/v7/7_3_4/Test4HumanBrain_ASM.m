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
   
	%% Generic
    %% "Human Recon3D biomass_maintenance"
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    model.c(ismember(model.rxns,'biomass_maintenance'))=1;
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Human Recon3D biomass_maintenance';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Mets for below rxns are present in media:
    %% BrainTasks_Old
    %% asp-L -> oaa / transamination pair
    model = modelOri;
    model.c(find(model.c)) = 0;
    [model] = addSinkReactions(model,{'asp-L(c)','oaa(c)'},[-10 -0.01; 0 1000]);
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
    [model] = addSinkReactions(model,{'ala-L(c)','pyr(c)'},[-10 -0.01; 0 1000]);
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
    model.lb(ismember(model.rxns,'EX_glu-L(e)'))=-10;model.ub(ismember(model.rxns,'EX_glu-L(e)'))=-0.01;
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
    [model] = addSinkReactions(model,{'gln-L(c)','nh4(c)'},[-10 -0.01; 0 1000]);
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
    [model] = addSinkReactions(model,{'glu-L(c)','4abut(c)'},[-10 -0.01; 0 1000]);
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
	[model] = addSinkReactions(model,{'cys-L(c)','glu-L(c)','gly(c)','gthrd(c)'},[-10 -0.01; -10 -0.01; -10 -0.01; 0 1000]);
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
    model.lb(ismember(model.rxns,'EX_glc(e)'))=-10;model.ub(ismember(model.rxns,'EX_glc(e)'))=-0.01;
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
	[model] = addSinkReactions(model,{'val-L(c)','succoa(m)'},[-10 -0.01; 0 1000]);
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
	[model] = addSinkReactions(model,{'leu-L(c)','accoa(c)'},[-10 -0.01; 0 1000]);
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
	[model] = addSinkReactions(model,{'ile-L(c)','accoa(c)'},[-10 -0.01; 0 1000]);
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
	[model] = addSinkReactions(model,{'tyr-L(c)','nrpphr(c)'},[-10 -0.01; 0 1000]);
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
	[model] = addSinkReactions(model,{'tyr-L(c)','adrnl(c)'},[-10 -0.01; 0 1000]);
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
	[model] = addSinkReactions(model,{'trp-L(c)','melatn(c)'},[-10 -0.01; 0 1000]);
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
	[model] = addSinkReactions(model,{'his-L(c)','hista(c)'},[-10 -0.01; 0 1000]);
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
	[model] = addSinkReactions(model,{'gly(c)','co2(c)','nh4(c)'},[-10 -0.01; 0.1 1000; 0.1 1000]);
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
	[model] = addSinkReactions(model,{'glu-L(c)','gln-L(c)'},[-10 -0.01; 0 1000]);
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
    
    model.lb(ismember(model.rxns,'EX_pyr(e)'))=-10;model.ub(ismember(model.rxns,'EX_pyr(e)'))=-0.01;
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
    
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addSinkReactions(model,{'pyr(m)','lac-L(m)'},[-10 -0.01; 0 1000]);
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
	[model] = addSinkReactions(model,{'cys-L(m)','coa(m)'},[-10 -0.01; 0 1000]);
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
	[model] = addSinkReactions(model,{'chol(c)','ach(c)'},[-10 -0.01; 0 1000]);
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
	[model] = addSinkReactions(model,{'pyr(m)','oaa(m)'},[-10 -0.01; 0 1000]);
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
    model.lb(ismember(model.rxns,'EX_glu-L(e)'))=-10;model.ub(ismember(model.rxns,'EX_glu-L(e)'))=-0.01;
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
    model.lb(ismember(model.rxns,'EX_met-L(e)'))=-10;model.ub(ismember(model.rxns,'EX_met-L(e)'))=-0.01;
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
    model.lb(ismember(model.rxns,'EX_arg-L(e)'))=-10;model.ub(ismember(model.rxns,'EX_arg-L(e)'))=-0.01;
    model.lb(ismember(model.rxns,'EX_gly(e)'))=-10;model.ub(ismember(model.rxns,'EX_gly(e)'))=-0.01;
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
    
    [model] = addSinkReactions(model,{'glu-L(s)','glu-L(c)'},[-10 -0.01; 0 1000]);
    model.c(ismember(model.rxns,'sink_glu-L(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glutamate synaptic uptake, aerobic';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
     %% glutathione demand synapse, aerobic
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'gthrd(s)'});
    model.c(ismember(model.rxns,'DM_gthrd(s)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'glutathione demand synapse, aerobic';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
     %% Potassium synaptic uptake, aerobic
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    
    [model] = addSinkReactions(model,{'k(s)','k(c)'},[-10 -0.01; 0 1000]);
    model.c(ismember(model.rxns,'sink_k(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Potassium synaptic uptake, aerobic';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
    
	%% Mets for below rxns NOT present in media:
    %% BrainTasks_Old
    %% ATP max, aerobic, citrate/TCA cycle
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_cit(e)'))=-10;model.ub(ismember(model.rxns,'EX_cit(e)'))=-0.01;
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
	[model] = addSinkReactions(model,{'g6p(c)','r5p(c)'},[-10 -0.01; 0 1000]);
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
	[model] = addSinkReactions(model,{'accoa(c)','pmtcoa(c)'},[-10 -0.01; 0 1000]);
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
    model.lb(ismember(model.rxns,'EX_adn(e)'))=-10;model.ub(ismember(model.rxns,'EX_adn(e)'))=-0.01;
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
    model.lb(ismember(model.rxns,'EX_hxan(e)'))=-10;model.ub(ismember(model.rxns,'EX_hxan(e)'))=-0.01;
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
    model.lb(ismember(model.rxns,'EX_gua(e)'))=-10;model.ub(ismember(model.rxns,'EX_gua(e)'))=-0.01;
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
	[model] = addSinkReactions(model,{'arachd(r)','txa2(r)'},[-10 -0.01; 0 1000]);
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
	[model] = addSinkReactions(model,{'bhb(m)','acac(m)'},[-10 -0.01; 0 1000]);
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
	[model] = addSinkReactions(model,{'mal-L(m)','pyr(m)'},[-10 -0.01; 0 1000]);
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
	[model] = addSinkReactions(model,{'occoa(m)','accoa(m)'},[-10 -0.01; 0 1000]);
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
	[model] = addSinkReactions(model,{'lnlncgcoa(c)','dlnlcgcoa(c)'},[-10 -0.01; 0 1000]);
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
    model.lb(ismember(model.rxns,'EX_pcreat(e)'))=-10;model.ub(ismember(model.rxns,'EX_pcreat(e)'))=-0.01;
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
    model.lb(ismember(model.rxns,'EX_creat(e)'))=-10;model.ub(ismember(model.rxns,'EX_creat(e)'))=-0.01;
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
	[model] = addSinkReactions(model,{'arachd(c)','prostgh2(c)'},[-10 -0.01; 0 1000]);
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
	[model] = addSinkReactions(model,{'arachd(c)','prostgd2(r)'},[-10 -0.01; 0 1000]);
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
	[model] = addSinkReactions(model,{'arachd(c)','prostge2(r)'},[-10 -0.01; 0 1000]);
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
	[model] = addSinkReactions(model,{'arachd(c)','prostgi2(r)'},[-10 -0.01; 0 1000]);
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
	[model] = addSinkReactions(model,{'arachd(c)','leuktrE4(c)'},[-10 -0.01; 0 1000]);
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
	[model] = addSinkReactions(model,{'arachd(c)','C06314(c)'},[-10 -0.01; 0 1000]);
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
	[model] = addSinkReactions(model,{'nrpphr(c)','3mox4hoxm(c)'},[-10 -0.01; 0 1000]);
	model.c(ismember(model.rxns,'sink_3mox4hoxm(c)'))=1;
	
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'nrpphr(c) -> 3mox4hoxm(c)/ degradation of norepinephrine';
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
          
end

TestSolutionName(:,2) = num2cell(TestSolution);
TestedRxns = unique(TestedRxns);
TestedRxns = intersect(modelOri.rxns,TestedRxns); % only those reactions that are also in modelOri not those that have been added to the network
PercTestedRxns = length(TestedRxns)*100/length(modelOri.rxns);


