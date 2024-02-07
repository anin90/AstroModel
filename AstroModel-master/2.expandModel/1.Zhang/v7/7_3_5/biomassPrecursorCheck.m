function [TestSolution,TestSolutionName,TestedRxns,PercTestedRxns] = test4HumanFctExt(model,test,optionSinks)
% tests the synthesis of all biomass_precursor metabolites
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
%       - ASC, 31/03/20
%% ASM_AC_START
%% define closed model (as per Sanity check!)
modelClosed = model;
modelexchanges1 = strmatch('Ex_',modelClosed.rxns);
modelexchanges4 = strmatch('EX_',modelClosed.rxns);
% modelexchanges2 = strmatch('DM_',modelClosed.rxns);
% modelexchanges3 = strmatch('sink_',modelClosed.rxns);
% selExc = (find( full((sum(abs(modelClosed.S)==1,1) ==1) & (sum(modelClosed.S~=0) == 1))))';

% set bounds of all exchanges to: lb=0 & ub=1000
modelexchanges = unique([modelexchanges1; modelexchanges4;]);
modelClosed.lb(find(ismember(modelClosed.rxns,modelClosed.rxns(modelexchanges))))=0;
% modelClosed.ub(find(ismember(modelClosed.rxns,modelClosed.rxns(modelexchanges))))=1000;
model = modelClosed;
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
    modelOri = modelOri;
    mediumCompounds = {'EX_ala_L(e)', 'EX_gln-L(e)', 'EX_arg-L(e)', 'EX_asn_L(e)', 'EX_asp_L(e)', 'EX_ca2(e)', 'EX_chol(e)', 'EX_cl(e)', 'EX_cys-L(e)', 'EX_fe3(e)', 'EX_fol(e)', 'EX_glc(e)', 'EX_glu-L(e)', 'EX_gly(e)', 'EX_hco3(e)', 'EX_his-L(e)', 'EX_ile-L(e)', 'EX_inost(e)', 'EX_k(e)', 'EX_leu_L(e)', 'EX_lys-L(e)', 'EX_met_L(e)', 'EX_mg2(e)', 'EX_na1(e)', 'EX_ncam(e)', 'EX_no3(e)', 'EX_phe_L(e)', 'EX_pi(e)', 'EX_pnto-R(e)', 'EX_prgstrn(e)', 'EX_pro-L(e)', 'EX_ptrc(e)', 'EX_pydx(e)', 'EX_pyr(e)', 'EX_ribflv(e)', 'EX_selni(e)', 'EX_ser-L(e)', 'EX_so4(e)', 'EX_thm(e)', 'EX_thr-L(e)', 'EX_trp_L(e)', 'EX_tyr-L(e)', 'EX_val_L(e)', 'EX_zn2(e)','EX_o2(e)','sink_chsterol(c)','sink_clpn-hs(c)','sink_ctp(c)','sink_gtp(c)','sink_ile-L(c)','sink_pail-hs(c)','sink_pe-hs(c)','sink_pglyc-hs(c)','sink_ps-hs(c)','sink_sphmyln-hs(c)','sink_utp(c)'}; 
 
    I = strmatch('EX_', modelOri.rxns);
    I = strmatch('sink_', modelOri.rxns);
    I = strmatch('DM_', modelOri.rxns);
    
    for i=1:length(I);
        Ex= I(i);
        modelOri.lb(Ex,1) = 0;
        if modelOri.ub(Ex,1) < 0;
            modelOri.ub(Ex,1)=1
        end
        modelOri.ub(Ex,1) = 1;% uncomment to run for tcell models
    end
     modelOri.lb(find(ismember(modelOri.rxns,mediumCompounds)))=-10;
    
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

    %% Test for synthesis of biomass_maintenance precursor metabolites:
    %% Alanine demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'ala-L(c)'});
    model.c(ismember(model.rxns,'DM_ala-L(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Alanine demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Arginine demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'arg-L(c)'});
    model.c(ismember(model.rxns,'DM_arg-L(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Arginine demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Asparagine demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'asn-L(c)'});
    model.c(ismember(model.rxns,'DM_asn-L(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Asparagine demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Aspartic acid demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'asp-L(c)'});
    model.c(ismember(model.rxns,'DM_asp-L(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Aspartic acid demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% ATPdemand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    % [model] = addDemandReaction(model, {'atp(c)'});
    model.c(ismember(model.rxns,'DM_atp(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'ATP demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Cholesterol demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'chsterol(c)'});
    model.c(ismember(model.rxns,'DM_chsterol(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Cholesterol demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Cardiolipin demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'clpn-hs(c)'});
    model.c(ismember(model.rxns,'DM_clpn-hs(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Cardiolipin demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% CTP demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'ctp(c)'});
    model.c(ismember(model.rxns,'DM_ctp(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'CTP demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Cysteine demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'cys-L(c)'});
    model.c(ismember(model.rxns,'DM_cys-L(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Cysteine demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Glucose 6-phosphate demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'g6p(c)'});
    model.c(ismember(model.rxns,'DM_g6p(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Glucose 6-phosphate demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Glutamine demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'gln-L(c)'});
    model.c(ismember(model.rxns,'DM_gln-L(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Glutamine demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Glutamate demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'glu-L(c)'});
    model.c(ismember(model.rxns,'DM_glu-L(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Glutamate demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Glycine demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'gly(c)'});
    model.c(ismember(model.rxns,'DM_gly(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Glycine demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% GTP demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'gtp(c)'});
    model.c(ismember(model.rxns,'DM_gtp(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'GTP demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% H2O demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'h2o(c)'});
    model.c(ismember(model.rxns,'DM_h2o(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'H2O demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Histidine demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'his-L(c)'});
    model.c(ismember(model.rxns,'DM_his-L(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Histidine demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Isoleucine demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'ile-L(c)'});
    model.c(ismember(model.rxns,'DM_ile-L(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Isoleucine demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Leucine demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'leu-L(c)'});
    model.c(ismember(model.rxns,'DM_leu-L(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Leucine demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Lysine demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'lys-L(c)'});
    model.c(ismember(model.rxns,'DM_lys-L(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Lysine demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Methionine demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'met-L(c)'});
    model.c(ismember(model.rxns,'DM_met-L(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Methionine demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% 1-Phosphatidyl-1D-Myo-Inositol demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'pail_hs(c)'});
    model.c(ismember(model.rxns,'DM_pail_hs(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = '1-Phosphatidyl-1D-Myo-Inositol demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Phosphatidylcholine demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'pchol-hs(c)'});
    model.c(ismember(model.rxns,'DM_pchol-hs(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Phosphatidylcholine demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Phosphatidylethanolamine demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'pe_hs(c)'});
    model.c(ismember(model.rxns,'DM_pe_hs(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Phosphatidylethanolamine demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Phosphatidylglycerol demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'pglyc-hs(c)'});
    model.c(ismember(model.rxns,'DM_pglyc-hs(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Phosphatidylglycerol demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Phenylalanine demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'phe-L(c)'});
    model.c(ismember(model.rxns,'DM_phe-L(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Phenylalanine demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Proline demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'pro-L(c)'});
    model.c(ismember(model.rxns,'DM_pro-L(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Proline demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Phosphatidylserine demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'ps-hs(c)'});
    model.c(ismember(model.rxns,'DM_ps-hs(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Phosphatidylserine demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Serine demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'ser-L(c)'});
    model.c(ismember(model.rxns,'DM_ser-L(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Serine demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Sphingomyelin demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'sphmyln_hs(c)'});
    model.c(ismember(model.rxns,'DM_sphmyln_hs(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Sphingomyelin demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Threonine demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'thr-L(c)'});
    model.c(ismember(model.rxns,'DM_thr-L(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Threonine demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Tryptophan demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'trp-L(c)'});
    model.c(ismember(model.rxns,'DM_trp-L(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Tryptophan demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% Tyrosine demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'tyr-L(c)'});
    model.c(ismember(model.rxns,'DM_tyr-L(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Tyrosine demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA

    %% UTP demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'utp(c)'});
    model.c(ismember(model.rxns,'DM_utp(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'UTP demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
 
    %% Valine demand
    model = modelOri;
    model.c(find(model.c)) = 0;
    model.lb(ismember(model.rxns,'EX_o2(e)'))=-40;model.ub(ismember(model.rxns,'EX_o2(e)'))=-1;
    [model] = addDemandReaction(model, {'val-L(c)'});
    model.c(ismember(model.rxns,'DM_val-L(c)'))=1;
    
    if find(model.c)>0
        FBA = optimizeCbModel(model,'max','zero');
        TestSolution(k,1) = FBA.f;
    else
        TestSolution(k,1) = NaN;
    end
    TestSolutionName{k,1} = 'Valine demand';
    if ~isnan(TestSolution(k,1)); TestedRxns = [TestedRxns; model.rxns(find(abs(FBA.x)>tol))]; end ;k = k +1;clear FBA
 
end

TestSolutionName(:,2) = num2cell(TestSolution);
TestedRxns = unique(TestedRxns);
TestedRxns = intersect(modelOri.rxns,TestedRxns); % only those reactions that are also in modelOri not those that have been added to the network
PercTestedRxns = length(TestedRxns)*100/length(modelOri.rxns);

