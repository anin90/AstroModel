tStart = tic;

%% Load Data
% Load iAstro_Models
FileName   = 'iAstro_Models.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName

%% add path to dependencies:
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/1.Vadodaria/MTA_BD/')
    
%% DATASET: SOURCE to TAREGT
%% GSE32323: Healthy to Cancer (Recon1/Noam/FUT9) - GSE32323_a

%     load('/home/anirudh/Downloads/papers_to_read/mta/MSB-13-956-s008/EV1/model.mat');
%     load('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/Analysis_MTA_TestRuns/ExampleFiles_GSE32323/GSE32323.mat');
%     model = defineHumanMediaDMEM(model, flag);
%     model = prepareModel(model);
%     exprs.source = GSE32323.healthy_tab;
%     exprs.target = GSE32323.cancer_tab;
%     source = exprs.source; target = exprs.target;
%     smean = mean(exprs.source,2); tmean = mean(exprs.target,2);
%     [bin] = get_binary_expH(model, source, 1);
%     [~, Vref] = sampleiMAT(model, bin, 'gurobi', 'gurobi');
%     [discrete_rxns_vector, ~] = createDiscreteRxns_human(source, target, smean, tmean, model, 1);
%     [mta_tbl] = MTA_ImNotaGit(model, Vref, discrete_rxns_vector, 1:length(model.rxns), 'gurobi');
%     [mta_tbl_filt] = mta_tbl([mta_tbl.mta_score] > prctile([mta_tbl.mta_score]', 80));
%     [del_rxnID] = model.rxns([mta_tbl_filt.del_rxn]');
%     [mta_score] = [mta_tbl_filt.mta_score]';
%     [del_rxnID_FUT9] = intersect(findRxnsActiveWithGenes(model, '10690')', del_rxnID);

%% SOURCE_MODEL_ABS
%% Vadodaria_Ctrl: iPS Ctrl to Primary Ctrl

%     model = iAstro_iPS_Ctrl_TP_abs; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat');
%     data_Primary = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t1/GSE73721_HMA_CTX.mat');
%     dat_source = data_Ctrl; % SOURCE PHENOTYPE
%     dat_target = data_Primary; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_abs/Vadodaria_Ctrl_to_Primary.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Vadodaria_Ctrl: Primary Ctrl to iPS Ctrl

%     model = iAstro_Primary_TP_abs; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat');
%     data_Primary = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t1/GSE73721_HMA_CTX.mat');
%     dat_source = data_Primary; % SOURCE PHENOTYPE
%     dat_target = data_Ctrl; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_abs/Vadodaria_Primary_to_Ctrl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Vadodaria_BD: Ctrl to BD

%     model = iAstro_iPS_Ctrl_TP_abs; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat');
%     data_BD = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t1/Vadodaria_BD_Untreated.mat');
%     dat_source = data_Ctrl; % SOURCE PHENOTYPE
%     dat_target = data_BD; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_abs/Vadodaria_Ctrl_to_BD.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Vadodaria_BD: BD to Ctrl

%     model = iAstro_iPS_BD_TP_abs; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat');
%     data_BD = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t1/Vadodaria_BD_Untreated.mat');
%     dat_source = data_BD; % SOURCE PHENOTYPE
%     dat_target = data_Ctrl; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_abs/Vadodaria_BD_to_Ctrl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Vadodaria_BD_R: Ctrl to BD_R
 
%     model = iAstro_iPS_Ctrl_TP_abs; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat');
%     data_BD_R = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_Responder.mat');
%     dat_source = data_Ctrl; % SOURCE PHENOTYPE
%     dat_target = data_BD_R; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_abs/Vadodaria_Ctrl_to_BD_R.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Vadodaria_BD_R: BD_R to Ctrl

%     model = iAstro_iPS_BD_R_TP_abs; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat');
%     data_BD_R = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_Responder.mat');
%     dat_source = data_BD_R; % SOURCE PHENOTYPE
%     dat_target = data_Ctrl; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_abs/Vadodaria_BD_R_to_Ctrl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Vadodaria_BD_NR: Ctrl to BD_NR

%     model = iAstro_iPS_Ctrl_TP_abs; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat');
%     data_BD_NR = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_NonResponder.mat');
%     dat_source = data_Ctrl; % SOURCE PHENOTYPE
%     dat_target = data_BD_NR; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_abs/Vadodaria_Ctrl_to_BD_NR.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Vadodaria_BD_NR: BD_NR to Ctrl

%     model = iAstro_iPS_BD_NR_TP_abs; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat');
%     data_BD_NR = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_NonResponder.mat');
%     dat_source = data_BD_NR; % SOURCE PHENOTYPE
%     dat_target = data_Ctrl; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_abs/Vadodaria_BD_NR_to_Ctrl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Vadodaria: BD_R to BD_NR

%     model = iAstro_iPS_BD_R_TP_abs; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_BD_R = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_Responder.mat');
%     data_BD_NR = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_NonResponder.mat');
%     dat_source = data_BD_R; % SOURCE PHENOTYPE
%     dat_target = data_BD_NR; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_abs/Vadodaria_R_to_NR.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Vadodaria: BD_NR to BD_R

%     model = iAstro_iPS_BD_NR_TP_abs; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_BD_R = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_Responder.mat');
%     data_BD_NR = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_NonResponder.mat');
%     dat_source = data_BD_NR; % SOURCE PHENOTYPE
%     dat_target = data_BD_R; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_abs/Vadodaria_NR_to_R.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% SOURCE_MODEL_NORM_T1
%% Vadodaria_Ctrl: iPS Ctrl to Primary Ctrl

%     model = iAstro_iPS_Ctrl_TP_norm_t1; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat');
%     data_Primary = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t1/GSE73721_HMA_CTX.mat');
%     dat_source = data_Ctrl; % SOURCE PHENOTYPE
%     dat_target = data_Primary; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t1/Vadodaria_Ctrl_to_Primary.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Vadodaria_Ctrl: Primary Ctrl to iPS Ctrl

%     model = iAstro_Primary_TP_norm_t1; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat');
%     data_Primary = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t1/GSE73721_HMA_CTX.mat');
%     dat_source = data_Primary; % SOURCE PHENOTYPE
%     dat_target = data_Ctrl; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t1/Vadodaria_Primary_to_Ctrl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Vadodaria_BD: Ctrl to BD

%     model = iAstro_iPS_Ctrl_TP_norm_t1; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat');
%     data_BD = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t1/Vadodaria_BD_Untreated.mat');
%     dat_source = data_Ctrl; % SOURCE PHENOTYPE
%     dat_target = data_BD; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t1/Vadodaria_Ctrl_to_BD.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Vadodaria_BD: BD to Ctrl

%     model = iAstro_iPS_BD_TP_norm_t1; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat');
%     data_BD = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t1/Vadodaria_BD_Untreated.mat');
%     dat_source = data_BD; % SOURCE PHENOTYPE
%     dat_target = data_Ctrl; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t1/Vadodaria_BD_to_Ctrl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Vadodaria_BD_R: Ctrl to BD_R
 
%     model = iAstro_iPS_Ctrl_TP_norm_t1; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat');
%     data_BD_R = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_Responder.mat');
%     dat_source = data_Ctrl; % SOURCE PHENOTYPE
%     dat_target = data_BD_R; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t1/Vadodaria_Ctrl_to_BD_R.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Vadodaria_BD_R: BD_R to Ctrl

%     model = iAstro_iPS_BD_R_TP_norm_t1; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat');
%     data_BD_R = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_Responder.mat');
%     dat_source = data_BD_R; % SOURCE PHENOTYPE
%     dat_target = data_Ctrl; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t1/Vadodaria_BD_R_to_Ctrl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Vadodaria_BD_NR: Ctrl to BD_NR

%     model = iAstro_iPS_Ctrl_TP_norm_t1; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat');
%     data_BD_NR = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_NonResponder.mat');
%     dat_source = data_Ctrl; % SOURCE PHENOTYPE
%     dat_target = data_BD_NR; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t1/Vadodaria_Ctrl_to_BD_NR.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Vadodaria_BD_NR: BD_NR to Ctrl

%     model = iAstro_iPS_BD_NR_TP_norm_t1; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat');
%     data_BD_NR = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_NonResponder.mat');
%     dat_source = data_BD_NR; % SOURCE PHENOTYPE
%     dat_target = data_Ctrl; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t1/Vadodaria_BD_NR_to_Ctrl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Vadodaria: BD_R to BD_NR

%     model = iAstro_iPS_BD_R_TP_norm_t1; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_BD_R = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_Responder.mat');
%     data_BD_NR = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_NonResponder.mat');
%     dat_source = data_BD_R; % SOURCE PHENOTYPE
%     dat_target = data_BD_NR; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t1/Vadodaria_R_to_NR.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Vadodaria: BD_NR to BD_R

%     model = iAstro_iPS_BD_NR_TP_norm_t1; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_BD_R = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_Responder.mat');
%     data_BD_NR = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_NonResponder.mat');
%     dat_source = data_BD_NR; % SOURCE PHENOTYPE
%     dat_target = data_BD_R; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t1/Vadodaria_NR_to_R.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% SOURCE_MODEL_NORM_T2
%% Vadodaria_Ctrl: iPS Ctrl to Primary Ctrl

%     model = iAstro_iPS_Ctrl_TP_norm_t2; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat');
%     data_Primary = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t1/GSE73721_HMA_CTX.mat');
%     dat_source = data_Ctrl; % SOURCE PHENOTYPE
%     dat_target = data_Primary; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t2/Vadodaria_Ctrl_to_Primary.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Vadodaria_Ctrl: Primary Ctrl to iPS Ctrl

%     model = iAstro_Primary_TP_norm_t2; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat');
%     data_Primary = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t1/GSE73721_HMA_CTX.mat');
%     dat_source = data_Primary; % SOURCE PHENOTYPE
%     dat_target = data_Ctrl; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t2/Vadodaria_Primary_to_Ctrl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Vadodaria_BD: Ctrl to BD

%     model = iAstro_iPS_Ctrl_TP_norm_t2; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat');
%     data_BD = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t1/Vadodaria_BD_Untreated.mat');
%     dat_source = data_Ctrl; % SOURCE PHENOTYPE
%     dat_target = data_BD; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t2/Vadodaria_Ctrl_to_BD.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Vadodaria_BD: BD to Ctrl

%     model = iAstro_iPS_BD_TP_norm_t2; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat');
%     data_BD = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/2.BD_Untreated/v3/norm_t1/Vadodaria_BD_Untreated.mat');
%     dat_source = data_BD; % SOURCE PHENOTYPE
%     dat_target = data_Ctrl; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t2/Vadodaria_BD_to_Ctrl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Vadodaria_BD_R: Ctrl to BD_R
 
%     model = iAstro_iPS_Ctrl_TP_norm_t2; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat');
%     data_BD_R = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_Responder.mat');
%     dat_source = data_Ctrl; % SOURCE PHENOTYPE
%     dat_target = data_BD_R; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t2/Vadodaria_Ctrl_to_BD_R.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Vadodaria_BD_R: BD_R to Ctrl

%     model = iAstro_iPS_BD_R_TP_norm_t2; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat');
%     data_BD_R = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_Responder.mat');
%     dat_source = data_BD_R; % SOURCE PHENOTYPE
%     dat_target = data_Ctrl; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t2/Vadodaria_BD_R_to_Ctrl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Vadodaria_BD_NR: Ctrl to BD_NR

%     model = iAstro_iPS_Ctrl_TP_norm_t2; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat');
%     data_BD_NR = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_NonResponder.mat');
%     dat_source = data_Ctrl; % SOURCE PHENOTYPE
%     dat_target = data_BD_NR; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t2/Vadodaria_Ctrl_to_BD_NR.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Vadodaria_BD_NR: BD_NR to Ctrl

%     model = iAstro_iPS_BD_NR_TP_norm_t2; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/1.Control_Untreated/v3/norm_t1/Vadodaria_Control_Untreated.mat');
%     data_BD_NR = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_NonResponder.mat');
%     dat_source = data_BD_NR; % SOURCE PHENOTYPE
%     dat_target = data_Ctrl; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t2/Vadodaria_BD_NR_to_Ctrl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Vadodaria: BD_R to BD_NR

%     model = iAstro_iPS_BD_R_TP_norm_t2; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_BD_R = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_Responder.mat');
%     data_BD_NR = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_NonResponder.mat');
%     dat_source = data_BD_R; % SOURCE PHENOTYPE
%     dat_target = data_BD_NR; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t2/Vadodaria_R_to_NR.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Vadodaria: BD_NR to BD_R

%     model = iAstro_iPS_BD_NR_TP_norm_t2; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_BD_R = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/3.BD_Responder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_Responder.mat');
%     data_BD_NR = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/2.Vadodaria/4.BD_NonResponder_Untreated/v3/norm_t1/Vadodaria_BD_Untreated_NonResponder.mat');
%     dat_source = data_BD_NR; % SOURCE PHENOTYPE
%     dat_target = data_BD_R; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t2/Vadodaria_NR_to_R.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%%
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));

%%
clearvars -except mta_tbl_filt
