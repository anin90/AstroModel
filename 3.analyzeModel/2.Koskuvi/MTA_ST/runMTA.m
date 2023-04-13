tStart = tic;

%% Load Data
% Load iAstro_Models
FileName   = 'iAstro_Models.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName

%% add path to dependencies:
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/2.Koskuvi/MTA_ST/')
    
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
%% Koskuvi_Ctrl: iPS Ctrl to Primary Ctrl

%     model = iAstro_iPS_SCZ_Ctrl_TP_abs; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/1.Control/v1/norm_t1/Koskuvi_Control.mat');
%     data_Primary = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t1/GSE73721_HMA_CTX.mat');
%     dat_source = data_Ctrl; % SOURCE PHENOTYPE
%     dat_target = data_Primary; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_abs/Koskuvi_Ctrl_to_Primary.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Koskuvi_Ctrl: Primary Ctrl to iPS Ctrl

%     model = iAstro_Primary_TP_abs; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/1.Control/v1/norm_t1/Koskuvi_Control.mat');
%     data_Primary = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t1/GSE73721_HMA_CTX.mat');
%     dat_source = data_Primary; % SOURCE PHENOTYPE
%     dat_target = data_Ctrl; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_abs/Koskuvi_Primary_to_Ctrl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Koskuvi_HT: iPS Ctrl to iPS HT

%     model = iAstro_iPS_SCZ_Ctrl_TP_abs; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/1.Control/v1/norm_t1/Koskuvi_Control.mat');
%     data_HT = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/2.HT/v1/norm_t1/Koskuvi_HT.mat');
%     dat_source = data_Ctrl; % SOURCE PHENOTYPE
%     dat_target = data_HT; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_abs/Koskuvi_Ctrl_to_HT.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Koskuvi_HT: iPS HT to iPS Ctrl

%     model = iAstro_iPS_SCZ_HT_TP_abs; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/1.Control/v1/norm_t1/Koskuvi_Control.mat');
%     data_HT = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/2.HT/v1/norm_t1/Koskuvi_HT.mat');
%     dat_source = data_HT; % SOURCE PHENOTYPE
%     dat_target = data_Ctrl; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_abs/Koskuvi_HT_to_Ctrl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Koskuvi_ST: iPS Ctrl to iPS ST

%     model = iAstro_iPS_SCZ_Ctrl_TP_abs; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/1.Control/v1/norm_t1/Koskuvi_Control.mat');
%     data_ST = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/3.ST/v1/norm_t1/Koskuvi_ST.mat');
%     dat_source = data_Ctrl; % SOURCE PHENOTYPE
%     dat_target = data_ST; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_abs/Koskuvi_Ctrl_to_ST.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Koskuvi_HT: iPS ST to iPS Ctrl

%     model = iAstro_iPS_SCZ_ST_TP_abs; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/1.Control/v1/norm_t1/Koskuvi_Control.mat');
%     data_ST = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/3.ST/v1/norm_t1/Koskuvi_ST.mat');
%     dat_source = data_ST; % SOURCE PHENOTYPE
%     dat_target = data_Ctrl; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_abs/Koskuvi_ST_to_Ctrl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% SOURCE_MODEL_NORM_T1
%% Koskuvi_Ctrl: iPS Ctrl to Primary Ctrl

%     model = iAstro_iPS_SCZ_Ctrl_TP_norm_t1; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/1.Control/v1/norm_t1/Koskuvi_Control.mat');
%     data_Primary = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t1/GSE73721_HMA_CTX.mat');
%     dat_source = data_Ctrl; % SOURCE PHENOTYPE
%     dat_target = data_Primary; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t1/Koskuvi_Ctrl_to_Primary.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Koskuvi_Ctrl: Primary Ctrl to iPS Ctrl

%     model = iAstro_Primary_TP_norm_t1; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/1.Control/v1/norm_t1/Koskuvi_Control.mat');
%     data_Primary = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t1/GSE73721_HMA_CTX.mat');
%     dat_source = data_Primary; % SOURCE PHENOTYPE
%     dat_target = data_Ctrl; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t1/Koskuvi_Primary_to_Ctrl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Koskuvi_HT: iPS Ctrl to iPS HT

%     model = iAstro_iPS_SCZ_Ctrl_TP_norm_t1; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/1.Control/v1/norm_t1/Koskuvi_Control.mat');
%     data_HT = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/2.HT/v1/norm_t1/Koskuvi_HT.mat');
%     dat_source = data_Ctrl; % SOURCE PHENOTYPE
%     dat_target = data_HT; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t1/Koskuvi_Ctrl_to_HT.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Koskuvi_HT: iPS HT to iPS Ctrl

%     model = iAstro_iPS_SCZ_HT_TP_norm_t1; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/1.Control/v1/norm_t1/Koskuvi_Control.mat');
%     data_HT = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/2.HT/v1/norm_t1/Koskuvi_HT.mat');
%     dat_source = data_HT; % SOURCE PHENOTYPE
%     dat_target = data_Ctrl; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t1/Koskuvi_HT_to_Ctrl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Koskuvi_ST: iPS Ctrl to iPS ST

%     model = iAstro_iPS_SCZ_Ctrl_TP_norm_t1; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/1.Control/v1/norm_t1/Koskuvi_Control.mat');
%     data_ST = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/3.ST/v1/norm_t1/Koskuvi_ST.mat');
%     dat_source = data_Ctrl; % SOURCE PHENOTYPE
%     dat_target = data_ST; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t1/Koskuvi_Ctrl_to_ST.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Koskuvi_HT: iPS ST to iPS Ctrl

%     model = iAstro_iPS_SCZ_ST_TP_norm_t1; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/1.Control/v1/norm_t1/Koskuvi_Control.mat');
%     data_ST = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/3.ST/v1/norm_t1/Koskuvi_ST.mat');
%     dat_source = data_ST; % SOURCE PHENOTYPE
%     dat_target = data_Ctrl; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t1/Koskuvi_ST_to_Ctrl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% SOURCE_MODEL_NORM_T2
%% Koskuvi_Ctrl: iPS Ctrl to Primary Ctrl

%     model = iAstro_iPS_SCZ_Ctrl_TP_norm_t2; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/1.Control/v1/norm_t1/Koskuvi_Control.mat');
%     data_Primary = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t1/GSE73721_HMA_CTX.mat');
%     dat_source = data_Ctrl; % SOURCE PHENOTYPE
%     dat_target = data_Primary; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t2/Koskuvi_Ctrl_to_Primary.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Koskuvi_Ctrl: Primary Ctrl to iPS Ctrl

%     model = iAstro_Primary_TP_norm_t2; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/1.Control/v1/norm_t1/Koskuvi_Control.mat');
%     data_Primary = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/1.Zhang/v12/norm_t1/GSE73721_HMA_CTX.mat');
%     dat_source = data_Primary; % SOURCE PHENOTYPE
%     dat_target = data_Ctrl; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t2/Koskuvi_Primary_to_Ctrl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Koskuvi_HT: iPS Ctrl to iPS HT

%     model = iAstro_iPS_SCZ_Ctrl_TP_norm_t2; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/1.Control/v1/norm_t1/Koskuvi_Control.mat');
%     data_HT = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/2.HT/v1/norm_t1/Koskuvi_HT.mat');
%     dat_source = data_Ctrl; % SOURCE PHENOTYPE
%     dat_target = data_HT; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t2/Koskuvi_Ctrl_to_HT.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Koskuvi_HT: iPS HT to iPS Ctrl

%     model = iAstro_iPS_SCZ_HT_TP_norm_t2; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/1.Control/v1/norm_t1/Koskuvi_Control.mat');
%     data_HT = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/2.HT/v1/norm_t1/Koskuvi_HT.mat');
%     dat_source = data_HT; % SOURCE PHENOTYPE
%     dat_target = data_Ctrl; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t2/Koskuvi_HT_to_Ctrl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Koskuvi_ST: iPS Ctrl to iPS ST

%     model = iAstro_iPS_SCZ_Ctrl_TP_norm_t2; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/1.Control/v1/norm_t1/Koskuvi_Control.mat');
%     data_ST = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/3.ST/v1/norm_t1/Koskuvi_ST.mat');
%     dat_source = data_Ctrl; % SOURCE PHENOTYPE
%     dat_target = data_ST; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t2/Koskuvi_Ctrl_to_ST.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% Koskuvi_HT: iPS ST to iPS Ctrl

%     model = iAstro_iPS_SCZ_ST_TP_norm_t2; % SOURCE MODEL
%     model = changeObjective(model, 'biomass_maintenance', 1);
%     data_Ctrl = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/1.Control/v1/norm_t1/Koskuvi_Control.mat');
%     data_ST = importdata('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/1.matrix2model/3.Koskuvi/3.ST/v1/norm_t1/Koskuvi_ST.mat');
%     dat_source = data_ST; % SOURCE PHENOTYPE
%     dat_target = data_Ctrl; % TARGET PHENOTYPE
%     [mta_tbl] = generateMTAscores(model, dat_source, dat_target);
%     writetable(mta_tbl, 'PlotResults/mta_tbl_norm_t2/Koskuvi_ST_to_Ctrl.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%%
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));

%%
clearvars -except mta_tbl_filt
