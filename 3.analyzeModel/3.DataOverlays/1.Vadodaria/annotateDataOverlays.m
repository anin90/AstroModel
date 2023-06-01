%%
tStart = tic;

%% add path to dependencies:
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/1.Vadodaria/PlotResults/')
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/1.Vadodaria/PlotResults/plotDataOverlaysHyper_Tbl/')
addpath('/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/3.DataOverlays/1.Vadodaria/PlotResults/plotDataOverlaysHyper_Tbl_Final/')

%% BD_Lumped

    % BD_S1
    Tbl = readtable('BD_S1.csv', 'Delimiter', '\t');
    subSystem = Tbl.subSystem;
    disrupted_module = Tbl.disrupted_module;
    rxnList = Tbl.rxnList;
    [rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers] = annotateRxnList(rxnList);
    Tbl = [subSystem, disrupted_module, rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers];
    Tbl = cell2table(Tbl);
    Tbl.Properties.VariableNames = {'subSystem', 'disrupted_module', 'rxnList', 'MetabolicUnits', 'Localization', 'rxnFormula', 'rxnName', 'rxnECNumbers'};
    Tbl(:,'Localization') = [];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_S1.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % BD_S2
    Tbl = readtable('BD_S2.csv', 'Delimiter', '\t');
    subSystem = Tbl.subSystem;
    disrupted_module = Tbl.disrupted_module;
    rxnList = Tbl.rxnList;
    [rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers] = annotateRxnList(rxnList);
    Tbl = [subSystem, disrupted_module, rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers];
    Tbl = cell2table(Tbl);
    Tbl.Properties.VariableNames = {'subSystem', 'disrupted_module', 'rxnList', 'MetabolicUnits', 'Localization', 'rxnFormula', 'rxnName', 'rxnECNumbers'};
    Tbl(:,'Localization') = [];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_S2.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % BD_S3
    Tbl = readtable('BD_S3.csv', 'Delimiter', '\t');
    subSystem = Tbl.subSystem;
    disrupted_module = Tbl.disrupted_module;
    rxnList = Tbl.rxnList;
    [rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers] = annotateRxnList(rxnList);
    Tbl = [subSystem, disrupted_module, rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers];
    Tbl = cell2table(Tbl);
    Tbl.Properties.VariableNames = {'subSystem', 'disrupted_module', 'rxnList', 'MetabolicUnits', 'Localization', 'rxnFormula', 'rxnName', 'rxnECNumbers'};
    Tbl(:,'Localization') = [];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_S3.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% BD_R

    % BD_R_S1
    Tbl = readtable('BD_R_S1.csv', 'Delimiter', '\t');
    subSystem = Tbl.subSystem;
    disrupted_module = Tbl.disrupted_module;
    rxnList = Tbl.rxnList;
    [rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers] = annotateRxnList(rxnList);
    Tbl = [subSystem, disrupted_module, rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers];
    Tbl = cell2table(Tbl);
    Tbl.Properties.VariableNames = {'subSystem', 'disrupted_module', 'rxnList', 'MetabolicUnits', 'Localization', 'rxnFormula', 'rxnName', 'rxnECNumbers'};
    Tbl(:,'Localization') = [];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_R_S1.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % BD_R_S2
    Tbl = readtable('BD_R_S2.csv', 'Delimiter', '\t');
    subSystem = Tbl.subSystem;
    disrupted_module = Tbl.disrupted_module;
    rxnList = Tbl.rxnList;
    [rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers] = annotateRxnList(rxnList);
    Tbl = [subSystem, disrupted_module, rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers];
    Tbl = cell2table(Tbl);
    Tbl.Properties.VariableNames = {'subSystem', 'disrupted_module', 'rxnList', 'MetabolicUnits', 'Localization', 'rxnFormula', 'rxnName', 'rxnECNumbers'};
    Tbl(:,'Localization') = [];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_R_S2.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % BD_R_S3
    Tbl = readtable('BD_R_S3.csv', 'Delimiter', '\t');
    subSystem = Tbl.subSystem;
    disrupted_module = Tbl.disrupted_module;
    rxnList = Tbl.rxnList;
    [rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers] = annotateRxnList(rxnList);
    Tbl = [subSystem, disrupted_module, rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers];
    Tbl = cell2table(Tbl);
    Tbl.Properties.VariableNames = {'subSystem', 'disrupted_module', 'rxnList', 'MetabolicUnits', 'Localization', 'rxnFormula', 'rxnName', 'rxnECNumbers'};
    Tbl(:,'Localization') = [];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_R_S3.csv', 'WriteVariableNames', true, 'Delimiter','\t');

%% BD_NR

    % BD_NR_S1
    Tbl = readtable('BD_NR_S1.csv', 'Delimiter', '\t');
    subSystem = Tbl.subSystem;
    disrupted_module = Tbl.disrupted_module;
    rxnList = Tbl.rxnList;
    [rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers] = annotateRxnList(rxnList);
    Tbl = [subSystem, disrupted_module, rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers];
    Tbl = cell2table(Tbl);
    Tbl.Properties.VariableNames = {'subSystem', 'disrupted_module', 'rxnList', 'MetabolicUnits', 'Localization', 'rxnFormula', 'rxnName', 'rxnECNumbers'};
    Tbl(:,'Localization') = [];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_NR_S1.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % BD_NR_S2
    Tbl = readtable('BD_NR_S2.csv', 'Delimiter', '\t');
    subSystem = Tbl.subSystem;
    disrupted_module = Tbl.disrupted_module;
    rxnList = Tbl.rxnList;
    [rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers] = annotateRxnList(rxnList);
    Tbl = [subSystem, disrupted_module, rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers];
    Tbl = cell2table(Tbl);
    Tbl.Properties.VariableNames = {'subSystem', 'disrupted_module', 'rxnList', 'MetabolicUnits', 'Localization', 'rxnFormula', 'rxnName', 'rxnECNumbers'};
    Tbl(:,'Localization') = [];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_NR_S2.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % BD_NR_S3
    Tbl = readtable('BD_NR_S3.csv', 'Delimiter', '\t');
    subSystem = Tbl.subSystem;
    disrupted_module = Tbl.disrupted_module;
    rxnList = Tbl.rxnList;
    [rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers] = annotateRxnList(rxnList);
    Tbl = [subSystem, disrupted_module, rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers];
    Tbl = cell2table(Tbl);
    Tbl.Properties.VariableNames = {'subSystem', 'disrupted_module', 'rxnList', 'MetabolicUnits', 'Localization', 'rxnFormula', 'rxnName', 'rxnECNumbers'};
    Tbl(:,'Localization') = [];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_NR_S3.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % BD_NR_S4
    Tbl = readtable('BD_NR_S4.csv', 'Delimiter', '\t');
    subSystem = Tbl.subSystem;
    disrupted_module = Tbl.disrupted_module;
    rxnList = Tbl.rxnList;
    [rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers] = annotateRxnList(rxnList);
    Tbl = [subSystem, disrupted_module, rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers];
    Tbl = cell2table(Tbl);
    Tbl.Properties.VariableNames = {'subSystem', 'disrupted_module', 'rxnList', 'MetabolicUnits', 'Localization', 'rxnFormula', 'rxnName', 'rxnECNumbers'};
    Tbl(:,'Localization') = [];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_NR_S4.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % BD_NR_S5
    Tbl = readtable('BD_NR_S5.csv', 'Delimiter', '\t');
    subSystem = Tbl.subSystem;
    disrupted_module = Tbl.disrupted_module;
    rxnList = Tbl.rxnList;
    [rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers] = annotateRxnList(rxnList);
    Tbl = [subSystem, disrupted_module, rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers];
    Tbl = cell2table(Tbl);
    Tbl.Properties.VariableNames = {'subSystem', 'disrupted_module', 'rxnList', 'MetabolicUnits', 'Localization', 'rxnFormula', 'rxnName', 'rxnECNumbers'};
    Tbl(:,'Localization') = [];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_NR_S5.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % BD_NR_S6
    Tbl = readtable('BD_NR_S6.csv', 'Delimiter', '\t');
    subSystem = Tbl.subSystem;
    disrupted_module = Tbl.disrupted_module;
    rxnList = Tbl.rxnList;
    [rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers] = annotateRxnList(rxnList);
    Tbl = [subSystem, disrupted_module, rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers];
    Tbl = cell2table(Tbl);
    Tbl.Properties.VariableNames = {'subSystem', 'disrupted_module', 'rxnList', 'MetabolicUnits', 'Localization', 'rxnFormula', 'rxnName', 'rxnECNumbers'};
    Tbl(:,'Localization') = [];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_NR_S6.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % BD_NR_S7
    Tbl = readtable('BD_NR_S7.csv', 'Delimiter', '\t');
    subSystem = Tbl.subSystem;
    disrupted_module = Tbl.disrupted_module;
    rxnList = Tbl.rxnList;
    [rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers] = annotateRxnList(rxnList);
    Tbl = [subSystem, disrupted_module, rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers];
    Tbl = cell2table(Tbl);
    Tbl.Properties.VariableNames = {'subSystem', 'disrupted_module', 'rxnList', 'MetabolicUnits', 'Localization', 'rxnFormula', 'rxnName', 'rxnECNumbers'};
    Tbl(:,'Localization') = [];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_NR_S7.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % BD_NR_S8
    Tbl = readtable('BD_NR_S8.csv', 'Delimiter', '\t');
    subSystem = Tbl.subSystem;
    disrupted_module = Tbl.disrupted_module;
    rxnList = Tbl.rxnList;
    [rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers] = annotateRxnList(rxnList);
    Tbl = [subSystem, disrupted_module, rxnList, MetabolicUnits, Localization, RxnFormula, RxnName, RxnECNumbers];
    Tbl = cell2table(Tbl);
    Tbl.Properties.VariableNames = {'subSystem', 'disrupted_module', 'rxnList', 'MetabolicUnits', 'Localization', 'rxnFormula', 'rxnName', 'rxnECNumbers'};
    Tbl(:,'Localization') = [];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/BD_NR_S8.csv', 'WriteVariableNames', true, 'Delimiter','\t');
    
%% DF_X_BD_Minerva

    % DF_1_BD_S1_Minerva
    Tbl = readtable('DF_1_BD_S1_Minerva.csv', 'Delimiter', '\t');
    Tbl.Var1=[];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/DF_1_BD_S1_Minerva.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % DF_2_BD_S2_Minerva
    Tbl = readtable('DF_2_BD_S2_Minerva.csv', 'Delimiter', '\t');
    Tbl.Var1=[];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/DF_2_BD_S2_Minerva.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % DF_3_BD_S3_Minerva
    Tbl = readtable('DF_3_BD_S3_Minerva.csv', 'Delimiter', '\t');
    Tbl.Var1=[];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/DF_3_BD_S3_Minerva.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % DF_4_BD_R_S1_Minerva
    Tbl = readtable('DF_4_BD_R_S1_Minerva.csv', 'Delimiter', '\t');
    Tbl.Var1=[];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/DF_4_BD_R_S1_Minerva.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % DF_5_BD_R_S2_Minerva
    Tbl = readtable('DF_5_BD_R_S2_Minerva.csv', 'Delimiter', '\t');
    Tbl.Var1=[];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/DF_5_BD_R_S2_Minerva.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % DF_6_BD_R_S3_Minerva
    Tbl = readtable('DF_6_BD_R_S3_Minerva.csv', 'Delimiter', '\t');
    Tbl.Var1=[];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/DF_6_BD_R_S3_Minerva.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % DF_7_BD_NR_S1_Minerva
    Tbl = readtable('DF_7_BD_NR_S1_Minerva.csv', 'Delimiter', '\t');
    Tbl.Var1=[];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/DF_7_BD_NR_S1_Minerva.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % DF_8_BD_NR_S2_Minerva
    Tbl = readtable('DF_8_BD_NR_S2_Minerva.csv', 'Delimiter', '\t');
    Tbl.Var1=[];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/DF_8_BD_NR_S2_Minerva.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % DF_9_BD_NR_S3_Minerva
    Tbl = readtable('DF_9_BD_NR_S3_Minerva.csv', 'Delimiter', '\t');
    Tbl.Var1=[];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/DF_9_BD_NR_S3_Minerva.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % DF_10_BD_NR_S4_Minerva
    Tbl = readtable('DF_10_BD_NR_S4_Minerva.csv', 'Delimiter', '\t');
    Tbl.Var1=[];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/DF_10_BD_NR_S4_Minerva.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % DF_11_BD_NR_S5_Minerva
    Tbl = readtable('DF_11_BD_NR_S5_Minerva.csv', 'Delimiter', '\t');
    Tbl.Var1=[];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/DF_11_BD_NR_S5_Minerva.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % DF_12_BD_NR_S6_Minerva
    Tbl = readtable('DF_12_BD_NR_S6_Minerva.csv', 'Delimiter', '\t');
    Tbl.Var1=[];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/DF_12_BD_NR_S6_Minerva.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % DF_13_BD_NR_S7_Minerva
    Tbl = readtable('DF_13_BD_NR_S7_Minerva.csv', 'Delimiter', '\t');
    Tbl.Var1=[];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/DF_13_BD_NR_S7_Minerva.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % DF_14_BD_NR_S8_Minerva
    Tbl = readtable('DF_14_BD_NR_S8_Minerva.csv', 'Delimiter', '\t');
    Tbl.Var1=[];
    writetable(Tbl, 'PlotResults/plotDataOverlaysHyper_Tbl_Final/DF_14_BD_NR_S8_Minerva.csv', 'WriteVariableNames', true, 'Delimiter','\t');
    
%% Final_BD_Minerva

    % BD_Rxns
    Tbl = readtable('BD_Rxns.csv', 'Delimiter', '\t');
    Tbl.Var1=[];
    writetable(Tbl, 'PlotResults/BD_Rxns.csv', 'WriteVariableNames', true, 'Delimiter','\t');
  
    % BD_R_Rxns
    Tbl = readtable('BD_R_Rxns.csv', 'Delimiter', '\t');
    Tbl.Var1=[];
    writetable(Tbl, 'PlotResults/BD_R_Rxns.csv', 'WriteVariableNames', true, 'Delimiter','\t');

    % BD_NR_Rxns
    Tbl = readtable('BD_NR_Rxns.csv', 'Delimiter', '\t');
    Tbl.Var1=[];
    writetable(Tbl, 'PlotResults/BD_NR_Rxns.csv', 'WriteVariableNames', true, 'Delimiter','\t');


%%
tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));    

%%
clearvars
