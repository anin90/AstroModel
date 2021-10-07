%% Load Data
tic;
% Load iPS_Models 
FileName   = 'iAstro_iPS_Models.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/';
File = fullfile(FolderName, FileName);
load(File);
% Load Primary_Models
FileName   = 'iAstro_Primary_Models.mat';
FolderName = '/media/anirudh/Work/ADBS_NIMHANS/Thesis/1.Science/Analysis/cobratoolbox/AstroModel/3.analyzeModel/';
File = fullfile(FolderName, FileName);
load(File);
clear File FileName FolderName

%% model_TP
[TestResults_iAstro_Primary_TP] = Test4HumanAstro(iAstro_Primary_TP);
[TestResults_iAstro_iPS_Ctrl_TP] = Test4HumanAstro(iAstro_iPS_Ctrl_TP);
[TestResults_iAstro_iPS_BD_TP] = Test4HumanAstro(iAstro_iPS_BD_TP);
[TestResults_iAstro_iPS_BD_R_TP] = Test4HumanAstro(iAstro_iPS_BD_R_TP);
[TestResults_iAstro_iPS_BD_NR_TP] = Test4HumanAstro(iAstro_iPS_BD_NR_TP);
VarNames_TP = {'TestSolutionName','Primary_TP','iPS_Ctrl_TP','iPS_BD_TP','iPS_BD_R_TP','iPS_BD_NR_TP'};
Test4HumanAstro_TP = table(TestResults_iAstro_Primary_TP.TestSolutionName, TestResults_iAstro_Primary_TP.TestSolution, TestResults_iAstro_iPS_Ctrl_TP.TestSolution, TestResults_iAstro_iPS_BD_TP.TestSolution, TestResults_iAstro_iPS_BD_R_TP.TestSolution, TestResults_iAstro_iPS_BD_NR_TP.TestSolution,'VariableNames', VarNames_TP);
writetable(Test4HumanAstro_TP,'Test4HumanAstro_TP.csv','WriteVariableNames', true, 'Delimiter','\t');

%% model_T
[TestResults_iAstro_Primary_T] = Test4HumanAstro(iAstro_Primary_T);
[TestResults_iAstro_iPS_Ctrl_T] = Test4HumanAstro(iAstro_iPS_Ctrl_T);
[TestResults_iAstro_iPS_BD_T] = Test4HumanAstro(iAstro_iPS_BD_T);
[TestResults_iAstro_iPS_BD_R_T] = Test4HumanAstro(iAstro_iPS_BD_R_T);
[TestResults_iAstro_iPS_BD_NR_T] = Test4HumanAstro(iAstro_iPS_BD_NR_T);
VarNames_T = {'TestSolutionName','Primary_T','iPS_Ctrl_T','iPS_BD_T','iPS_BD_R_T','iPS_BD_NR_T'};
Test4HumanAstro_T = table(TestResults_iAstro_Primary_T.TestSolutionName, TestResults_iAstro_Primary_T.TestSolution, TestResults_iAstro_iPS_Ctrl_T.TestSolution, TestResults_iAstro_iPS_BD_T.TestSolution, TestResults_iAstro_iPS_BD_R_T.TestSolution, TestResults_iAstro_iPS_BD_NR_T.TestSolution,'VariableNames', VarNames_T);
writetable(Test4HumanAstro_T,'Test4HumanAstro_T.csv','WriteVariableNames', true, 'Delimiter','\t');

%%
clearvars -except Test4HumanAstro_TP Test4HumanAstro_T
save('Test4HumanAstro.mat');

%%
toc;

