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
Test4HumanAstro_TP = table(TestResults_iAstro_Primary_TP.TestSolutionName, TestResults_iAstro_Primary_TP.TestSolution, TestResults_iAstro_iPS_Ctrl_TP.TestSolution, TestResults_iAstro_iPS_BD_TP.TestSolution, TestResults_iAstro_iPS_BD_R_TP.TestSolution, TestResults_iAstro_iPS_BD_NR_TP.TestSolution);
% writetable(Test4HumanAstro_TP,'Test4HumanAstro_TP.csv','WriteVariableNames', true, 'Delimiter','\t');

%%
clearvars -except Test4HumanAstro_TP
save('Test4HumanAstro_TP.mat');

%%
toc;

