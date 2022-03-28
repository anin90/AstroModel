mta_path = '/home/fountain/Documents/Projects/ourtools/MTA';
cplex_path = '/home/fountain/ibm/ILOG/CPLEX_Studio128/cplex/matlab/x86-64_linux';
tomlab_path = '/home/fountain/tomlab';
gurobi_path = '/home/fountain/gurobi800/linux64/matlab';
cobra_path = '/home/fountain/biocomp_tools/cobratoolbox';
cur_dir = pwd;

% add cplex path
addpath(cplex_path)
% initiate tomlab
cd(tomlab_path)
startup
% initiate gurobi
cd(gurobi_path)
gurobi_setup
% initiate cobratoolbox
cd(cobra_path)
initCobraToolbox
% add MTA paths
addpath(fullfile(mta_path, 'matlab', 'iMAT'))
addpath(fullfile(mta_path, 'matlab', 'MTA'))
addpath(fullfile(mta_path, 'models'))
% return to current dir
cd(cur_dir)

clear mta_path cplex_path tomlab_path gurobi_path cobra_path cur_dir;
