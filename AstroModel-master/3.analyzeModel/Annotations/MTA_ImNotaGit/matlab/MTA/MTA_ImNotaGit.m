%MTA:

%Input - 1. The Model iwht a given media
%        2. v_ref -  a reference vector describing the source state as derived from an iMAT-based analysis
%        3. discrete_rxns_vector - A discrete vector (1,-1,0), describing reactions that their
%           flux should increase/decrease (1,-1, respectively) or
%           remain constant (0) in order to transform to the
%           target state
%        4. rxns_to_delete - A set of reactions to KO, for all reactions, set it to 0:3788 for recon1 (model.mat); 0 means also run the control, where no reaction is deleted

%Output - 1. score - The score obtained by MTA following the KO of each of
%                    the reactions in rxns_to_delete
%         2. stat -  The solver's returned status for each KO
%         3. v_res - The optimal flux vector returned by the solver for each KO, thus a vector of the number of reactions by the number of KO's

function[y] = MTA_ImNotaGit(model, v_ref, discrete_rxns_vector, rxns_to_delete, solver)

DEFINE_PARAM;
alpha = 0.9;
thr = 0.01; %threshold defining the the minimum required change for the integer constraints
e_rxns = find(discrete_rxns_vector==1);
r_rxns = find(discrete_rxns_vector==-1);
s_rxns = find(discrete_rxns_vector==0 & model.c~=1);

fwd = [intersect(find(v_ref >= 0),e_rxns); intersect(find(v_ref < 0),r_rxns)];
bck = [intersect(find(v_ref <= 0),e_rxns); intersect(find(v_ref > 0),r_rxns)];

%Building the constraints matrix
[~,n_s] = size(model.S);
for i=1:length(fwd)
    model = addFWDCons(fwd(i),v_ref(fwd(i)),model,thr);
end


for i=1:length(bck)
    if v_ref(bck(i))==0 && model.lb(bck(i))==0
        continue;
    end
    model = addBCKCons(bck(i),v_ref(bck(i)),model,thr);
end

cons_rxns = [fwd; bck]';

%Constructing the optimization integer vector, QP vector and matrix
[~,n] = size(model.S);
model.int_vars(n_s+1:n) = 1;
model.c = zeros(n,1);
model.c(n_s+2:2:n) = alpha/2;%discrete_rxns_p(cons_rxns);
model.c(s_rxns) = -2*v_ref(s_rxns).*(1-alpha);
vec = zeros(n,1);
vec(s_rxns) = 2*(1-alpha);
model.F = diag(vec);

switch solver
    case 'tomlab'
        solv = @RunTomlabMIQP;
    case 'cplex'
        solv = @RunCplexMIQP;
    case 'gurobi'
        solv = @RunGurobiMIQP;
end

%v_res = zeros([n length(rxns_to_delete)]);

parfor i=1:length(rxns_to_delete)
    tmp = model;
    
    if rxns_to_delete(i)~=0
        % When rxns_to_delete(i)==0, it means the control where no reaction is deleted;
        % when it's not 0, perform the KO
        tmp.lb(rxns_to_delete(i)) = 0;
        tmp.ub(rxns_to_delete(i)) = 0;
    end
    
    %Running the optimization
    Res = solv(tmp, 1);
    
    %stat(i) = Res.result_status;
    %v_res(:,i) = Res.result_vector;

    if (isnan(Res.result_vector))
        %score(i) = NaN;
        x = struct();
        x.del_rxn = rxns_to_delete(i);
        x.solv_stat = NaN;
        x.obj_opt = NaN;
        x.v_opt = NaN;
        x.int_opt = NaN;
        x.rxns_change_yes = NaN;
        x.rxns_change_no = NaN;
        x.rxns_change_overdo = NaN;
        x.advs_change_yes = NaN;
        x.advs_change_no = NaN;
        x.advs_change_overdo = NaN;
        x.advs_steady = NaN;
        x.score_change = NaN;
        x.score_steady = NaN;
        x.mta_score = NaN;
        x.alt_score = NaN;
    else
        %Calculating the score
        %diff_change = calculateDiff(v_ref,fwd,bck,cons_rxns,Res.result_vector);
        %diff_steady = sum(abs(Res.result_vector(s_rxns)-v_ref(s_rxns)));
        %score(i) = diff_change/diff_steady;
        [success, unsuccess, wrong, dvs, dvu, dvw, diffc] = calculateDiff(v_ref, fwd, bck, cons_rxns, Res.result_vector);
        x = struct();
        x.del_rxn = rxns_to_delete(i);
        x.solv_stat = Res.result_status;
        x.obj_opt = Res.result_opt;
        x.v_opt = Res.result_vector(1:length(v_ref));
        x.int_opt = Res.result_vector(length(v_ref)+1:length(Res.result_vector));
        x.rxns_change_yes = success;
        x.rxns_change_no = unsuccess;
        x.rxns_change_overdo = wrong;
        x.advs_change_yes = dvs;
        x.advs_change_no = dvu;
        x.advs_change_overdo = dvw;
        x.advs_steady = abs(Res.result_vector(s_rxns)-v_ref(s_rxns));
        x.score_change = diffc;
        x.score_steady = sum(x.advs_steady);
        x.mta_score = diffc/x.score_steady;
        x.alt_score = diffc/(length(success)+length(unsuccess)+length(wrong)) - x.score_steady/length(s_rxns);
    end

    y(i) = x;
end
