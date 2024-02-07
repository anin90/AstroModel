function [sample_points, warmup_points, state]= model_sample(model, samples_count, lp_solver, warmup_points, state)
    if nargin < 4
        warmup_points = model_sampling_warmup(model, lp_solver, 5000, 1);
    end;
    if nargin < 5
        [~, state] = model_sampling_ACHR(model, warmup_points,  1000, 400);
    end
    [sample_points, state] = model_sampling_ACHR(model, warmup_points,  samples_count, 400, state);
end
