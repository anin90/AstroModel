function[mets, uniqueMetNames, uniqueMetIDs] = printUniqueMets(model, rxnList);
%%
tic;
%%
mets = findMetsFromRxns(model, rxnList);

[intis,k] = ismember(model.mets,(mets));
[~,p] = sort(k(intis)); idx = find(intis); idx = idx(p);
uniqueMetNames = unique(model.metNames(idx));

[intis,k] = ismember(model.metNames,(uniqueMetNames));
[~,p] = sort(k(intis)); idx = find(intis); idx = idx(p);
metIDmaps = model.mets(idx);

newA = cellfun(@(x) strsplit(x, '['), metIDmaps, 'UniformOutput', false);
newA = vertcat(newA{:}); 
uniqueMetIDs = unique(newA(:,1), 'stable');

%%
toc;
end