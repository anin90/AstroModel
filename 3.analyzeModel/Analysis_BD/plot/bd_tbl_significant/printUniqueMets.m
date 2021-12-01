function[mets_tbl_unique] = printUniqueMets(model, rxnList)
%%
tic;
%%
mets = findMetsFromRxns(model, rxnList);
[intis,k] = ismember(model.mets,(mets));
[~,p] = sort(k(intis)); idx = find(intis); idx = idx(p);
metNames = model.metNames(idx);
mets = [mets metNames];
newA = cellfun(@(x) strsplit(x, '['), mets(:,1), 'UniformOutput', false);
newA = vertcat(newA{:});
mets_newA = [newA(:,1) mets(:,2)];
mets_tbl = cell2table(mets_newA,'VariableNames', {'metID','metName'});
[~,idu] = unique(mets_tbl(:,1));
mets_tbl_unique = mets_tbl(idu,:);

%%
toc;
end