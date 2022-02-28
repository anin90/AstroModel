
a = metabolicUnits;

for ix=2:length(a)
    if isempty(a{ix})
        a{ix} = a{ix-1};
    end
end

MetabolicUnits = a;
clear a ix

save('MetabolicUnits.mat');

% MU = table(MetabolicUnits, HarvettaRxns);
% [C,ia] = unique(MU(:,2:end),'rows');
% MetabolicUnits_ACS = MU(ia,:); clear MU C ia

% List1 = Recon3DModel.rxns;
% List2 = MetabolicUnits_ACS.HarvettaRxns;
% List3 = MetabolicUnits_ACS.MetabolicUnits;
% nameidx_L = getnameidx(List2, List1)';
% C = repmat({''},size(nameidx_L));
% C(nameidx_L~=0) = List3(nonzeros(nameidx_L));
% Recon3DModel.MetabolicUnits = C;
% Recon3DModel_MetabolicUnits = Recon3DModel;
% clear List1 List2 List3 nameidx_L C

% R = Recon3DModel_MetabolicUnits;
% unique(R.MetabolicUnits(~cellfun('isempty',R.MetabolicUnits)));


