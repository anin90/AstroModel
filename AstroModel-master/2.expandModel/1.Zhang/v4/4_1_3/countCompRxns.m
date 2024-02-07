tic;
model = MBA_model_Synapse;
%comp = {'[e]','[s]','[c]','[g]','[l]','[m]','[n]','[r]','[x]','[i]'};
e = length(findRxnFromCompartment(model, '[e]'));
s = length(findRxnFromCompartment(model, '[s]'));
c = length(findRxnFromCompartment(model, '[c]'));
g = length(findRxnFromCompartment(model, '[g]'));
l = length(findRxnFromCompartment(model, '[l]'));
m = length(findRxnFromCompartment(model, '[m]'));
n = length(findRxnFromCompartment(model, '[n]'));
r = length(findRxnFromCompartment(model, '[r]'));
x = length(findRxnFromCompartment(model, '[x]'));
i = length(findRxnFromCompartment(model, '[i]'));
toc;
clear model
clear e
clear s