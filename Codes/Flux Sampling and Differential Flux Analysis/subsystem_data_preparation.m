uirxns=modelUI.rxns;
irxns=modelI.rxns;
unrxns=intersect(uirxns,irxns);

subnames={};
k=0;
for i=1:length(unrxns)

    subs=ihuman.subSystems(ismember(ihuman.rxns,unrxns(i)));
    subnames=[subnames;vertcat(subs{:})];
    k=k+1;
end

data=table(unrxns,subnames);

writetable(data,'ReactionSetNHBE.csv')