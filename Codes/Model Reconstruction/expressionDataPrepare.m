function [val_UI,val_I]=expressionDataPrepare(key)

if(key==1)
    range_UI=3:5;
    range_I=6:8;
elseif(key==2)
    range_UI=50:53;
    range_I=54:57;
elseif(key==3)
    range_UI=68:69;
    range_I=70:71;
end

expr=readtable('MetTPMAll.csv');

geneExpr_UI=str2double(table2array(expr(:,range_UI))); %Gene expression in uninfected/mock cell line
geneExpr_I=str2double(table2array(expr(:,range_I))) ;%Gene expression in infected cell line

%Average gene expression accross replicates

val_UI=mean(geneExpr_UI,2);
val_I=mean(geneExpr_I,2);




end
