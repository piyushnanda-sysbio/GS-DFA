function ihsarscov2=ReactionAdd(model)

%Adding reactions to closed model Human_GEM.mat

%     {'Extracellular'        } 1
%     {'Peroxisome'           } 2
%     {'Mitochondria'         } 3
%     {'Cytosol'              } 4
%     {'Lysosome'             } 5
%     {'Endoplasmic reticulum'} 6
%     {'Golgi apparatus'      } 7
%     {'Nucleus'              } 8
%     {'Boundary'             } 9
%     {'Inner mitochondria'   } 10

ih=model; %Crete a duplicate model

met='sarscov2c';

ihsarscov2=addMetabolite(ih,met);
ihsarscov2.metComps(ismember(ihsarscov2.mets,met))=4; %Add to cytosol
ihsarscov2.metNames(ismember(ihsarscov2.mets,met))={'VBOF'};
%Adding biomass equation

vbof=readtable('VBOF.csv');
metList=table2array(vbof(:,1));
stoicList=table2array(vbof(:,2));

ihsarscov2=addReaction(ihsarscov2,'Cov2VBOF','metaboliteList',metList,'stoichCoeffList',stoicList,'reversible',false);

%Add extracellular dummy metabolites

met='sarscov2s';
ihsarscov2=addMetabolite(ihsarscov2,met);
ihsarscov2.metComps(ismember(ihsarscov2.mets,met))=1; %Add to extracellular space
ihsarscov2.metNames(ismember(ihsarscov2.mets,met))={'VBOF'};

%Add transport reaction from cytosol to extracellular

ihsarscov2=addReaction(ihsarscov2,'VBOFt','reactionFormula','sarscov2c -> sarscov2s');


%Add transport reaction from extracellular to boundary

ihsarscov2=addExchangeRxn(ihsarscov2,'sarscov2s');

end