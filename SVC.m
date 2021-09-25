function [P,d] = SVC(Global,t)
    
    M = Global.M;
    Np = Global.N;
    D = Global.D;
    
    
    %%DataGeneration
    [P,DS] = DataGeneration(Np,t);

    %%Data set division
    V(1:M,:) = eye(M);
    
    empty_probmodel.dd=[];
    probmodel = repmat(empty_probmodel,M,1);
    [~,replace]=sort(pdist2((DS.objs-repmat(min(DS.objs),size(DS,2),1))./(repmat(max(DS.objs)-min(DS.objs),size(DS,2),1)),V,'cosine'),2);
    u = replace(:,1);
    for i  = 1 : M
        x = u==i;
        PP = DS(x);
        if size(PP,2)>1
            probmodel(i).dd = (max(PP.decs)-min(PP.decs))./(Global.upper-Global.lower);
        end
    end
    
    %%Variable property classification
    H1 = [];
    for i = 1 : M
        s = probmodel(i).dd;
        H1 = [H1;s];
    end
    [~,d1] = sort(min(H1));
    d = d1(1:M-1);
end

