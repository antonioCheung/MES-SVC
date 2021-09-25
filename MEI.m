function P = MEI(Global,xPrimeList,P,u1,d,Ns)

    AllP = [];
    re = [];
    Nc = size(xPrimeList,1);
    for i = 1 : Nc
        x = find(u1==i);
        re = [re,x'];
        PP = P(x);
        for j = 1 : size(PP,2)
            p1 = PP(j);
            p2 = xPrimeList(i);
            Arc= OptimizationEND(Global,p1,p2,Ns);
            Arc2= OptimizationEND(Global,p2,p1,10);
            Arc2(:,d) = repmat(p1.dec(d),20,1);
            NewP = INDIVIDUAL([Arc;Arc2]);
            [~,r] = min(sum(NewP.objs,2));
            P1 = NewP(r);
            AllP = [AllP,P1];     
        end
    end
    P(re) = AllP;
end

