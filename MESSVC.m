function MESSVC(Global)
% <algorithm> <M>

    %% Statistical variable classification and initialization
    t = 10000; Nc = Global.M+1; T = ceil(Global.N/10); gamma = 4; Ns = 10; threshold = 0.03;
    [P,d] = SVC(Global,t);

    while Global.NotTermination(P)
        Old = sum(sum(P.objs,2));
        %%Memetic exploitation and exploration
        [xPrimeList,P,u1] = MEE(Global,P,Nc,T,d,gamma);
        %%Memetic elite imitation
        P = MEI(Global,xPrimeList,P,u1,d,Ns);
        New = sum(sum(P.objs,2));
        if Old-New<threshold && Global.evaluated<0.6*Global.evaluation 
            %%Diversity optimization strategy
            for i = 1 : Global.M
                Z = min(P.objs);
                Znad = max(P.objs);
                PopObj = P.objs;
                PopObj = (PopObj-repmat(Z,size(PopObj,1),1))./repmat(Znad-Z,size(PopObj,1),1);
                [~,re1] = sort(pdist2(PopObj,PopObj,'cosine'),2);
                B = re1(:,1:T);
                re2=flip(re1,2);
                B = [B,re2(:,1:T)];
                P = DistributionOptimization(P,d,B,T,Global.N);
            end
        end
    end
end



 function Population = DistributionOptimization(Population,PV,B,T,N)
    NewP = [];
    for i = 1 : size(Population,2)
        P = Population(B(i,:));
        MatingPool = randi(T,1,T);
        OffDec = P.decs;
        NewDec       = GA(P(MatingPool).decs);
        OffDec(1:size(NewDec,1),PV) = NewDec(:,PV);
        NewP = [NewP,INDIVIDUAL(OffDec)];
    end
    Population   = VaEAEnvironmentalSelection([Population,NewP],N);
 end