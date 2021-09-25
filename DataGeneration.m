function [P,DS] = DataGeneration(Np,t)

    Global  = GLOBAL.GetObj();
    P= Global.Initialization();
    count_off = floor(Np/2)*2;
    i = 1;
    start = (i-1)*Np+1;
    last = i*Np;
    DS(start:last) = P;
    i = i + 1;
    while Global.evaluated < t
        MatingPool = randi(Np,1,Np);
        Offspring  = GA(P(MatingPool)); 
        start = (i-1)*count_off+1;
        last = i*count_off;
        DS(start:last) = Offspring;
        i = i + 1;
        P= VaEAEnvironmentalSelection([P,Offspring],Np);
    end
end

