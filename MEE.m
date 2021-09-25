function [xPrimeList,P,u1] = MEE(Global,P,Nc,T,d,gamma)

    
    Np = Global.N;
    M = Global.M;
    D = Global.D;

    %% Niche memetic exploitation strategy
    %obtain the reference points
    [~,xPrimeList] = AngleSelection(P,Nc); 
    
    %Normalization
    Z = min(P.objs); Znad = max(P.objs);
    PopObj = P.objs;
    PopObj = (PopObj-repmat(Z,size(PopObj,1),1))./repmat(Znad-Z,size(PopObj,1),1);
    xPrimeListObj = xPrimeList.objs;
    xPrimeListObj = (xPrimeListObj-repmat(Z,size(xPrimeListObj,1),1))./repmat(Znad-Z,size(xPrimeListObj,1),1);
    
    %Association and Randomly select one individual from each cluster
    [~,re] = sort(pdist2(PopObj,xPrimeListObj,'cosine'),2);
    u1 = re(:,1);
    Index = zeros(1,Nc);
    xPrimeList = repmat(INDIVIDUAL,1,Nc);
    for i = 1 :  Nc
        x = find(u1==i);
        PP= P(x);
        a=size(PP,2);
        r=randi(a);
        xPrimeList(i) = PP(r);
        Index(i) = x(r);
    end

    %Association
    PP= repmat(INDIVIDUAL,1,T*size(xPrimeList,2));
    for i = 1 : size(xPrimeList,2)
         xPrimeListObj = xPrimeList(i).objs;
         xPrimeListObj = (xPrimeListObj-repmat(Z,size(xPrimeListObj,1),1))./repmat(Znad-Z,size(xPrimeListObj,1),1);
         [~,re] = sort(pdist2(PopObj,xPrimeListObj,'cosine'));
         start = (i-1)*T+1;
         last = i*T;
         PP(start:last) = P(re(1:T));
    end
    
    %Exploitation
    bestVal = zeros(Nc,D);
    for i = 1 : Nc
        start = (i-1)*T+1;
        last = i*T;
        PopDec = PP(start:last).decs;
        dec = xPrimeList(i).dec;
        PopDec(:,d) = repmat(dec(d),T,1);
        Pop = INDIVIDUAL(PopDec);
        [~,re] = min(sum(Pop.objs,2));
        bestVal(i,:) = Pop(re).dec;
    end
    
    
    
    %% Associative memetic exploration strategy
    decs = xPrimeList.decs;
    divs = decs(:,d);
    n_task = size(xPrimeList,2);
    
    TaskIndividual.Population = INDIVIDUAL();
    TaskIndividual.skill = [];
    TaskIndividual.rank = [];
    
    %RandomGroup
    groups = randomGroup(Global.D, gamma, d);
   
    
    for gg = 1 : gamma
        %SmallPopuConstruct
        TaskPopulation = repmat(TaskIndividual,1,n_task*T);
        dim = groups{gg};
        for i = 1 : n_task
            PopDec = repmat(bestVal(i,:),T-1,1);
            PopDec2 = unifrnd(repmat(Global.lower(dim),T-1,1),repmat(Global.upper(dim),T-1,1));
            dec = xPrimeList(i).dec;
            PopDec(:,dim) = PopDec2;
            PopDec = [PopDec;dec];
            PopDec(:,d) = repmat(dec(d),T,1);
            PP = INDIVIDUAL(PopDec);
            [~,re] = sort(sum(PP.objs,2));
            taskP = repmat(TaskIndividual,1,T);
            for j = 1 : T
                taskP(j).Population = PP(j);
                taskP(j).skill      = i;
                taskP(j).rank       = re(j);
            end
            start = (i-1)*T+1;
            last = i*T;
            TaskPopulation(start:last) = taskP;
        end

        %RandomGroup
        r_task = randomGroup1(n_task, 3);
        
     
        RTaskPopulation = TaskPopulation;
        for r = 1 : size(r_task,2)
            d1 = r_task{r};
            flag = zeros(1,size(RTaskPopulation,2));
            for j = 1 : size(d1,2)
                flag = flag | [RTaskPopulation.skill]==d1(j);
            end
            TaskPopulation = RTaskPopulation(flag);
            
            %AssortativeMating
            rmp = 0.3;
            pop = size(TaskPopulation,2);

            for generation=1:size(dim,2)
                rndlist=randperm(pop);
                TaskPopulation=TaskPopulation(rndlist);
                parent = repmat(TaskIndividual,1,pop);
                for i = 1:pop 
                    p1=1+round(rand(1)*(pop-1));
                    p2=1+round(rand(1)*(pop-1));
                    if TaskPopulation(p1).rank < TaskPopulation(p2).rank
                        parent(i).Population=TaskPopulation(p1).Population;
                        parent(i).skill=TaskPopulation(p1).skill;
                    elseif TaskPopulation(p1).rank == TaskPopulation(p2).rank
                        if rand(1) <= 0.5
                            parent(i).Population=TaskPopulation(p1).Population;
                            parent(i).skill=TaskPopulation(p1).skill;
                        else
                            parent(i).Population=TaskPopulation(p2).Population;
                            parent(i).skill=TaskPopulation(p2).skill;
                        end
                    else
                        parent(i).Population=TaskPopulation(p2).Population;
                        parent(i).skill=TaskPopulation(p2).skill;
                    end
                end

                %VerticalCulturalTransmission
                count=1;
                child = repmat(TaskIndividual,1,pop);
                for i=1:2:pop-1 % Create offspring population via mutation and crossover
                    p1=i;
                    p2=i+1;
                    if parent(p1).skill==parent(p2).skill
                        P1 = parent(p1).Population; dec1 = P1.decs;
                        P2 = parent(p2).Population; dec2 = P2.decs;
                        OffDec = GA34([dec1(dim);dec2(dim)],dim);
                        child(count).skill = parent(p1).skill;
                        child(count+1).skill = parent(p1).skill; 
                        OffDec2 = [bestVal(child(count).skill,:);bestVal(child(count+1).skill,:)];
                        OffDec2(:,dim) = OffDec;
                        OffDec = OffDec2;
                        OffDec(:,d) = repmat(dec1(d),2,1);
                        Offspring = INDIVIDUAL(OffDec);
                        child(count).Population = Offspring(1);
                        child(count+1).Population = Offspring(2);

                    else
                        if rand(1)<rmp
                            P1 = parent(p1).Population; dec1 = P1.decs;
                            P2 = parent(p2).Population; dec2 = P2.decs;
                            OffDec = GA34([dec1(dim);dec2(dim)],dim);
                            child(count).skill = d(randi(size(d,2)));
                            child(count+1).skill = d(randi(size(d,2)));
                            OffDec2 = [bestVal(child(count).skill,:);bestVal(child(count+1).skill,:)];
                            OffDec2(:,dim) = OffDec;
                            OffDec = OffDec2;
                            OffDec(:,d) = [divs(child(count).skill,:);divs(child(count+1).skill,:)];
                            Offspring = INDIVIDUAL(OffDec);
                            child(count).Population = Offspring(1);
                            child(count+1).Population = Offspring(2);
                        else
                            P1 = parent(p1).Population; dec1 = P1.decs;
                            P2 = parent(p2).Population; dec2 = P2.decs;
                            OffDec = GA24([dec1(dim);dec2(dim)],dim);
                            child(count).skill = parent(p1).skill;
                            child(count+1).skill = parent(p2).skill;
                            OffDec2 = [bestVal(child(count).skill,:);bestVal(child(count+1).skill,:)];
                            OffDec2(:,dim) = OffDec;
                            OffDec = OffDec2;
                            OffDec(:,d) = [divs(child(count).skill,:);divs(child(count+1).skill,:)];
                            Offspring = INDIVIDUAL(OffDec);
                            child(count).Population = Offspring(1);
                            child(count+1).Population = Offspring(2);  
                        end
                    end
                    count=count+2;
                end

                %Select best T individuals for each selected task;
                Comibine = [TaskPopulation,child];
                newP = repmat(TaskIndividual,1,pop);
                for k = 1 : size(d1,2)
                    intpopulation = Comibine([Comibine.skill]==d1(k));
                    F = zeros(1,size(intpopulation,2));
                    for m = 1 : size(intpopulation,2)
                        p = intpopulation(m).Population;
                        F(m) = sum(p.obj); 
                    end
                    [~,re] = sort(F);
                    intpopulation = intpopulation(re);
                    for m = 1 : size(intpopulation,2)
                        intpopulation(m).rank = m;
                    end
                    start = (k-1)*T+1;
                    last = k*T;
                    newP(start:last) = intpopulation(1:T);
                end
                for k = 1 : pop     
                     TaskPopulation(k) = newP(k);
                end

            end
            
            
            for j = 1 : size(d1,2)
                flag = [RTaskPopulation.skill]==d1(j);
                flag1 = [TaskPopulation.skill]==d1(j);
                RTaskPopulation(flag) = TaskPopulation(flag1);
            end

        end
       TaskPopulation = RTaskPopulation;

       for s = 1 : n_task
           PP = TaskPopulation((s-1)*T+1).Population;
           if sum(PP.obj)<sum(xPrimeList(s).obj)
               xPrimeList(s) = PP;
               bestVal(s,:) = PP.dec;
           end
       end
    end

    %For each task, keep a best individual
    PP= repmat(INDIVIDUAL,1,size(xPrimeList,2));
    for i = 1 : size(xPrimeList,2)
        p= TaskPopulation((i-1)*T+1).Population;
        PP(i) = p;
    end
    xPrimeList=PP;
    P(Index) = PP;  



end

