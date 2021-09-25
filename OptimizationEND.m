function Arc= OptimizationEND(Global,P1,P2,N)

    Arc = [];
    %% Calculate the reference directions
	lower       = Global.lower;upper = Global.upper;
    delta = upper - lower;
    decs1 = (P1.decs - lower)./delta;
    decs2 = (P2.decs - lower)./delta;
    
    
	Direction   = sum((decs1-decs2).^2,2).^(0.5);
	Direct      = (decs1-decs2)./repmat(Direction,1,Global.D);    
    
    decs = decs1;
    w = [];
    for i = 1 : Global.D
        if Direct(i)==0
            continue;
        end
        if Direct(i)>0
            w = [w,(decs(i))/Direct(i)];
        else 
            w = [w,(1-decs(i))/abs(Direct(i))];
        end
    end
    wmax1 = min(w);
    d = wmax1/N;
%     w0 = rand(N,1).*wmax1;                                % Initialize the population
    w0 = (1:N)'.*d;
    PopNew1 = fitfunc(w0,Direct,Global,decs1,delta);	% Calculate the fitness and store the solutions
  	
    wmax2        = mean(w)-min(w);
    d = wmax2/N;
    w0 = (1:N)'.*d + min(w);
%     w0 = rand(N,1).*wmax2 + min(w);
    PopNew2 = fitfunc(w0,Direct,Global,decs1,delta);
%     
%     wmax3        = max(w)-mean(w);
%     d = wmax3/N;
%     w0 = (1:N)'.*d + mean(w);
% %     w0 = rand(N,1).*wmax3 + mean(w);
%     PopNew3 = fitfunc(w0,Direct,Global,decs1,delta);

    Arc = [Arc;PopNew1;PopNew2];
    
end


%% 首先还原个体的decs，然后变为个体计算PBI值
function OffSpring = fitfunc(w0,direct,Global,P,delta)
    [SubN,~] = size(w0); 
%     Obj   	= zeros(SubN,1);
    OffSpring  = [];
    for i = 1 : SubN 
        PopDec  = Global.lower + (P - repmat(w0(i),1,Global.D).*direct).*delta;
        PopDec(PopDec<Global.lower) = Global.lower(PopDec<Global.lower);
        PopDec(PopDec>Global.upper) = Global.upper(PopDec>Global.upper);
%         OffWPop = INDIVIDUAL(PopDec);
%         Obj(i) = sum(OffWPop.objs);
        OffSpring = [OffSpring;PopDec];
    end
end

