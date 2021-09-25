function [Index,Population] = AngleSelection(Population,N)
% The environmental selection of VaEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Association operation
    Next   = [];
    Choose = Association([],Population.objs,N);
    Index = find(Choose==1);
    % Population for next generation
    Population = Population(Choose);
end

function Choose = Association(PopObj1,PopObj2,N)
% Association operation in the algorithm

    [N1,~] = size(PopObj1);
    [N2,M] = size(PopObj2);
    PopObj = [PopObj1;PopObj2];

    %% Normalization
    Zmin   = min(PopObj,[],1);
    Zmax   = max(PopObj,[],1);
    PopObj = (PopObj-repmat(Zmin,size(PopObj,1),1))./repmat(Zmax-Zmin,size(PopObj,1),1);
    
    
    %% Angle between each two solutions
    angle = acos(1-pdist2(PopObj,PopObj,'cosine'));
    
    %% Niching
    Choose = [true(1,N1),false(1,N2)];
    if ~any(Choose)
        % Select the extreme solutions first
        [~,extreme]        = min(pdist2(PopObj2,eye(M),'cosine'),[],1);
        Choose(N1+extreme) = true;
%         % Select the first M best converged solutions
%         [~,rank] = sort(fit(N1+1:end));
%         Choose(N1+rank(1:min(M,length(rank)))) = true;
    end
    while sum(Choose) < N
        % Maximum vector angle first
        Select  = find(Choose);
        Remain  = find(~Choose);
        [~,rho] = max(min(angle(Remain,Select),[],2));
        Choose(Remain(rho)) = true;
    end
end