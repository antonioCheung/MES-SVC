function groups = randomGroup(dimension, gamma, realDiv)
GS = ceil((dimension-size(realDiv,2)) / gamma);
d = 1: dimension;
d = d(randperm(length(d)));
for j = 1: size(realDiv,2)
    d(d==realDiv(j))=[];
end
groups = {};
for i = 1: gamma
    begin_index = (i - 1) * GS + 1;
    end_index = min(i * GS, dimension-size(realDiv,2));
    groups{i} = d(begin_index: end_index);
end
end