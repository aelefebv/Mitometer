function NN = getNN3d(stats)

numMito = size(stats,1);

%initialize
xTrack = zeros(1,numMito);
yTrack = zeros(1,numMito);
zTrack = zeros(1,numMito);
xDiff = zeros(1,numMito);
yDiff = zeros(1,numMito);
zDiff = zeros(1,numMito);

for mitoNum = 1:numMito
    xTrack(mitoNum) = stats.WeightedCentroid(mitoNum,1);
    yTrack(mitoNum) = stats.WeightedCentroid(mitoNum,2);
    zTrack(mitoNum) = stats.WeightedCentroid(mitoNum,3);   
end

for mitoNum = 1:numMito
    xDiff(mitoNum,:) = xTrack'-stats.WeightedCentroid(mitoNum,1);
    yDiff(mitoNum,:) = yTrack'-stats.WeightedCentroid(mitoNum,2);
    zDiff(mitoNum,:) = zTrack'-stats.WeightedCentroid(mitoNum,3);
end

distDiff = sqrt((xDiff.^2)+(yDiff.^2)+(zDiff.^2));
distDiff(distDiff==0) = Inf;

NN = min(distDiff);


end