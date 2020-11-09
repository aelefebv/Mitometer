function NN = getNN(stats)

numMito = length(stats);

%initialize
xTrack = zeros(1,numMito);
yTrack = zeros(1,numMito);
xDiff = zeros(1,numMito);
yDiff = zeros(1,numMito);

for mitoNum = 1:numMito
    xTrack(mitoNum) = stats(mitoNum).WeightedCentroid(1);
    yTrack(mitoNum) = stats(mitoNum).WeightedCentroid(2);
end

for mitoNum = 1:numMito
    xDiff(mitoNum,:) = xTrack'-stats(mitoNum).WeightedCentroid(1);
    yDiff(mitoNum,:) = yTrack'-stats(mitoNum).WeightedCentroid(2);
end

distDiff = sqrt((xDiff.^2)+(yDiff.^2));
distDiff(distDiff==0) = Inf;

NN = min(distDiff);


end