function distanceMatrix = getDistanceMatrix(mito,track)

numMito = length(mito);
numTracks = length(track);

%initialize
xTrack = zeros(1,numMito);
yTrack = zeros(1,numMito);
xDiff = zeros(numMito,numTracks);
yDiff = zeros(numMito,numTracks);

for mitoNum = 1:numMito
    xTrack(mitoNum) = mito(mitoNum).WeightedCentroid(1);
    yTrack(mitoNum) = mito(mitoNum).WeightedCentroid(2);
end

for trackNum = 1:numTracks
    xDiff(:,trackNum) = xTrack'-track(trackNum).WeightedCentroid(end,1);
    yDiff(:,trackNum) = yTrack'-track(trackNum).WeightedCentroid(end,2);
end

distanceMatrix = sqrt((xDiff.^2)+(yDiff.^2));


end