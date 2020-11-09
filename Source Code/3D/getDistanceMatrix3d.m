function distanceMatrix = getDistanceMatrix3d(mito,track)

numMito = length(mito);
numTracks = length(track);

%initialize
xTrack = zeros(1,numMito);
yTrack = zeros(1,numMito);
zTrack = zeros(1,numMito);

xDiff = zeros(numMito,numTracks);
yDiff = zeros(numMito,numTracks);
zDiff = zeros(numMito,numTracks);

%In the case of distance differences, we take the differences in x, y, and
%z (in 3D) intensity weighted centroid positions. We then use the square
%root of the squared sums as the distance.
for mitoNum = 1:numMito
    xTrack(mitoNum) = mito(mitoNum).WeightedCentroid(1);
    yTrack(mitoNum) = mito(mitoNum).WeightedCentroid(2);
    zTrack(mitoNum) = mito(mitoNum).WeightedCentroid(3);
end

for trackNum = 1:numTracks
    xDiff(:,trackNum) = xTrack'-track(trackNum).WeightedCentroid(end,1);
    yDiff(:,trackNum) = yTrack'-track(trackNum).WeightedCentroid(end,2);
    zDiff(:,trackNum) = zTrack'-track(trackNum).WeightedCentroid(end,3);

end

distanceMatrix = sqrt((xDiff.^2)+(yDiff.^2)+(zDiff.^2));


end