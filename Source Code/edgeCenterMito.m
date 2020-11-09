function [edgeMito,centerMito] = edgeCenterMito(track,mitoCoM,numFrames,is3D)

numTracks = length(track);

centroidDistance = nan(numTracks,numFrames);

for frameNum = 1:numFrames
    for trackNum = 1:numTracks
        index = find(track(trackNum).frame==frameNum);
        if sum(index)
            if is3D
                trackCentroid = track(trackNum).WeightedCentroid((index*3-2):(index*3));
            else
                trackCentroid = track(trackNum).WeightedCentroid((index*2-1):(index*2));
            end
            centroidDistance(trackNum,frameNum) = norm(mitoCoM(frameNum).centroid - trackCentroid);
        end
    end
end

maxCentroidDistance = max(centroidDistance,[],2);
normCent = normalize(maxCentroidDistance,'range');
lengthHisto = round(length(normCent)*0.75);
mitoHisto = histcounts(normCent,(min(normCent):(max(normCent)-min(normCent))/lengthHisto:max(normCent)));
thresh = otsuthresh(mitoHisto);

edgeMito = find(normCent>thresh)';
centerMito = find(normCent<=thresh)';

end