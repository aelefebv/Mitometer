function [MSD,MSDstd] = getMSD(trackMSD,micronPerPixel)

numFrames = length(trackMSD.frame);
maxDeltaT = trackMSD.frame(end)-trackMSD.frame(1);

MSDmatrix = nan(numFrames,maxDeltaT);

for frameNum = 1:numFrames
    coordsCurrent = [trackMSD.WeightedCentroid(frameNum*2-1),trackMSD.WeightedCentroid(frameNum*2)];
    for deltaT = 1:maxDeltaT
        index = find(trackMSD.frame==trackMSD.frame(frameNum)+deltaT);
        if ~sum(index)
            continue
        end
        coordsDeltaT = [trackMSD.WeightedCentroid((index)*2-1),trackMSD.WeightedCentroid((index)*2)];

        MSDmatrix(frameNum,deltaT) = norm(coordsDeltaT*micronPerPixel-coordsCurrent*micronPerPixel).^2;
    end
end

MSD = nanmean(MSDmatrix);
MSDstd = nanstd(MSDmatrix);

end