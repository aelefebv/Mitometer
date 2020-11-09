function extremaMatrix = getExtremaMatrixFission(newTrack,track,frameNum)

numTracks = length(track);

%initialize
extremaMatrix = zeros(1,numTracks);

for trackNum = 1:numTracks
    frameIdx = find(track(trackNum).frame==frameNum);
    if ~frameIdx
        extremaMatrix(1,trackNum) = Inf;
    else
        xExtrema = track(trackNum).Extrema(:,frameIdx*2-1);
        yExtrema = track(trackNum).Extrema(:,frameIdx*2);
        DistMatX = inf(size(track(trackNum).Extrema,1),size(newTrack.Extrema,1));
        DistMatY = inf(size(track(trackNum).Extrema,1),size(newTrack.Extrema,1));
        for newExtremaNum = 1:size(newTrack.Extrema,1)
            DistMatX(:,newExtremaNum) = xExtrema - newTrack.Extrema(newExtremaNum,1);
            DistMatY(:,newExtremaNum) = yExtrema - newTrack.Extrema(newExtremaNum,2);
        end
        DistMat = sqrt(DistMatX.^2+DistMatY.^2);
%         
%         for numExtrema1 = 1:size(track(trackNum).Extrema,1)
%             for numExtrema2 = 1:size(newTrack.Extrema,1)
%                 DistMat(numExtrema1,numExtrema2) = norm((track(trackNum).Extrema(numExtrema1,frameIdx*2-1:frameIdx*2))-newTrack.Extrema(numExtrema2,1:2));
%             end
%         end
        extremaMatrix(1,trackNum) = min(DistMat,[],'all');    
    end
end
end
