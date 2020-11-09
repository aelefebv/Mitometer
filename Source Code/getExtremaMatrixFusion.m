function extremaMatrix = getExtremaMatrixFusion(lostTrack,track,frameNum)

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
        DistMatX = inf(size(track(trackNum).Extrema,1),size(lostTrack.Extrema,1));
        DistMatY = inf(size(track(trackNum).Extrema,1),size(lostTrack.Extrema,1));
        for lostExtremaNum = 1:size(lostTrack.Extrema,1)
            DistMatX(:,lostExtremaNum) = xExtrema - lostTrack.Extrema(lostExtremaNum,1);
            DistMatY(:,lostExtremaNum) = yExtrema - lostTrack.Extrema(lostExtremaNum,2);
        end
        DistMat = sqrt(DistMatX.^2+DistMatY.^2);   
%         
%         for numExtrema1 = 1:size(track(trackNum).Extrema,1)
%             for numExtrema2 = 1:size(lostTrack.Extrema,1)
%                 DistMat(numExtrema1,numExtrema2) = norm((track(trackNum).Extrema(numExtrema1,frameIdx*2-1:frameIdx*2))-lostTrack.Extrema(numExtrema2,end-1:end));
%             end
%         end
        extremaMatrix(1,trackNum) = min(DistMat,[],'all');    
    end
end
end
