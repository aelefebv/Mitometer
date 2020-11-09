function extremaMatrix = getExtremaMatrixFission3d(newTrack,track,frameNum)

numTracks = length(track);

%initialize
extremaMatrix = zeros(1,numTracks);

for trackNum = 1:numTracks
    frameIdx = find(track(trackNum).frame==frameNum);
    if ~frameIdx
        extremaMatrix(1,trackNum) = Inf;
    else
        xExtrema = track(trackNum).ConvexHull{frameIdx}(:,1);
        yExtrema = track(trackNum).ConvexHull{frameIdx}(:,2);
        zExtrema = track(trackNum).ConvexHull{frameIdx}(:,3);
        DistMatX = inf(size(track(trackNum).ConvexHull{frameIdx},1),size(newTrack.ConvexHull{1},1));
        DistMatY = inf(size(track(trackNum).ConvexHull{frameIdx},1),size(newTrack.ConvexHull{1},1));
        DistMatZ = inf(size(track(trackNum).ConvexHull{frameIdx},1),size(newTrack.ConvexHull{1},1));
        for newConvexNum = 1:size(newTrack.ConvexHull,1)
            DistMatX(:,newConvexNum) = xExtrema - newTrack.ConvexHull{1}(newConvexNum,1);
            DistMatY(:,newConvexNum) = yExtrema - newTrack.ConvexHull{1}(newConvexNum,2);
            DistMatZ(:,newConvexNum) = zExtrema - newTrack.ConvexHull{1}(newConvexNum,3);
        end
        DistMat = sqrt(DistMatX.^2+DistMatY.^2+DistMatZ.^2);
        
%         DistMat = inf(size(track(trackNum).ConvexHull{frameIdx},1),size(newTrack.ConvexHull{1},1));
%         for numPoly1 = 1:size(track(trackNum).ConvexHull{frameIdx},1)
%             for numPoly2 = 1:size(newTrack.ConvexHull{1},1)
%                 DistMat(numPoly1,numPoly2) = norm(track(trackNum).ConvexHull{frameIdx}(numPoly1,:)-newTrack.ConvexHull{1}(numPoly2,:));
%             end
%         end
        extremaMatrix(1,trackNum) = min(DistMat,[],'all');    
    end
end
end
