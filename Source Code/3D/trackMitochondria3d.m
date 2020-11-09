function [track,mitoCoM,extra] = trackMitochondria3d(Im,ImOG,weights,micronPerPixel,secondsPerFrame,distThreshMicrons,frameThreshSeconds,zDistanceMicrons,stdRange,skipDyn)

numFrames = size(Im,3);

zRatio = zDistanceMicrons/micronPerPixel;

% User Parameters
distThresh = (distThreshMicrons/micronPerPixel)*secondsPerFrame;
frameThresh = round(frameThreshSeconds/secondsPerFrame);
% Other Parameters
wVolume = weights(1);
wMajAx = weights(2);
wMinAx = weights(3);
wSol = weights(4);
wSurf = weights(5);
wInt = weights(6);
wZAx = weights(7);

startCost = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Assign Tracks
% Initialize
extra = struct();
CC = cell(1,numFrames);
L = cell(1,numFrames);
mito = cell(1,numFrames);
stats = cell(1,numFrames);
numMito = zeros(1,numFrames);
NN = cell(1,numFrames);
NNExtrema = cell(1,numFrames);
firstUsableFrame = true;
secondUsableFrame = false;
mitoCoM = struct();
distanceMatrix = cell(1,numFrames);
volumeMatrix = cell(1,numFrames);
majAxMatrix = cell(1,numFrames);
minAxMatrix = cell(1,numFrames);
zAxMatrix = cell(1,numFrames);
solMatrix = cell(1,numFrames);
surfMatrix = cell(1,numFrames);
intMatrix = cell(1,numFrames);
frameMatrix = cell(1,numFrames);
NPA = cell(1,numFrames);
distanceMatrixThresholded = cell(1,numFrames);
distanceMatrixSq = cell(1,numFrames);
volumeMatrixSq = cell(1,numFrames);
majAxMatrixSq = cell(1,numFrames);
minAxMatrixSq = cell(1,numFrames);
zAxMatrixSq = cell(1,numFrames);
solMatrixSq = cell(1,numFrames);
surfMatrixSq = cell(1,numFrames);
intMatrixSq = cell(1,numFrames);
distanceMatrixZ = cell(1,numFrames);
volumeMatrixZ = cell(1,numFrames);
majAxMatrixZ = cell(1,numFrames);
minAxMatrixZ = cell(1,numFrames);
zAxMatrixZ = cell(1,numFrames);
solMatrixZ = cell(1,numFrames);
surfMatrixZ = cell(1,numFrames);
intMatrixZ = cell(1,numFrames);
differenceMatrix = cell(1,numFrames);
assignTrackMatrix = cell(1,numFrames);
matchedTrack = cell(1,numFrames);
confident = cell(1,numFrames);


firstFrame = 1;

runningConfidentCosts = [];

loading = waitbar(0,'Please wait...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(loading,'canceling',0);
pause(.0002)
frameIm = zeros(size(Im,1),size(Im,2),size(Im,4));

for frameNum = 1:numFrames
    %
    if getappdata(loading,'canceling')
        delete(loading)
        break
    end
    waitbar(frameNum/numFrames,loading,sprintf('Tracking frame %d of %d.',frameNum,numFrames));
    pause(.0002)
    %
    
    frameIm(:,:,:) = Im(:,:,frameNum,:);
    %Get stats for each mitochondrion in each frame.
    CC{frameNum} = bwconncomp(frameIm);
    L{frameNum} = labelmatrix(CC{frameNum});
    stats{frameNum} = regionprops3(CC{frameNum},frameIm,'ConvexHull','BoundingBox','Orientation','PrincipalAxisLength','Solidity','SurfaceArea','Volume','VoxelIdxList','MeanIntensity','WeightedCentroid');
    numMito(frameNum) = size(stats{frameNum},1);
    
    for mitoNum = 1:numMito(frameNum)
        mito{frameNum}(mitoNum).BoundingBox = stats{frameNum}.BoundingBox(mitoNum,:);
        mito{frameNum}(mitoNum).Orientation = stats{frameNum}.Orientation(mitoNum,:);
        mito{frameNum}(mitoNum).MajorAxisLength = stats{frameNum}.PrincipalAxisLength(mitoNum,1);
        mito{frameNum}(mitoNum).MinorAxisLength = stats{frameNum}.PrincipalAxisLength(mitoNum,2);
        mito{frameNum}(mitoNum).ZAxisLength = stats{frameNum}.PrincipalAxisLength(mitoNum,3);
        mito{frameNum}(mitoNum).Solidity = stats{frameNum}.Solidity(mitoNum,:);
        mito{frameNum}(mitoNum).SurfaceArea = stats{frameNum}.SurfaceArea(mitoNum,:);
        mito{frameNum}(mitoNum).ConvexHull = stats{frameNum}.ConvexHull(mitoNum,:);
        mito{frameNum}(mitoNum).Volume = stats{frameNum}.Volume(mitoNum,:);
        mito{frameNum}(mitoNum).VoxelIdxList = stats{frameNum}.VoxelIdxList(mitoNum,:);
        mito{frameNum}(mitoNum).MeanIntensity = stats{frameNum}.MeanIntensity(mitoNum,:);
        mito{frameNum}(mitoNum).WeightedCentroid = stats{frameNum}.WeightedCentroid(mitoNum,:);
    end
    
    %Keep iterating through frames until one with mitochondria appears.
    if ~sum(numMito)
        firstFrame = firstFrame+1;
        continue
    end
    
    %Get the nearest neighbor distances based on intensity weighted
    %centroids and based on extrema.
    NN{frameNum} = getNN3d(stats{frameNum});
%     NNExtrema{frameNum} = getNNExtrema3d(mito{frameNum});
    
    %Correct the angles of the mitochondria and assign NN, frame, and
    %initiliaze other parameters.
    for mitoNum = 1:numMito(frameNum)
%         mito{frameNum}(mitoNum).MaxFeretAngle = angle0to180(180-mito{frameNum}(mitoNum).MaxFeretAngle);
        mito{frameNum}(mitoNum).NN = NN{frameNum}(mitoNum);
%         mito{frameNum}(mitoNum).NNExtrema = NNExtrema{frameNum}(mitoNum);
        mito{frameNum}(mitoNum).frame = frameNum;
        mito{frameNum}(mitoNum).NPA = 0;
        mito{frameNum}(mitoNum).confident = 1;
        mito{frameNum}(mitoNum).lost = 0;
        mito{frameNum}(mitoNum).fission = 0;
        mito{frameNum}(mitoNum).fusion = 0;
        mito{frameNum}(mitoNum).label = mitoNum;
    end

    if firstUsableFrame
        firstUsableFrame = false;
        secondUsableFrame = true;
        
        %Create a track from all mito in the first frame
        track = mito{frameNum};
        
        %Keep track of mito center of mass in the first frame
        tempOG(:,:,:) = Im(:,:,frameNum,:);
        binaryOG = true(size(tempOG));
        labeledOG = logical(binaryOG);
        CoMstats = regionprops3(labeledOG,tempOG,'WeightedCentroid');
        mitoCoM(frameNum).centroid = CoMstats.WeightedCentroid;
        mitoCoM(frameNum).x = CoMstats.WeightedCentroid(1);
        mitoCoM(frameNum).y = CoMstats.WeightedCentroid(2);
        mitoCoM(frameNum).z = CoMstats.WeightedCentroid(3);
   
        continue
    end
        
    %Keep track of mito center of mass
    tempOG(:,:,:) = Im(:,:,frameNum,:);
    binaryOG = true(size(tempOG));
    labeledOG = logical(binaryOG);
    CoMstats = regionprops3(labeledOG,tempOG,'WeightedCentroid');
    mitoCoM(frameNum).centroid = CoMstats.WeightedCentroid;
    mitoCoM(frameNum).x = CoMstats.WeightedCentroid(1);
    mitoCoM(frameNum).y = CoMstats.WeightedCentroid(2);
    mitoCoM(frameNum).z = CoMstats.WeightedCentroid(3);
   
    %Get the costs for the cost matrix
    distanceMatrix{frameNum} = getDistanceMatrix3d(mito{frameNum},track);
    volumeMatrix{frameNum} = getDeltaMatrix(mito{frameNum},track,'Volume');
    majAxMatrix{frameNum} = getDeltaMatrix(mito{frameNum},track,'MajorAxisLength');
    minAxMatrix{frameNum} = getDeltaMatrix(mito{frameNum},track,'MinorAxisLength');
    zAxMatrix{frameNum} = getDeltaMatrix(mito{frameNum},track,'ZAxisLength');
    solMatrix{frameNum} = getDeltaMatrix(mito{frameNum},track,'Solidity');
    surfMatrix{frameNum} = getDeltaMatrix(mito{frameNum},track,'SurfaceArea');
    intMatrix{frameNum} = getDeltaMatrix(mito{frameNum},track,'MeanIntensity');
    
    %Get the search time threshold matrix
    frameMatrix{frameNum} = getFrameMatrix(mito{frameNum},track);
    
    %Threshold out distances if it surpasses the max velocity or search time
    distanceMatrixThresholded{frameNum} = distanceMatrix{frameNum};
    
    distanceMatrixThresholded{frameNum}(distanceMatrix{frameNum}>distThresh) = Inf;
    distanceMatrixThresholded{frameNum}(frameMatrix{frameNum}>frameThresh) = Inf;
    
    %Get number of possible assignments for each mito.
    NPA{frameNum} = sum(~isinf(distanceMatrixThresholded{frameNum}),2);
    for mitoNum = 1:numMito(frameNum)
        mito{frameNum}(mitoNum).NPA = NPA{frameNum}(mitoNum);
    end
    
    
    %Square the difference matrices
    distanceMatrixSq{frameNum} = distanceMatrixThresholded{frameNum}.^2;
    volumeMatrixSq{frameNum} = volumeMatrix{frameNum}.^2;
    majAxMatrixSq{frameNum} = majAxMatrix{frameNum}.^2;
    minAxMatrixSq{frameNum} = majAxMatrix{frameNum}.^2;
    zAxMatrixSq{frameNum} = majAxMatrix{frameNum}.^2;
    solMatrixSq{frameNum} = solMatrix{frameNum}.^2;
    surfMatrixSq{frameNum} = majAxMatrix{frameNum}.^2;
    intMatrixSq{frameNum} = intMatrix{frameNum}.^2;
    
%     solMatrixSq{frameNum}(isinf(solMatrixSq{frameNum})) = max(solMatrixSq{frameNum}(~isinf(solMatrixSq{frameNum})));
%     solMatrixSq{frameNum}(isnan(solMatrixSq{frameNum})) = max(solMatrixSq{frameNum}(~isnan(solMatrixSq{frameNum})));

    %Z-score of the difference matrices to normalize variables
    distanceMatrixZ{frameNum} = distanceMatrixSq{frameNum};
    volumeMatrixZ{frameNum} = distanceMatrixSq{frameNum};
    majAxMatrixZ{frameNum} = distanceMatrixSq{frameNum};
    minAxMatrixZ{frameNum} = distanceMatrixSq{frameNum};
    zAxMatrixZ{frameNum} = distanceMatrixSq{frameNum};
    solMatrixZ{frameNum} = distanceMatrixSq{frameNum};
    surfMatrixZ{frameNum} = distanceMatrixSq{frameNum};
    intMatrixZ{frameNum} = distanceMatrixSq{frameNum};

    distanceMatrixZ{frameNum}(~isinf(distanceMatrixSq{frameNum})) = zscore(distanceMatrixSq{frameNum}(~isinf(distanceMatrixSq{frameNum})),0,'all');
    volumeMatrixZ{frameNum}(~isinf(distanceMatrixSq{frameNum})) = zscore(volumeMatrixSq{frameNum}(~isinf(distanceMatrixSq{frameNum})),0,'all');
    majAxMatrixZ{frameNum}(~isinf(distanceMatrixSq{frameNum})) = zscore(majAxMatrixSq{frameNum}(~isinf(distanceMatrixSq{frameNum})),0,'all');
    minAxMatrixZ{frameNum}(~isinf(distanceMatrixSq{frameNum})) = zscore(minAxMatrixSq{frameNum}(~isinf(distanceMatrixSq{frameNum})),0,'all');
    zAxMatrixZ{frameNum}(~isinf(distanceMatrixSq{frameNum})) = zscore(zAxMatrixSq{frameNum}(~isinf(distanceMatrixSq{frameNum})),0,'all');
    solMatrixZ{frameNum}(~isinf(distanceMatrixSq{frameNum})) = zscore(solMatrixSq{frameNum}(~isinf(distanceMatrixSq{frameNum})),0,'all');
    surfMatrixZ{frameNum}(~isinf(distanceMatrixSq{frameNum})) = zscore(surfMatrixSq{frameNum}(~isinf(distanceMatrixSq{frameNum})),0,'all');
    intMatrixZ{frameNum}(~isinf(distanceMatrixSq{frameNum})) = zscore(intMatrixSq{frameNum}(~isinf(distanceMatrixSq{frameNum})),0,'all');

    
    %Create a matrix for assigning mito as new tracks instead
%     newTrackMatrix = Inf(numMito(frameNum),numMito(frameNum));
%     for mitoNum = 1:numMito(frameNum)
%         newTrackMatrix(mitoNum,mitoNum) = 1;
%     end
    
    %Create a cost matrix of the z-scores and append newTrackMatrix
    differenceMatrix{frameNum} = distanceMatrixZ{frameNum} + (wVolume * volumeMatrixZ{frameNum}) + (wMajAx * majAxMatrixZ{frameNum}) + (wMinAx * minAxMatrixZ{frameNum}) + (wZAx * zAxMatrixZ{frameNum}) + (wSol * solMatrixZ{frameNum}) + (wSurf * surfMatrixZ{frameNum}) + (wInt * intMatrixZ{frameNum});
    tempDifferenceMatrix = differenceMatrix{frameNum};
    tempDifferenceMatrix(tempDifferenceMatrix==Inf) = nan;
%     newTrackCost = quantile(tempDifferenceMatrix,quantileCost,'all');
    if ~secondUsableFrame
%         startCost = quantile(runningConfidentCosts,1-quantileCost,'all'); 
        startCost = quantile(runningConfidentCosts,0.98,'all');
    end
    secondUsableFrame = false;
    if isnan(startCost)
        startCost = 1;
    end
    newTrackMatrix = diag(ones(size(tempDifferenceMatrix,1),1))*startCost;
    newTrackMatrix(newTrackMatrix==0) = Inf;

            %     for i = 1:length(newTrackMatrix)
% %         newTrackMatrix(i,i) = quantile(tempDifferenceMatrix,min(0.2,size(differenceMatrix{frameNum},2)/mito{frameNum}(i).NPA),'all');
%         newTrackMatrix(i,i) = quantile(tempDifferenceMatrix(i,:),quantileCost,'all')+0.001;%,quantile(tempDifferenceMatrix(i,:),0.9,'all'));
% 
%     end
%     newTrackMatrix(isnan(newTrackMatrix))=min(tempDifferenceMatrix(:));

    assignTrackMatrix{frameNum} = [differenceMatrix{frameNum},newTrackMatrix];
    
    %Assign mitochondria to tracks
%     frameNum
%     assignTrackMatrix{frameNum}

%     newTrackMatrix

    for assignCol = 1:size(assignTrackMatrix{frameNum},2)
        if ~sum(~isinf(assignTrackMatrix{frameNum}(:,assignCol))) || ~sum(~isnan(assignTrackMatrix{frameNum}(:,assignCol)))
            assignTrackMatrix{frameNum}(:,assignCol) = 100;
        end
    end
    matchedTrack{frameNum} = matchpairs(assignTrackMatrix{frameNum},1e9);
    
    %Find high confidence assignments, where minimum assignment value per
    %row is the same as the cost matrix minimization assignment
    [M,I] = min(assignTrackMatrix{frameNum},[],2);
    [~,idx] = sort(matchedTrack{frameNum}(:,1));
    sortedMatches = matchedTrack{frameNum}(idx,:);
    
    confident{frameNum} = (sortedMatches(:,2) == I);
    for mitoNum = 1:numMito(frameNum)
        mito{frameNum}(mitoNum).confident = confident{frameNum}(mitoNum);
    end
    allConfidentCosts = M(confident{frameNum});
    nonNewConfidentCosts = allConfidentCosts(I(confident{frameNum})<length(track));
    runningConfidentCosts = [runningConfidentCosts;nonNewConfidentCosts];
    %add mito to tracks / make new tracks
    tempTrack = track;
    for matchNum = 1:length(matchedTrack{frameNum}(:,1))
        mitoMatch = sortedMatches(matchNum,1);
        trackMatch = sortedMatches(matchNum,2);
        if trackMatch <= length(tempTrack)
            track(trackMatch) = addToTrack(mito{frameNum}(mitoMatch),track(trackMatch));
        else %prepare to check for fission and create a new track
            newTrack = length(track)+1;
            track(newTrack) = mito{frameNum}(mitoMatch);
        end
    end
end

%Get Dynamics of tracks
[mitoDynamics,~] = getDynamics3d(micronPerPixel,secondsPerFrame,track,zRatio);
speedStats = getDynamicStats(mitoDynamics.speed);
xyAngleStats = getDynamicStats(mitoDynamics.xyAngle);
xzAngleStats = getDynamicStats(mitoDynamics.xzAngle);
yzAngleStats = getDynamicStats(mitoDynamics.yzAngle);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Clean Up Track Gaps

checkClean = 0;
tracktotal = length(track);
cleanNum = 0;
cleanCoV = 0.1;

while ~checkClean
    cleanNum = cleanNum+1;
    lostNum = 1;
    newNum = 1;
    
    for trackNum = 1:length(track)
        
        %
        if getappdata(loading,'canceling')
            delete(loading)
            break
        end
        waitbar(trackNum/length(track),loading,sprintf('Cleaning track %d of %d.',trackNum,length(track)));
        pause(.0002)
        %

        %Make a list of lost tracks
        if (track(trackNum).frame(end) == max([track.frame]))
            continue
        end
        lostList(lostNum,:) = [trackNum,track(trackNum).frame(end),speedStats.mean(trackNum),speedStats.std(trackNum),speedStats.CoV(trackNum),xyAngleStats.mean(trackNum),xyAngleStats.std(trackNum),xyAngleStats.CoV(trackNum),xzAngleStats.mean(trackNum),xzAngleStats.std(trackNum),xzAngleStats.CoV(trackNum),yzAngleStats.mean(trackNum),yzAngleStats.std(trackNum),yzAngleStats.CoV(trackNum)];
        lostCoords(lostNum,:) = [trackNum,track(trackNum).WeightedCentroid(end-2),track(trackNum).WeightedCentroid(end-1),track(trackNum).WeightedCentroid(end)];
        lostNum = lostNum+1;

        %Make a list of new tracks
        if (track(trackNum).frame(1) == firstFrame)
            continue
        end
        newList(newNum,:) = [trackNum,track(trackNum).frame(1),speedStats.mean(trackNum),speedStats.std(trackNum),speedStats.CoV(trackNum),xyAngleStats.mean(trackNum),xyAngleStats.std(trackNum),xyAngleStats.CoV(trackNum),xzAngleStats.mean(trackNum),xzAngleStats.std(trackNum),xzAngleStats.CoV(trackNum),yzAngleStats.mean(trackNum),yzAngleStats.std(trackNum),yzAngleStats.CoV(trackNum)];
        newCoords(newNum,:) = [trackNum,track(trackNum).WeightedCentroid(1),track(trackNum).WeightedCentroid(2),track(trackNum).WeightedCentroid(3)];
        newNum = newNum+1;

    end

    if exist('newList') && exist('lostList') && sum(sum(~isnan(newList))) && sum(sum(~isnan(lostList))) 
        
        
        %Threshold based on search time
        newLostFrameMatrix = zeros(size(newList,1),size(lostList,1));

        for newTrackNum = 1:newNum-1
            newLostFrameMatrix(newTrackNum,:) = newList(newTrackNum,2)-lostList(:,2)';
        end

        newLostFrameMatrix(newLostFrameMatrix<1) = Inf;
        newLostFrameMatrix(newLostFrameMatrix>frameThresh) = Inf;

        %Threshold based on maximum travel distance
        newLostDistanceMatrix = inf(size(newList,1),size(lostList,1));

        for newTrackNum = 1:newNum-1
            for lostTrackNum = 1:lostNum-1
                if isinf(newLostFrameMatrix(newTrackNum,lostTrackNum))
                    continue
                end
                newLostDistanceMatrix(newTrackNum,lostTrackNum) = norm(newCoords(newTrackNum,2:4) - lostCoords(lostTrackNum,2:4));
            end
        end

        newLostDistanceMatrix(newLostDistanceMatrix>distThresh) = Inf;
        for newTrackNum = 1:newNum-1
            holdRow = newLostDistanceMatrix(newTrackNum,:);
            if numel(holdRow(~isinf(holdRow)))==1 && track(newList(newTrackNum,1)).NPA(1)==1
                newLostDistanceMatrix(newTrackNum,:) = newLostDistanceMatrix(newTrackNum,:)*1E-9;
            end
        end

        %Combine tracks and look at coefficient of variation
        newLostCoVMatrix = inf(size(newList,1),size(lostList,1));

        for newTrackNum = 1:newNum-1
            for lostTrackNum = 1:lostNum-1
                if isinf(newLostDistanceMatrix(newTrackNum,lostTrackNum))
                    continue
                elseif newList(newTrackNum) == lostList(lostTrackNum)
                    continue
                end
                
                tempTrack = addToTrack(track(newList(newTrackNum,1)),track(lostList(lostTrackNum,1)));
                tempDynamics = getDynamics3d(micronPerPixel,secondsPerFrame,tempTrack,zRatio);
                % maybe just replace this with 
%                 tempDynamicsSpeed = getDynamicStats(tempDynamics.speed);
                tempDynamicsxyAngle = getDynamicStats(tempDynamics.xyAngle);
                tempDynamicsxzAngle = getDynamicStats(tempDynamics.xzAngle);
                tempDynamicsyzAngle = getDynamicStats(tempDynamics.yzAngle);
                
%                 tempCoVMean = mean([tempDynamicsSpeed.CoV]);
                tempxyCoVMean = mean([tempDynamicsxyAngle.CoV]);
                tempxzCoVMean = mean([tempDynamicsxzAngle.CoV]);
                tempyzCoVMean = mean([tempDynamicsyzAngle.CoV]);

        %         tempCoVMean = tempDynamicsSpeed.CoV;

                comboCoV = 1/3*(tempxyCoVMean+tempxzCoVMean+tempyzCoVMean);
                
%                 if 1/3*(tempxyCoVMean+tempxzCoVMean+tempyzCoVMean) == 0
%                     newLostCoVMatrix(newTrackNum,lostTrackNum) = Inf;
%                     continue
%                 end

                newLostCoVMatrix(newTrackNum,lostTrackNum) = comboCoV;
            end
            holdRow = newLostCoVMatrix(newTrackNum,:);
            if numel(holdRow(~isinf(holdRow)))==1 && track(newList(newTrackNum,1)).NPA(1)==1
                newLostCoVMatrix(newTrackNum,:) = newLostCoVMatrix(newTrackNum,:)*-100;
            end
        end
        
        newLostCoVMatrix(newLostCoVMatrix == 0) = Inf;
        
%         newLostDistanceMatrix(newLostCoVMatrix>2) = Inf;
        newLostDistanceMatrix(newLostCoVMatrix>cleanCoV) = Inf;

        if numel(newLostCoVMatrix(~isinf(newLostCoVMatrix)))
            
            extra.newLostDistanceMatrix{cleanNum} = newLostDistanceMatrix;
            matchedGaps = matchpairs(newLostDistanceMatrix,distThresh);
            if numel(matchedGaps)

                %Assign tracks in descending order to avoid messing up indices
                [~,idx] = sort(matchedGaps(:,2),'descend'); % sort just the first column
                sortedGaps = matchedGaps(idx,:);   % sort the whole matrix using the sort indices
                for gapNum = 1:size(matchedGaps,1)
                    newTrackGap = newList(sortedGaps(gapNum,1));
                    lostTrackGap = lostList(sortedGaps(gapNum,2));
                    track(lostTrackGap) = addToTrack(track(newTrackGap),track(lostTrackGap));
                end

                %Delete tracks in descending order to avoid messing up indices.
                removedTracks = sort(matchedGaps(:,1),'descend');
                for trackRemove = 1:length(removedTracks)
                    track(newList(removedTracks(trackRemove))) = [];
                end
            end
        end
    end
    if length(track) == tracktotal
        checkClean = 1;
    else
        tracktotal = length(track);
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
if nargin<10
    %%Fission
    numTracks = length(track);

    %Parameters
    if nargin > 8
        stdRangeFission = stdRange;
    else
        stdRangeFission = 1;
    end

    CVpercFission = 0.01;

    %Initialize
    fissionTrack = zeros(length(track)-length(mito{firstFrame}),length(track));
    extremaFission = inf(length(track)-length(mito{firstFrame}),length(track));
    fissionMatrix = zeros(length(track)-length(mito{firstFrame}),length(track));

    %Only check track after the first frame, since no fission can occur b4 then
    newTrackIdx = ((length(mito{firstFrame})+1):numTracks);

    %Check for fission
    for trackNum = 1:length(newTrackIdx)

        %
        if getappdata(loading,'canceling')
            delete(loading)
            break
        end
        waitbar(trackNum/length(newTrackIdx),loading,sprintf('Checking fission of new track %d of %d.',trackNum,length(newTrackIdx)));
        pause(.0002)
        %

        fissionTrackFirst = zeros(frameThresh,length(track));
        extremaFissionFirst = nan(frameThresh,length(track));

        frameDifference = track(newTrackIdx(trackNum)).frame(1)-frameThresh;
        frameCheck = (max(1,frameDifference):track(newTrackIdx(trackNum)).frame(1)-1);
    %     frameCheck
        for frameNum = 1:length(frameCheck)

            for checkTrack = 1:(newTrackIdx(trackNum)-1)
                if length(track(newTrackIdx(trackNum)).Volume)>2 && sum(track(checkTrack).frame==frameCheck(frameNum)) && ~nnz(~track(newTrackIdx(trackNum)).confident) && ~nnz(~track(checkTrack).confident)
    %                 if isnan(extremaFissionFirst(frameNum,checkTrack))
    %                     continue
    %                 end

                    extremaFissionFirst(frameNum,checkTrack) = getExtremaMatrixFission3d(track(newTrackIdx(trackNum)),track(checkTrack),frameCheck(frameNum));

                    %only check track which are within search time from fission event
                    fissionTrackFirst(frameNum,checkTrack) = 1;
                end
            end
        end

        %Can only undergo fission with tracks that have a value within the
        %search time, this is our search threshold. The rows are the new tracks
        %(index corresponds to values in newTrackIdx). The cols are the
        %existing tracks. Anything above 0 is checked against the new track.
        fissionTrack(trackNum,:) = sum(fissionTrackFirst);

        extremaFission(trackNum,:) = nanmean(extremaFissionFirst);
        extremaFission(isnan(extremaFission)) = Inf;

        %Fission events must create a track that has the area of the old track
        %before fission minus the area of the old track after fission, within
        %one standard deviation and a coefficient of variation above 1.
        [fissionMatrix(trackNum,:),fissionExtra(trackNum)] = checkFissionVolume(newTrackIdx(trackNum),track,stdRangeFission,CVpercFission);

        maskFission = ((fissionTrack>0).*fissionMatrix);
    %     fissionTrack
        maskFission(maskFission == 0) = Inf;

        extremaFissionThresholded = extremaFission.*maskFission;
        extremaFissionThresholded(extremaFissionThresholded>distThresh) = Inf;


    end
    % extremaFissionThresholded
    %Assign fission
    if ~isempty(newTrackIdx)
        [fissionMin,fissionIndex] = min(extremaFissionThresholded,[],2);
        for trackNum = 1:length(newTrackIdx)
            if ~isinf(fissionMin(trackNum))
                track(newTrackIdx(trackNum)).fission(1) = fissionIndex(trackNum);
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%Fusion

    %Parameters
    if nargin > 8
        stdRangeFusion = stdRange;
    else
        stdRangeFusion = 1;
    end

    CVpercFusion = 0.01;

    lostTracks = zeros(1,numTracks);

    lostTrackNum = 0;
    for trackNum = 1:numTracks
        if track(trackNum).frame(end)~=numFrames
            lostTrackNum = lostTrackNum + 1;
            lostTracks(trackNum) = 1;
        end
    end

    %Initialize
    fusionTrack = zeros(lostTrackNum,length(track));
    extremaFusion = inf(lostTrackNum,length(track));
    fusionMatrix = zeros(lostTrackNum,length(track));

    lostTrackIdx = find(lostTracks);

    %Check for fusion
    for trackNum = 1:length(lostTrackIdx)

        %
        if getappdata(loading,'canceling')
            delete(loading)
            break
        end
        waitbar(trackNum/length(lostTrackIdx),loading,sprintf('Checking fusion of lost track %d of %d.',trackNum,length(lostTrackIdx)));
        pause(.0002)
        %

        fusionTrackFirst = zeros(frameThresh,length(track));
        extremaFusionFirst = nan(frameThresh,length(track));

        frameDifference = frameThresh+track(lostTrackIdx(trackNum)).frame(end);
        frameCheck = (track(lostTrackIdx(trackNum)).frame(end)+1:min(frameDifference,numFrames));

        for frameNum = 1:length(frameCheck)
            for checkTrack = (lostTrackIdx(trackNum)+1):length(track)%1:(lostTrackIdx(trackNum)-1)
                if length(track(lostTrackIdx(trackNum)).Volume)>2 && sum(track(checkTrack).frame==frameCheck(frameNum)) && ~nnz(~track(lostTrackIdx(trackNum)).confident) && ~nnz(~track(checkTrack).confident)
                    %Now check the extrema and assign a fission event to the 
                    %closest existing track that withstood thresholding
                    extremaFusionFirst(frameNum,checkTrack) = getExtremaMatrixFusion3d(track(lostTrackIdx(trackNum)),track(checkTrack),frameCheck(frameNum));

                    %only check track which are within search time from fission event
                    fusionTrackFirst(frameNum,checkTrack) = 1;
                end
            end
        end

        %Can only undergo fusion with tracks that have a value within the
        %search time, this is our search threshold. The rows are the lost tracks
        %(index corresponds to values in lostTrackIdx). The cols are the
        %existing tracks. Anything above 0 is checked against the lost track.
        fusionTrack(trackNum,:) = sum(fusionTrackFirst);


        extremaFusion(trackNum,:) = nanmean(extremaFusionFirst);
        extremaFusion(isnan(extremaFusion)) = Inf;

        %Fusion events must lose a track that has the area of the fused track
        %after fusion plus the area of the fused track before fission, within
        %one standard deviation and a coefficient of variance above 1.
        [fusionMatrix(trackNum,:),fusionExtra(trackNum)] = checkFusionVolume(lostTrackIdx(trackNum),track,stdRangeFusion,CVpercFusion);

        maskFusion = ((fusionTrack>0).*fusionMatrix);
        maskFusion(maskFusion == 0) = Inf;

        extremaFusionThresholded = extremaFusion.*maskFusion;
        extremaFusionThresholded(extremaFusionThresholded>distThresh) = Inf;


    end

    %Assign fission
    if ~isempty(lostTrackIdx)
        [fusionMin,fusionIndex] = min(extremaFusionThresholded,[],2);
        for trackNum = 1:length(lostTrackIdx)
            if ~isinf(fusionMin(trackNum))
                track(lostTrackIdx(trackNum)).fusion(end) = fusionIndex(trackNum);
            end
        end
    end
    
    
    extra.fissionTrack = fissionTrack;
    extra.extremaFission = extremaFission;
    extra.fissionMatrix = fissionMatrix;
    if ~isempty(newTrackIdx)
        extra.fissionExtra = fissionExtra;
        extra.fissionMin = fissionMin;
        extra.fissionIndex = fissionIndex;
        extra.newTrackIndex = newTrackIdx;
    end
    extra.fusionTrack = fusionTrack;
    extra.extremaFusion = extremaFusion;
    extra.fusionMatrix = fusionMatrix;
    if ~isempty(lostTrackIdx)
        extra.fusionExtra = fusionExtra;
        extra.fusionMin = fusionMin;
        extra.fusionIndex = fusionIndex;
        extra.lostTrackIndex = lostTrackIdx;
    end
end

for trackNum = 1:length(track)
    track(trackNum).MajorAxisLength = track(trackNum).MajorAxisLength*micronPerPixel;
    track(trackNum).MinorAxisLength = track(trackNum).MinorAxisLength*micronPerPixel;
    track(trackNum).ZAxisLength = track(trackNum).ZAxisLength*zDistanceMicrons;
    track(trackNum).SurfaceArea = track(trackNum).SurfaceArea*micronPerPixel*zDistanceMicrons;
    track(trackNum).Volume = track(trackNum).Volume*micronPerPixel^2*zDistanceMicrons;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
extra.CC = CC;
extra.L = L;
extra.mito = mito;
extra.numMito = numMito;
extra.NN = NN;
extra.NNExtrema = NNExtrema;
extra.distanceMatrix = distanceMatrix;
extra.areaMatrix = volumeMatrix;
extra.majAxMatrix = majAxMatrix;
extra.minAxMatrix = minAxMatrix;
extra.eccMatrix = solMatrix;
extra.periMatrix = surfMatrix;
extra.intMatrix = intMatrix;
extra.frameMatrix = frameMatrix;
extra.NPA = NPA;
extra.distanceMatrixThresholded = distanceMatrixThresholded;
extra.distanceMatrixSq = distanceMatrixSq;
extra.areaMatrixSq = volumeMatrixSq;
extra.majAxMatrixSq = majAxMatrixSq;
extra.minAxMatrixSq = minAxMatrixSq;
extra.eccMatrixSq = solMatrixSq;
extra.periMatrixSq = surfMatrixSq;
extra.intMatrixSq = intMatrixSq;
extra.distanceMatrixZ = distanceMatrixZ;
extra.areaMatrixZ = volumeMatrixZ;
extra.majAxMatrixZ = majAxMatrixZ;
extra.minAxMatrixZ = minAxMatrixZ;
extra.eccMatrixZ = solMatrixZ;
extra.periMatrixZ = surfMatrixZ;
extra.intMatrixZ = intMatrixZ;
extra.differenceMatrix = differenceMatrix;
extra.assignTrackMatrix = assignTrackMatrix;
extra.matchedTrack = matchedTrack;
extra.confident = confident;
delete(loading)
end

