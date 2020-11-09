function [track,mitoCoM,extra] = trackMitochondria(Im,ImOG,weights,micronPerPixel,secondsPerFrame,distThreshMicrons,frameThreshSeconds,stdRange,skipDyn)


numFrames = size(Im,3);

% User Parameters
distThresh = (distThreshMicrons/micronPerPixel)*secondsPerFrame;
frameThresh = round(frameThreshSeconds/secondsPerFrame);
% Other Parameters
wArea = weights(1);
wMajAx = weights(2);
wMinAx = weights(3);
wSol = weights(4);
wPeri = weights(5);
wInt = weights(6);
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
numMito = zeros(1,numFrames);
NN = cell(1,numFrames);
NNExtrema = cell(1,numFrames);
firstUsableFrame = true;
secondUsableFrame = false;
mitoCoM = struct();
distanceMatrix = cell(1,numFrames);
areaMatrix = cell(1,numFrames);
majAxMatrix = cell(1,numFrames);
minAxMatrix = cell(1,numFrames);
solMatrix = cell(1,numFrames);
periMatrix = cell(1,numFrames);
intMatrix = cell(1,numFrames);
frameMatrix = cell(1,numFrames);
NPA = cell(1,numFrames);
distanceMatrixThresholded = cell(1,numFrames);
distanceMatrixSq = cell(1,numFrames);
areaMatrixSq = cell(1,numFrames);
majAxMatrixSq = cell(1,numFrames);
minAxMatrixSq = cell(1,numFrames);
solMatrixSq = cell(1,numFrames);
periMatrixSq = cell(1,numFrames);
intMatrixSq = cell(1,numFrames);
distanceMatrixZ = cell(1,numFrames);
areaMatrixZ = cell(1,numFrames);
majAxMatrixZ = cell(1,numFrames);
minAxMatrixZ = cell(1,numFrames);
solMatrixZ = cell(1,numFrames);
periMatrixZ = cell(1,numFrames);
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

for frameNum = 1:numFrames
    
    %
    if getappdata(loading,'canceling')
        delete(loading)
        break
    end
    waitbar(frameNum/numFrames,loading,sprintf('Tracking frame %d of %d.',frameNum,numFrames));
    pause(.0002)
    %
    
    %Get stats for each mitochondrion in each frame.
    CC{frameNum} = bwconncomp(Im(:,:,frameNum));
    L{frameNum} = labelmatrix(CC{frameNum});
    mito{frameNum} = regionprops(CC{frameNum},Im(:,:,frameNum),'Area','Centroid','Solidity','MaxFeretProperties','Extrema','MeanIntensity','WeightedCentroid','MajorAxisLength','MinorAxisLength','Perimeter','PixelIdxList');
    numMito(frameNum) = length(mito{frameNum});
    
    %Keep iterating through frames until one with mitochondria appears.
    if ~sum(numMito)
        firstFrame = firstFrame+1;
        continue
    end
    
    %Get the nearest neighbor distances based on intensity weighted
    %centroids and based on extrema.
    NN{frameNum} = getNN(mito{frameNum});
    NNExtrema{frameNum} = getNNExtrema(mito{frameNum});
    
    %Correct the angles of the mitochondria and assign NN, frame, and
    %initiliaze other parameters.
    for mitoNum = 1:numMito(frameNum)
        mito{frameNum}(mitoNum).MaxFeretAngle = angle0to180(180-mito{frameNum}(mitoNum).MaxFeretAngle);
        mito{frameNum}(mitoNum).NN = NN{frameNum}(mitoNum);
        mito{frameNum}(mitoNum).NNExtrema = NNExtrema{frameNum}(mitoNum);
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
        binaryOG = true(size(ImOG(:,:,frameNum)));
        labeledOG = logical(binaryOG);
        CoMstats = regionprops(labeledOG,ImOG(:,:,frameNum),'WeightedCentroid');
        mitoCoM(frameNum).centroid = CoMstats.WeightedCentroid;
        mitoCoM(frameNum).x = CoMstats.WeightedCentroid(1);
        mitoCoM(frameNum).y = CoMstats.WeightedCentroid(2);
   
        continue
    end
        
    %Keep track of mito center of mass
    binaryOG = true(size(ImOG(:,:,frameNum)));
    labeledOG = logical(binaryOG);
    CoMstats = regionprops(labeledOG,ImOG(:,:,frameNum),'WeightedCentroid');
    mitoCoM(frameNum).centroid = CoMstats.WeightedCentroid;
    mitoCoM(frameNum).x = CoMstats.WeightedCentroid(1);
    mitoCoM(frameNum).y = CoMstats.WeightedCentroid(2);
   
    %Get the costs for the cost matrix
    distanceMatrix{frameNum} = getDistanceMatrix(mito{frameNum},track);
    areaMatrix{frameNum} = getDeltaMatrix(mito{frameNum},track,'Area');
    majAxMatrix{frameNum} = getDeltaMatrix(mito{frameNum},track,'MajorAxisLength');
    minAxMatrix{frameNum} = getDeltaMatrix(mito{frameNum},track,'MinorAxisLength');
    solMatrix{frameNum} = getDeltaMatrix(mito{frameNum},track,'Solidity');
    periMatrix{frameNum} = getDeltaMatrix(mito{frameNum},track,'Perimeter');
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
    areaMatrixSq{frameNum} = areaMatrix{frameNum}.^2;
    majAxMatrixSq{frameNum} = majAxMatrix{frameNum}.^2;
    minAxMatrixSq{frameNum} = majAxMatrix{frameNum}.^2;
    solMatrixSq{frameNum} = solMatrix{frameNum}.^2;
    periMatrixSq{frameNum} = majAxMatrix{frameNum}.^2;
    intMatrixSq{frameNum} = intMatrix{frameNum}.^2;
    
    solMatrixSq{frameNum}(isinf(solMatrixSq{frameNum})) = max(solMatrixSq{frameNum}(~isinf(solMatrixSq{frameNum})));
    solMatrixSq{frameNum}(isnan(solMatrixSq{frameNum})) = max(solMatrixSq{frameNum}(~isnan(solMatrixSq{frameNum})));

    %Z-score of the difference matrices to normalize variables
    distanceMatrixZ{frameNum} = distanceMatrixSq{frameNum};
    areaMatrixZ{frameNum} = distanceMatrixSq{frameNum};
    majAxMatrixZ{frameNum} = distanceMatrixSq{frameNum};
    minAxMatrixZ{frameNum} = distanceMatrixSq{frameNum};
    solMatrixZ{frameNum} = distanceMatrixSq{frameNum};
    periMatrixZ{frameNum} = distanceMatrixSq{frameNum};
    intMatrixZ{frameNum} = distanceMatrixSq{frameNum};

    distanceMatrixZ{frameNum}(~isinf(distanceMatrixSq{frameNum})) = zscore(distanceMatrixSq{frameNum}(~isinf(distanceMatrixSq{frameNum})),0,'all');
    areaMatrixZ{frameNum}(~isinf(distanceMatrixSq{frameNum})) = zscore(areaMatrixSq{frameNum}(~isinf(distanceMatrixSq{frameNum})),0,'all');
    majAxMatrixZ{frameNum}(~isinf(distanceMatrixSq{frameNum})) = zscore(majAxMatrixSq{frameNum}(~isinf(distanceMatrixSq{frameNum})),0,'all');
    minAxMatrixZ{frameNum}(~isinf(distanceMatrixSq{frameNum})) = zscore(minAxMatrixSq{frameNum}(~isinf(distanceMatrixSq{frameNum})),0,'all');
    solMatrixZ{frameNum}(~isinf(distanceMatrixSq{frameNum})) = zscore(solMatrixSq{frameNum}(~isinf(distanceMatrixSq{frameNum})),0,'all');
    periMatrixZ{frameNum}(~isinf(distanceMatrixSq{frameNum})) = zscore(periMatrixSq{frameNum}(~isinf(distanceMatrixSq{frameNum})),0,'all');
    intMatrixZ{frameNum}(~isinf(distanceMatrixSq{frameNum})) = zscore(intMatrixSq{frameNum}(~isinf(distanceMatrixSq{frameNum})),0,'all');

    
    %Create a matrix for assigning mito as new tracks instead
%     newTrackMatrix = Inf(numMito(frameNum),numMito(frameNum));
%     for mitoNum = 1:numMito(frameNum)
%         newTrackMatrix(mitoNum,mitoNum) = 1;
%     end
    
    %Create a cost matrix of the z-scores and append newTrackMatrix
    differenceMatrix{frameNum} = distanceMatrixZ{frameNum} + (wArea * areaMatrixZ{frameNum}) + (wMajAx * majAxMatrixZ{frameNum}) + (wMinAx * minAxMatrixZ{frameNum}) + (wSol * solMatrixZ{frameNum}) + (wPeri * periMatrixZ{frameNum}) + (wInt * intMatrixZ{frameNum});
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
[mitoDynamics,~] = getDynamics(micronPerPixel,secondsPerFrame,track);
speedStats = getDynamicStats(mitoDynamics.speed);
angleStats = getDynamicStats(mitoDynamics.travelAngle);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Clean Up Track Gaps

checkClean = 0;
tracktotal = length(track);
cleanNum = 0;
cleanCoV = 0.2;

run = 1;

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
        lostList(lostNum,:) = [trackNum,track(trackNum).frame(end),speedStats.mean(trackNum),speedStats.std(trackNum),speedStats.CoV(trackNum),angleStats.mean(trackNum),angleStats.std(trackNum),angleStats.CoV(trackNum)];
        lostCoords(lostNum,:) = [trackNum,track(trackNum).WeightedCentroid(end-1),track(trackNum).WeightedCentroid(end)];
        lostNum = lostNum+1;

        %Make a list of new tracks
        if (track(trackNum).frame(1) == firstFrame)
            continue
        end
        newList(newNum,:) = [trackNum,track(trackNum).frame(1),speedStats.mean(trackNum),speedStats.std(trackNum),speedStats.CoV(trackNum),angleStats.mean(trackNum),angleStats.std(trackNum),angleStats.CoV(trackNum)];
        newCoords(newNum,:) = [trackNum,track(trackNum).WeightedCentroid(1),track(trackNum).WeightedCentroid(2)];
        newNum = newNum+1;

    end

    if exist('newList') && exist('lostList') && sum(sum(~isnan(newList))) && sum(sum(~isnan(lostList))) 
        %Threshold based on search time
        newLostFrameMatrix = zeros(size(newList,1),size(lostList,1));

        for newTrackNum = 1:newNum-1
            newLostFrameMatrix(newTrackNum,:) = newList(newTrackNum,2)-lostList(:,2)';
        end
        newLostFrameMatrixHold = newLostFrameMatrix;
        newLostFrameMatrix(newLostFrameMatrix<1) = Inf;
        newLostFrameMatrix(newLostFrameMatrix>frameThresh) = Inf;

        %Threshold based on maximum travel distance
        newLostDistanceMatrix = inf(size(newList,1),size(lostList,1));

        for newTrackNum = 1:newNum-1
            for lostTrackNum = 1:lostNum-1
                if isinf(newLostFrameMatrix(newTrackNum,lostTrackNum))
                    continue
                end
                newLostDistanceMatrix(newTrackNum,lostTrackNum) = norm(newCoords(newTrackNum,2:3) - lostCoords(lostTrackNum,2:3));
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
                tempDynamics = getDynamics(micronPerPixel,secondsPerFrame,tempTrack);
                % maybe just replace this with 
%                 tempDynamicsSpeed = getDynamicStats(tempDynamics.speed);
                tempDynamicsAngle = getDynamicStats(tempDynamics.travelAngle);

%                 tempCoVMean = mean([tempDynamicsSpeed.CoV]);
                tempCoVMean = mean([tempDynamicsAngle.CoV]);

        %         tempCoVMean = tempDynamicsSpeed.CoV;

                if tempCoVMean == 0
                    newLostCoVMatrix(newTrackNum,lostTrackNum) = Inf;
                    continue
                end

                newLostCoVMatrix(newTrackNum,lostTrackNum) = tempCoVMean;
            end
                            
            holdRow = newLostCoVMatrix(newTrackNum,:);
            if numel(holdRow(~isinf(holdRow)))==1 && track(newList(newTrackNum,1)).NPA(1)==1
                newLostCoVMatrix(newTrackNum,:) = newLostCoVMatrix(newTrackNum,:)*-100;
            end
        end

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
    trackCount(run) = tracktotal;
    run = run+1;
    if length(track) == tracktotal
        checkClean = 1;
    else
        tracktotal = length(track);
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
if nargin<9
    %%Fission
    numTracks = length(track);

    %Parameters
    if nargin > 7
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
                if length(track(newTrackIdx(trackNum)).Area)>2 && sum(track(checkTrack).frame==frameCheck(frameNum)) && ~nnz(~track(newTrackIdx(trackNum)).confident) && ~nnz(~track(checkTrack).confident)
    %                 if isnan(extremaFissionFirst(frameNum,checkTrack))
    %                     continue
    %                 end

                    %Now check the extrema and assign a fission event to the 
                    %closest existing track that withstood thresholding
                    extremaFissionFirst(frameNum,checkTrack) = getExtremaMatrixFission(track(newTrackIdx(trackNum)),track(checkTrack),frameCheck(frameNum));

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
        %one standard deviation and a coefficient of variance above 1.
        [fissionMatrix(trackNum,:),fissionExtra(trackNum)] = checkFissionArea(newTrackIdx(trackNum),track,stdRangeFission,CVpercFission);

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
    if nargin > 7
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
                if length(track(lostTrackIdx(trackNum)).Area)>2 && sum(track(checkTrack).frame==frameCheck(frameNum)) && ~nnz(~track(lostTrackIdx(trackNum)).confident) && ~nnz(~track(checkTrack).confident)
                    %Now check the extrema and assign a fission event to the 
                    %closest existing track that withstood thresholding
                    extremaFusionFirst(frameNum,checkTrack) = getExtremaMatrixFusion(track(lostTrackIdx(trackNum)),track(checkTrack),frameCheck(frameNum));

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
        [fusionMatrix(trackNum,:),fusionExtra(trackNum)] = checkFusionArea(lostTrackIdx(trackNum),track,stdRangeFusion,CVpercFusion);

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
    track(trackNum).Perimeter = track(trackNum).Perimeter*micronPerPixel;
    track(trackNum).Area = track(trackNum).Area*micronPerPixel^2;
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
extra.areaMatrix = areaMatrix;
extra.majAxMatrix = majAxMatrix;
extra.minAxMatrix = minAxMatrix;
extra.solMatrix = solMatrix;
extra.periMatrix = periMatrix;
extra.intMatrix = intMatrix;
extra.frameMatrix = frameMatrix;
extra.NPA = NPA;
extra.distanceMatrixThresholded = distanceMatrixThresholded;
extra.distanceMatrixSq = distanceMatrixSq;
extra.areaMatrixSq = areaMatrixSq;
extra.majAxMatrixSq = majAxMatrixSq;
extra.minAxMatrixSq = minAxMatrixSq;
extra.solMatrixSq = solMatrixSq;
extra.periMatrixSq = periMatrixSq;
extra.intMatrixSq = intMatrixSq;
extra.distanceMatrixZ = distanceMatrixZ;
extra.areaMatrixZ = areaMatrixZ;
extra.majAxMatrixZ = majAxMatrixZ;
extra.minAxMatrixZ = minAxMatrixZ;
extra.solMatrixZ = solMatrixZ;
extra.periMatrixZ = periMatrixZ;
extra.intMatrixZ = intMatrixZ;
extra.differenceMatrix = differenceMatrix;
extra.assignTrackMatrix = assignTrackMatrix;
extra.matchedTrack = matchedTrack;
extra.confident = confident;
delete(loading)
end

