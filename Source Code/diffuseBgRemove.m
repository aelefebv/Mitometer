function [ImBgRemoved,ImMinMed,ImMedFiltered] = diffuseBgRemove(Im,minCircleFiltSize,maxCircleFiltSize)

ImMinMed = zeros(size(Im),class(Im));
ImBgRemoved = zeros(size(Im),class(Im));

minCircleFilt = max(2,minCircleFiltSize);


loading = waitbar(0,'Please wait...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(loading,'canceling',0);
pause(.0002)

for frameNum = 1:size(Im,3)
    
    %We add the original image into the stack as well.
    ImMedFiltered(:,:,1) = Im(:,:,frameNum);
    circleFilt = cell(1,round(maxCircleFiltSize) - round(minCircleFilt));

    circleFiltNum = 1;
    
    %We iterate through this process to create multiple median disk filters
    %in radius steps of 1 pixel from the minimum radius to the maximum
    %radius to create a stack of median disk filters.
    for circleSize = round(minCircleFilt):round(maxCircleFiltSize)
        %
        if getappdata(loading,'canceling')
            delete(loading)
            break
        end
        waitbar(circleSize/round(maxCircleFiltSize),loading,sprintf('Median filter %d of %d for frame %d of %d.',circleSize,round(maxCircleFiltSize),frameNum,size(Im,3)));
        pause(.0002)
        %
        circleFiltNum = circleFiltNum + 1;
        
        %We begin by creating a two-dimensional disk filters using MATLAB’s
        %fspecial function. We use this disk filter in MATLAB’s ordfilt2
        %function using the median value as the order, thereby creating a
        %median disk filter.
       
        circleFilt{circleFiltNum-1} = fspecial('disk',circleSize)>0;
        ImMedFiltered(:,:,circleFiltNum) = ordfilt2(ImMedFiltered(:,:,1),round(0.5*numel(find(circleFilt{circleFiltNum-1}))),circleFilt{circleFiltNum-1},'symmetric');
    end

    %We then find the stack-wide minimum value for each pixel to create the
    %diffuse background image and subtract this from the original image.
    ImMinMed(:,:,frameNum) = min(ImMedFiltered,[],3);

    ImBgRemoved(:,:,frameNum) = ImMedFiltered(:,:,1)-ImMinMed(:,:,frameNum);

    %We repeat this process for each frame and each z-slice in the case of
    %3D to remove the diffuse background image from the entire fluorescence
    %image stack.
end
delete(loading)

end