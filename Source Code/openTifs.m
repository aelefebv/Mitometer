function [importedImage,numFrames,fileName,path] = openTifs(numZ,fileName,path)

if nargin <2
    [fileName,path] = uigetfile('*.tif','Select one or more thresholded files','MultiSelect','on');
end



if isequal(fileName,0)
    disp('User canceled');
    return
elseif ~iscell(fileName)
    fileName = {fileName};
% else
%     showNumberOfFiles = sprintf('Number of files selected: %d\nDo you want to continue processing,\nor Cancel to abort processing?',size(fileName,2));
%     buttonContinue = questdlg(showNumberOfFiles, 'Continue', 'Continue', 'Cancel', 'Continue');
%     if ~strcmpi(buttonContinue, 'Continue')
%         disp('User canceled');
%         return; 
%     end
end

numFrames = cell(1,length(fileName));
importedImage = cell(1,length(fileName));

for i = 1:length(fileName)
    %Get image info
    imInfo = imfinfo(strcat(path,fileName{i}));

    imCols = imInfo(1).Width;
    imRows = imInfo(1).Height;

    numFrames{i} = length(imInfo);

    %Read in each tif frame and create a stack
    tifIm = Tiff(strcat(path,fileName{i}),'r');
    
    if nargin>0 %if there is also a z component
        importedImage{i} = zeros(imRows,imCols,numFrames{i}/numZ,numZ,'uint8');
        for zNum = 1:numZ
            for frameNum = 1:numFrames{i}/numZ
                tifIm.setDirectory((zNum)+(frameNum*numZ-numZ));
                if imInfo(1).BitDepth ~= 8
                    importedImage{i}(:,:,frameNum,zNum) = uint8(255 * mat2gray(tifIm.read()));
                else
                    importedImage{i}(:,:,frameNum,zNum) = tifIm.read();
                end
            end
        end
        tifIm.close();
    else
        importedImage{i} = zeros(imRows,imCols,numFrames{i},'uint8');

        for frameNum = 1:numFrames{i}
            tifIm.setDirectory(frameNum);
            if imInfo(1).BitDepth ~= 8
                importedImage{i}(:,:,frameNum) = uint8(255 * mat2gray(tifIm.read()));
            else
                importedImage{i}(:,:,frameNum) = tifIm.read();
            end
        end
        tifIm.close();
    end
end

end

