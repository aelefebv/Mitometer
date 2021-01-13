function [micronPerPixel, secondsPerFrame] = start2D()

prompt2D = {'What is your pixel size, in microns per pixel? ','What is your time between frames, in seconds? '};
dlgtitle2D = 'Start 2D';
dims2D = [1 35];
definput2D = {'0.0879','1'};
answer2D = inputdlg(prompt2D,dlgtitle2D,dims2D,definput2D);
micronPerPixel = str2double(answer2D{1});
secondsPerFrame = str2double(answer2D{2});

end