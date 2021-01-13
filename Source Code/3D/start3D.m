function [micronPerPixel, secondsPerFrame, zPlane, zDistance] = start3D()

prompt3D = {'What is your pixel size, in microns per pixel? ','What is your time between frames, in seconds? ', 'How many z-planes in the stack? ', 'What is the axial distance, in microns? '};
dlgtitle3D = 'Start 2D';
dims3D = [1 35];
definput3D = {'0.105','10', '21', '0.45'};
answer3D = inputdlg(prompt3D,dlgtitle3D,dims3D,definput3D);
micronPerPixel = str2double(answer3D{1});
secondsPerFrame = str2double(answer3D{2});
zPlane = str2double(answer3D{3});
zDistance = str2double(answer3D{4});

end