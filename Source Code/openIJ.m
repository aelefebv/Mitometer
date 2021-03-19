function openIJ()

if ismac
    path = uigetdir("~","Select your Fiji app's folder");
    addpath(strcat(path,'/Fiji.app/scripts')) % Update for your ImageJ installation as appropriate
elseif ispc
    path = uigetdir("~","Select your Fiji app's scripts folder");
    addpath(path);
else
    disp('Platform not supported')
end
ImageJ;

end