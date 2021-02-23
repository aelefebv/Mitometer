function updates = checkUpdates()

thisVersion = textread('version.txt');

[newVersion,status] = urlread('https://raw.githubusercontent.com/aelefebv/Mitometer/Revisions/version.txt');

if status~=0 && thisVersion<str2double(newVersion)
    msgbox('There is a new version available.', 'Note')
    updates = 1;
else
    updates = 0;
end

end