function updates = checkUpdates()

thisVersion = textread('version.txt');

[newVersion,status] = urlread('https://raw.githubusercontent.com/aelefebv/Mitometer/main/Source%20Code/version.txt');

if status~=0 && thisVersion<str2double(newVersion)
    updates = 1;
else
    updates = 0;
end

end