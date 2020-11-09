function correctedAngle = angle0to180(angle)

add180 = find(angle<0);
sub180 = find(angle>179);
angle(add180) = angle(add180)+180;
angle(sub180) = angle(sub180)-180;

% while angle < 0
%     angle = angle+180;
% end
% 
% while angle > 179
%     angle = angle-180;
% end

correctedAngle = angle;

end