function dist = gc_distance(coors,coorsi)

% a = sin((lat2-lat1)/2).^2 + cos(lat1) .* cos(lat2) .* sin((lon2-lon1)/2).^2;
% 
% % Ensure that a falls in the closed interval [0 1].
% a(a < 0) = 0;
% a(a > 1) = 1;
% 
% rng = r * 2 * atan2(sqrt(a),sqrt(1 - a));
coors = coors/180*pi;
coorsi = coorsi/180*pi;

a = sin((coorsi(:,2)-coors(:,2))/2).^2 + cos(coors(:,2)) .* cos(coorsi(:,2)) .* sin((coorsi(:,1)-coors(:,1))/2).^2;

% Ensure that a falls in the closed interval [0 1].
a(a < 0) = 0;
a(a > 1) = 1;

dist = (1 * 2 * atan2(sqrt(a),sqrt(1 - a)))/pi*180;
end