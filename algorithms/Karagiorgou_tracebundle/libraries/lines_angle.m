function p = lines_angle(P,R,S,T)
% Returns the angle between two lines PR, ST in degrees.
p = (rad2deg(angle2Points(P,R))-rad2deg(angle2Points(S,T)));
