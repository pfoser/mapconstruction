function a = line_angle(varargin)
% Computes the angle of the line segments with respect to 
% the four quadrant coordinate system.

% process input arguments
if length(varargin)==2
    p1 = varargin{1};
    p2 = varargin{2};
elseif length(varargin)==1
    var = varargin{1};
    p1 = var(1,:);
    p2 = var(2,:);
end    

% ensure data have same size
if size(p1, 1)==1
    p1 = p1(ones(size(p2,1), 1), :);
elseif size(p2, 1)==1
    p2 = p2(ones(size(p1,1), 1), :);
elseif size(p1, 1)~=size(p2, 1)
    error('angle2Points: wrong size for inputs');
end

% angle of line (P2 P1)
dp = p2-p1;
theta = mod(atan2(dp(:,2), dp(:,1)) + 2*pi, 2*pi);

a=theta*180/pi;