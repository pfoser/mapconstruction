function cm = DiscreteFrechetDist(P,Q)
% Calculates the discrete Frechet distance between curves P and Q
%
% cm = DiscreteFrechetDist(P,Q)
%
% P and Q are two sets of points that define polygonal curves with rows of
% vertices (data points) and columns of dimensionality. The points along
% the curves are taken to be in the order as they appear in P and Q.
%
% Returned in cm is the discrete Frechet distance, aka the coupling
% measure, which is zero when P equals Q and grows positively as the curves
% become more dissimilar.
%
% Explanation:
% The Frechet distance is a measure of similarity between to curves, P and
% Q. It is defined as the minimum cord-length sufficient to join a point
% traveling forward along P and one traveling forward along Q, although the
% rate of travel for either point may not necessarily be uniform.
%
% The Frechet distance, FD, is not in general computable for any given
% continuous P and Q. However, the discrete Frechet Distance, also called
% the coupling measure, cm, is a metric that acts on the endpoints of
% curves represented as polygonal chains. The magnitude of the coupling
% measure is bounded by FD plus the length of the longest segment in either
% P or Q,  and approaches FD in the limit of sampling P and Q.
%
% This function implements the algorithm to calculate discrete Frechet
% distance outlined in:
% T. Eiter and H. Mannila. Computing discrete Frechet distance. Technical
% Report 94/64, Christian Doppler Laboratory, Vienna University of
% Technology, 1994.
% 
%
%
% EXAMPLE:
% t = 0:pi/8:2*pi;
% y = linspace(1,3,6);
% P = [(2:7)' y']+0.3.*randn(6,2);
% Q = [t' sin(t')]+2+0.3.*randn(length(t),2);
% cm = DiscreteFrechetDist(P,Q);
% figure
% plot(Q(:,1),Q(:,2),'o-r','linewidth',2,'markerfacecolor','r')
% hold on
% plot(P(:,1),P(:,2),'o-b','linewidth',2,'markerfacecolor','b')
% title(['Discrete Frechet Distance of curves P and Q: ' num2str(cm)])
% legend('Q','P','location','best')
% line([2 cm+2],[1 1],'color','k')
% text(2,0.9,'dFD length')
% 
% %%% ZCD June 2011 %%%
%

% keep storage variable in memory
persistent CA

% size of the data curves
sP = size(P);
sQ = size(Q);

% check validity of inputs
if sP(2)~=sQ(2)
    error('Curves P and Q must be of the same dimension')
elseif sP(1)==0
    cm = 0;
    return;
end

% initialize ca to a matrix of -1s
CA = ones(sP(1),sQ(1)).*-1;

% distance function
dfcn = @(u,v) sqrt(sum( (u-v).^2 ));

% final coupling measure value
cm = c(sP(1),sQ(1),P,Q,dfcn);




function CAij = c(i,j,P,Q,dfcn)
    % coupling search function
    if CA(i,j)>-1
        % don't update CA in this case
        CAij = CA(i,j);
    elseif i==1 && j==1
        CA(i,j) = dfcn(P(1,:),Q(1,:));     % update the CA perminent
        CAij = CA(i,j);                    % set the current relevant value
    elseif i>1 && j==1
        CA(i,j) = max( c(i-1,1,P,Q,dfcn), dfcn(P(i,:),Q(1,:)) );
        CAij = CA(i,j);
    elseif i==1 && j>1
        CA(i,j) = max( c(1,j-1,P,Q,dfcn), dfcn(P(1,:),Q(j,:)) );
        CAij = CA(i,j);
    elseif i>1 && j>1
        CA(i,j) = max( min([c(i-1,j,P,Q,dfcn), ...
            c(i-1,j-1,P,Q,dfcn), c(i,j-1,P,Q,dfcn)]),...
            dfcn(P(i,:),Q(j,:)) );
        CAij = CA(i,j);
    else
        CA(i,j) = inf;
    end
end     % end function, c


end     % end main function

