% Shortest-path map construction evaluation 1.0
% Copyright 2013 Sophia Karagiorgou and Dieter Pfoser
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
% 
% ------------------------------------------------------------------------
% 
% This software is based on the following article(s). Please cite this
% article when using this code as part of a research publication:
% 
% S. Karagiorgou and D. Pfoser. On Vehicle Tracking Data-Based Road Network Generation. 
% In Proc. 20th ACM SIGSPATIAL GIS Conference, pp. 89-98, 2012. 
% 
% S. Karagiorgou and D. Pfoser. Segmentation-Based Road Network Construction. 
% In Proc. 21th ACM SIGSPATIAL GIS Conference, pp. 470-473, 2013.   
% 
% ------------------------------------------------------------------------
% 
% Author: Sophia Karagiorgou (karagior@gmail.com)
% Filename: shortest_paths_comparison.m

currentFolder = pwd;
cd(currentFolder);
fprintf('Starting at %s\n', timestamp);

vertices = load('vertices_original_osm.txt');
edges = load('edges_original_osm.txt');


hold on;
line = zeros(2);

numofshortestpaths=500;
dy=50;
bbox=[];
rbbox=[];
rids=[];
new_edges=[];
new_vertices=[];
fprintf('Loading ground truth road network\n');
for k=1:length(edges)
    rn1=find(vertices(:,1)==edges(k,2));
    rn2=find(vertices(:,1)==edges(k,3));
    
    ndegrees=line_angle([vertices(rn1,2),vertices(rn1,3)],[vertices(rn2,2),vertices(rn2,3)]);
    x1=vertices(rn1,2)+(dy)*cosd(ndegrees-90);
    y1=vertices(rn1,3)+(dy)*sind(ndegrees-90);
    x2=vertices(rn1,2)+(dy)*cosd(ndegrees+90);
    y2=vertices(rn1,3)+(dy)*sind(ndegrees+90);
    
    x3=vertices(rn2,2)+(dy)*cosd(ndegrees-90);
    y3=vertices(rn2,3)+(dy)*sind(ndegrees-90);
    x4=vertices(rn2,2)+(dy)*cosd(ndegrees+90);
    y4=vertices(rn2,3)+(dy)*sind(ndegrees+90);
    
    line(1,1)= vertices(rn1,2);
    line(1,2)= vertices(rn1,3);
    line(2,1)= vertices(rn2,2);
    line(2,2)= vertices(rn2,3);
    
    node1=[x1 y1;x3 y3];
    node2=[x4 y4;x2 y2];
    bbox=vertcat(bbox,[node1(:,1) node1(:,2)]);
    bbox=vertcat(bbox,[node2(:,1) node2(:,2)]);
    bbox=vertcat(bbox,[x1 y1]);
    bbox=vertcat(bbox,[NaN NaN]);
    new_edges=vertcat(new_edges,[edges(k,1) edges(k,2) edges(k,3)]);
    
    new_vertices=vertcat(new_vertices,[vertices(rn1,1) vertices(rn1,2) vertices(rn1,3)]);
    new_vertices=vertcat(new_vertices,[vertices(rn2,1) vertices(rn2,2) vertices(rn2,3)]);
end


nedges=load('tracebundle_edges.txt');
nvertices=load('tracebundle_vertices.txt');
fprintf('Loading generated road network\n');
fprintf('Extracting correspondent road network links\n');
hold on;
trids=[];
nedges=unique(nedges,'rows');
for k=1:length(nedges(:,1))
    c1=find(nvertices(:,1)==nedges(k,2));
    c2=find(nvertices(:,1)==nedges(k,3));
    
    xy=[];
    xy=[nvertices(c1,2) nvertices(c1,3)];
    xy=vertcat(xy,[nvertices(c2,2) nvertices(c2,3)]);
    line=createLine([nvertices(c1,2) nvertices(c1,3)], [nvertices(c2,2) nvertices(c2,3)]);
    rc=isnan(bbox(:,1));
    irc=find(rc==1);
    
    for m=1:length(irc)
        tbbox=[];
        if m==1
            rc1=1;
        else
            rc1=irc(m-1)+1;
        end
        
        rc2=irc(m)-1;
        tbbox=bbox(rc1:rc2,:);
        [in on]=inpolygon(xy(:,1),xy(:,2),tbbox(:,1),tbbox(:,2));
        poly=[tbbox(1,1) tbbox(1,2);tbbox(2,1) tbbox(2,2);tbbox(3,1) tbbox(3,2);tbbox(4,1) tbbox(4,2)];
        
        [xi, yi, ii] = polyxpoly(xy(:,1), xy(:,2), poly(:,1), poly(:,2));
        ri=find(in==1);
        ro=find(on==1);
        
        if ~isempty(xi)&&~isempty(yi)||(~isempty(ri)||~isempty(ro))
            if ~isempty(ri)
                for n=1:length(ri)
                    trids=vertcat(trids,[nedges(k,1) m ri(n)]);
                end
            end
            if ~isempty(ro)
                for n=1:length(ro)
                    trids=vertcat(trids,[nedges(k,1) m ro(n)]);
                end
            end
            if ~isempty(xi)
                for n=1:length(xi)
                    trids=vertcat(trids,[nedges(k,1) m ii(n,1)]);
                end
            end
        end
    end
end
dtrids=unique(trids(:,2));

redges=[];
rnodes=[];
line = zeros(2);
cnt=1;

for k=1:length(dtrids(:,1))
    redges=vertcat(redges,[cnt edges(dtrids(k,1),2) edges(dtrids(k,1),3)]);
    rn1=find(vertices(:,1)==edges(dtrids(k,1),2));
    rn2=find(vertices(:,1)==edges(dtrids(k,1),3));
    rnodes=vertcat(rnodes, [vertices(rn1,1) vertices(rn1,2) vertices(rn1,3)]);
    rnodes=vertcat(rnodes, [vertices(rn2,1) vertices(rn2,2) vertices(rn2,3)]);
    cnt=cnt+1;
end

fprintf('Creating adjacency matrices\n');
redges=unique(redges,'rows');
rnodes=unique(rnodes,'rows');
matrix_r=inf(length(rnodes(:,1)), length(rnodes(:,1)));
for i=1:length(rnodes(:,1))
    for j=1:length(rnodes(:,1))
        if i==j
            matrix_r(i,j)=1;
        end
    end
end

for k=1:length(redges(:,1))
    ri=find(rnodes(:,1)==redges(k,2));
    rj=find(rnodes(:,1)==redges(k,3));
    if ~isempty(ri)&&~isempty(rj)
        cst=sqrt((rnodes(ri,2)-rnodes(rj,2))^2+(rnodes(ri,3)-rnodes(rj,3))^2);
        matrix_r(ri, rj) = cst;
        matrix_r(rj, ri) = cst;
    end
end


minx=min(nvertices(:,2))-500;
maxx=max(nvertices(:,2))+500;
miny=min(nvertices(:,3))-500;
maxy=max(nvertices(:,3))+500;

rc=find(nvertices(:,2)>=minx&nvertices(:,2)<=maxx&nvertices(:,3)>=miny&nvertices(:,3)<=maxy);
cvertices_tmp=nvertices(rc,:);
pedges=[];
for g=1:length(cvertices_tmp(:,1))
    rs=find(nedges(:,2)==cvertices_tmp(g,1));
    if ~isempty(rs)
        pedges=vertcat(pedges,nedges(rs,1));
    end
    rs=find(nedges(:,3)==cvertices_tmp(g,1));
    if ~isempty(rs)
        pedges=vertcat(pedges,nedges(rs,1));
    end
end

pedges=unique(pedges,'rows');
nodes=[];

for i=1:length(pedges(:,1))
    rc=find(nedges(:,1)==pedges(i,1));
    if ~isempty(rc)
        rc1=find(nvertices(:,1)==nedges(rc,2));
        rc2=find(nvertices(:,1)==nedges(rc,3));
        nodes=vertcat(nodes,[nvertices(rc1,1) nvertices(rc1,2) nvertices(rc1,3)]);
        nodes=vertcat(nodes,[nvertices(rc2,1) nvertices(rc2,2) nvertices(rc2,3)]);
    end
end

nodes=unique(nodes,'rows');
matrix_p=inf(length(nodes(:,1)), length(nodes(:,1)));

for i=1:length(nodes(:,1))
    for j=1:length(nodes(:,1))
        if i==j
            matrix_p(i,j)=1;
        end
    end
end

for k=1:length(pedges(:,1))
    sres=find(nedges(:,1)==pedges(k,1));
    c1=find(nvertices(:,1)==nedges(sres,2));
    c2=find(nvertices(:,1)==nedges(sres,3));
    
    ri=find(nodes(:,1)==nvertices(c1,1));
    rj=find(nodes(:,1)==nvertices(c2,1));
    cst=sqrt((nodes(ri,2)-nodes(rj,2))^2+(nodes(ri,3)-nodes(rj,3))^2);
    if ~isempty(ri)&&~isempty(rj)
        matrix_p(ri, rj) = cst;
        matrix_p(rj, ri) = cst;
    end
end

activeNodes = [];
for i = 1:length(nodes(:,1))
    farthestPreviousHop(i) = i;
    farthestNextHop(i) = i;
end


fprintf('Generating uniformly distributed shortest paths\n');
path_p=[];
dpath_p=[];
dpath_p=inf(length(nodes(:,1)), length(nodes(:,1)));

for k=1:numofshortestpaths    
    rr=randi(length(nodes(:,1)),1,1);
    circle = [[nodes(rr,2) nodes(rr,3)] 5000];
    x1=nodes(rr,2)+5000;
    y1=nodes(rr,3)+5000;
    x2=nodes(rr,2)-5000;
    y2=nodes(rr,3)+5000;
    x3=nodes(rr,2)-5000;
    y3=nodes(rr,3)-5000;
    x4=nodes(rr,2)+5000;
    y4=nodes(rr,3)-5000;
    node=[x1 y1;x2 y2;x3 y3;x4 y4;x1 y1];
    rx=find(nodes(:,2)~=nodes(rr,2)&nodes(:,3)~=nodes(rr,3));
    
    in=inpolygon(nodes(rx,2), nodes(rx,3),node(:,1), node(:,2));
    rin=find(in~=0);
    
    Mdist=[nodes(rr,2) nodes(rr,3)];
    Mdist=vertcat(Mdist,[nodes(rx(rin),2) nodes(rx(rin),3)]);
    Y = pdist(Mdist,'euclid');
    SQ=squareform(Y);
    
    maxd=max(SQ(:,1));
    r=find(SQ(:,1)==maxd);
    
    rc1=find(nodes(:,1)==nodes(rr,1));
    rrc2=find(nodes(:,2)==Mdist(r,1)&nodes(:,3)==Mdist(r,2));
    rc2=find(nodes(:,1)==nodes(rrc2,1));
    [path, totalCost, farthestPreviousHop, farthestNextHop] = dijkstra_1(length(nodes(:,1)), matrix_p, rc1, rc2, farthestPreviousHop, farthestNextHop);
    d=0;
    if length(path) ~= 0
        d=0;
        for i = 1:(length(path)-1)
            d=d+sqrt((nodes(path(i),2)-nodes(path(i+1),2))^2+(nodes(path(i),3)-nodes(path(i+1),3))^2);
        end
    end
    if ~isempty(k)&& ~isempty(nodes(rr,1)) && ~isempty(nodes(rrc2,1)) && ~isempty(maxd) && ~isempty(d) && ~isempty(totalCost)
        if ~isinf(totalCost)
            path_p=vertcat(path_p,[k nodes(rr,1) nodes(rrc2,1) maxd d totalCost 0 0 0 0 0 0]);
            for i=1:length(path)
                dpath_p(k,i)=path(i);
            end
        end
    end
end


activeNodes = [];
for i = 1:length(rnodes(:,1))    
    farthestPreviousHop(i) = i;
    farthestNextHop(i) = i;
end

path_r=[];
dpath_r=[];
dpath_r=inf(length(rnodes(:,1)), length(rnodes(:,1)));
rs=~isinf(dpath_p(:,1));
dpath_p=dpath_p(rs,:);
fprintf('Calculating corresponding shortest paths for the ground truth road network\n');
for k=1:length(path_p(:,1))
    rr1=find(nodes(:,1)==path_p(k,2));
    rr2=find(nodes(:,1)==path_p(k,3));
    
    x1=nodes(rr1,2)+100;
    y1=nodes(rr1,3)+100;
    x2=nodes(rr1,2)-100;
    y2=nodes(rr1,3)+100;
    x3=nodes(rr1,2)-100;
    y3=nodes(rr1,3)-100;
    x4=nodes(rr1,2)+100;
    y4=nodes(rr1,3)-100;
    node=[x1 y1;x2 y2;x3 y3;x4 y4;x1 y1];
    
    if ~isnan(dpath_p(k,2))
        ndegrees=line_angle([nodes(rr1,2),nodes(rr1,3)],[nodes(dpath_p(k,2),2),nodes(dpath_p(k,2),3)]);
        rpd=find(dpath_p(k,:)~=inf);
        pdegrees=line_angle([nodes(dpath_p(k,length(rpd)-1),2),nodes(dpath_p(k,length(rpd)-1),3)],[nodes(rr2,2),nodes(rr2,3)]);
        
        in=inpolygon(rnodes(:,2), rnodes(:,3),node(:,1), node(:,2));
        rin1=find(in~=0);
        mindist=1000;
        minvdist=1000;
        idx1=0;
        for h=1:length(rin1)
            adf=find(matrix_r(rin1(h),:)~=inf);
            for n=1:length(adf)
                if adf(n)~=rin1(h)
                    ndegrees1=line_angle([rnodes(rin1(h),2),rnodes(rin1(h),3)],[rnodes(adf(n),2),rnodes(adf(n),3)]);                    
                    dp=distancePoints([nodes(rr1,2),nodes(rr1,3)], [rnodes(rin1(h),2),rnodes(rin1(h),3)],2);
                    ledge=[nodes(rr1,2),nodes(rr1,3), nodes(dpath_p(k,2),2),nodes(dpath_p(k,2),3)];                    
                    if dp<mindist&&(abs(ndegrees-ndegrees1)<50||abs(ndegrees-ndegrees1)>130&&abs(ndegrees-ndegrees1)<230)
                        mindist=dp;
                        idx1=h;
                    end
                end
            end
        end
        
        x1=nodes(rr2,2)+100;
        y1=nodes(rr2,3)+100;
        x2=nodes(rr2,2)-100;
        y2=nodes(rr2,3)+100;
        x3=nodes(rr2,2)-100;
        y3=nodes(rr2,3)-100;
        x4=nodes(rr2,2)+100;
        y4=nodes(rr2,3)-100;
        node=[x1 y1;x2 y2;x3 y3;x4 y4;x1 y1];
        
        in=inpolygon(rnodes(:,2), rnodes(:,3),node(:,1), node(:,2));
        rin2=find(in~=0);
        mindist=1000;
        minvdist=1000;
        idx2=0;
        
        for h=1:length(rin2)
            adf=find(matrix_r(rin2(h),:)~=inf);
            for n=1:length(adf)
                if adf(n)~=rin2(h)
                    ndegrees2=line_angle([rnodes(adf(n),2),rnodes(adf(n),3)],[rnodes(rin2(h),2),rnodes(rin2(h),3)]);
                    
                    dp=distancePoints([nodes(rr2,2),nodes(rr2,3)], [rnodes(rin2(h),2),rnodes(rin2(h),3)],2);
                    ledge=[nodes(dpath_p(k,length(rpd)-1),2),nodes(dpath_p(k,length(rpd)-1),3),nodes(rr2,2),nodes(rr2,3)];
                    if dp<mindist&&abs(pdegrees-ndegrees2)<60
                        mindist=dp;                        
                        idx2=h;
                    end
                end
            end
        end
        rin11=rin1(idx1);
        rin22=rin2(idx2);
        maxd=sqrt((rnodes(rin11,2)-rnodes(rin22,2))^2+(rnodes(rin11,3)-rnodes(rin22,3))^2);
        [path, totalCost, farthestPreviousHop, farthestNextHop] = dijkstra_1(length(rnodes(:,1)), matrix_r, rin11, rin22, farthestPreviousHop, farthestNextHop);
        d=0;
        if length(path) ~= 0
            d=0;
            for i = 1:(length(path)-1)
                d=d+sqrt((rnodes(path(i),2)-rnodes(path(i+1),2))^2+(rnodes(path(i),3)-rnodes(path(i+1),3))^2);                
            end
        end
        if ~isempty(k)&& ~isempty(rnodes(rin11,1)) && ~isempty(rnodes(rin22,1)) && ~isempty(maxd) && ~isempty(d) && ~isempty(totalCost)
            path_r=vertcat(path_r,[k rnodes(rin11,1) rnodes(rin22,1) maxd d totalCost 0 0 0 0 0 0]);
            for i=1:length(path)
                dpath_r(k,i)=path(i);
            end
        end
    end
end

path_p(:,7)=0;
path_r(:,7)=0;
path_p(:,8)=0;
path_r(:,8)=0;
path_p(:,9)=0;
path_r(:,9)=0;
path_p(:,10)=0;
path_r(:,10)=0;
fprintf('Calculating distance measures\n');
for k=1:length(path_p(:,1))
    ra=find(~isinf(dpath_p(k,:)));
    rb=find(~isinf(dpath_r(k,:)));
    A=[nodes(dpath_p(k,ra),2), nodes(dpath_p(k,ra),3)];
    B=[rnodes(dpath_r(k,rb),2), rnodes(dpath_r(k,rb),3)];
    if ~isempty(ra)&&~isempty(rb)
        path_p(k,7)=HausdorffDist(A,B);
        path_p(k,8)=ModHausdorffDist(A,B);
        path_p(k,9)=DiscreteFrechetDist(A,B);        
        path_p(k,10)=ModAvgDistance(A,B);       
        
        path_r(k,7)=HausdorffDist(A,B);
        path_r(k,8)=ModHausdorffDist(A,B);
        path_r(k,9)=DiscreteFrechetDist(A,B);        
        path_r(k,10)=ModAvgDistance(A,B);
    else
        path_p(k,7)=inf;
        path_p(k,8)=inf;
        path_p(k,9)=inf;        
        path_p(k,10)=inf;
        
        
        path_r(k,7)=inf;
        path_r(k,8)=inf;
        path_r(k,9)=inf;        
        path_r(k,10)=inf;        
    end
end


fprintf('Writing calculations into file\n');
wf=['shortest_path_distance_measures.txt'];
fileID = fopen(wf,'w');
fprintf(fileID,'SP_id,DirectDist,NetworkDist,HausdorffDist,ModHausdorffDist,DiscreteFrechetDist,AvgVerticalDistance');
for k=1:length(path_p(:,1))
    fprintf(fileID,'%d,%f,%f,%f,%f,%f,%f,%f\n',path_p(k,1),path_p(k,4),path_p(k,5),path_p(k,7),path_p(k,8),path_p(k,9),path_p(k,10));
end
fclose(fileID);

fprintf('Finishing at %s\n', timestamp);