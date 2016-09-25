% TraceBundle map construction 1.0
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
% Filename: intersection_nodes_extraction.m

clear;
currentFolder = pwd;
cd(currentFolder);
fprintf('Starting at %s\n', timestamp);

filetype=1;
total=0;

rz=[];
angles=[];
pts_ln=[];
turn_pts=[];
turn_seeds=[];
centre_pts=[];
clusterseeds=[];
pts_to_inter=[];
inter_to_cluster=[];

dist=40; %clustering distance range
angle_tol = 12.5; % angle tolerance after which a point is considered an interstion or turn

if filetype==0
    [fnames,path]=uigetfile('*.dat','...select trajectory files!','MultiSelect','on');
    lgth=length(fnames);
    n=ceil(sqrt(lgth));
    m=floor(lgth/n);
else
    [fnames,path]=uigetfile('*.txt','...select trajectory files!','MultiSelect','on');
    lgth=1;
    n=ceil(sqrt(lgth));
    m=floor(lgth/n);
end

% ---------------------------
% loop over trajectory files
% ---------------------------
ss=0;
dd=0;

rectr=[];
my_corner_pts=[];
mcnt=1;


fprintf('Looping over trajectories\n');


for g=1:length(fnames(1,:))
    mat=load([path,fnames{1,g}]);
    sm=size(mat);
    len=sm(1);
    
    x=mat(:,1);
    y=mat(:,2);
    t=mat(:,3);
    
    subtr=1;
    isubtr=1;
    
    %produce delta times for pt2-pt1, etc.
    dt=t(2:len) - t(1:len-1); % times for pt2 and on
    
    %distance between points pt2-pt1, etc.
    ds=sqrt((x(2:len)-x(1:len-1)).^2+(y(2:len)-y(1:len-1)).^2);
    dv=(ds./dt)/1000*3600;
    
    dx=(x(2:len)-x(1:len-1)); % diff p2-p1, p3-p2, etc.
    dy=(y(2:len)-y(1:len-1));
    
    ta=dy./dx; %tan a
    
    rx=dx./dt;
    ry=dy./dt;
    
    xy=[];
    imem=0;
    jmem=0;
    
    for j=2:len-1 % starts for second point, but index for ta is len-1, hence the -1
        
        if ds(j-1)>=1000 || (dt(j-1)>=45 && ds(j-1)>=1000) || (imem == g && abs(j-jmem)>1)
            subtr=subtr+1;
            isubtr=1;
        end
        
        if j==len-1
            p1=[x(j-1),y(j-1)];
            p2=[x(j),y(j)];
            a1=rad2deg(angle2Points(p1,p2));
            
            p1=[x(j),y(j)];
            p2=[x(j+1),y(j+1)];
            a2=rad2deg(angle2Points(p1,p2));
            
            pts_ln=vertcat(pts_ln,[x(j) y(j) g j t(j) a1 a2 subtr isubtr]);
            imem=g;
            jmem=j;
            
            if ds(j)>=1000 || (dt(j)>=45 && ds(j)>=1000)
                subtr=subtr+1;
                isubtr=1;
                pts_ln=vertcat(pts_ln,[x(j+1) y(j+1) g j+1 t(j+1) a2 -1 subtr isubtr]);
                imem=g;
                jmem=j+1;
            else
                isubtr=isubtr+1;
                pts_ln=vertcat(pts_ln,[x(j+1) y(j+1) g j+1 t(j+1) a2 -1 subtr isubtr]);
                imem=g;
                jmem=j+1;
            end
            
            rectr=vertcat(rectr,[x(j) y(j) g j t(j) x(j+1) y(j+1) g j+1 t(j+1) subtr isubtr]);
            
        elseif j==2
            
            p1=[x(j-1),y(j-1)];
            p2=[x(j),y(j)];
            a1=rad2deg(angle2Points(p1,p2));
            
            p1=[x(j),y(j)];
            p2=[x(j+1),y(j+1)];
            a2=rad2deg(angle2Points(p1,p2));
            
            if ds(j-1)>=1000 || (dt(j-1)>=45 && ds(j-1)>=1000)
                pts_ln=vertcat(pts_ln,[x(j-1) y(j-1) g j-1 t(j-1) -1 a1 subtr isubtr]);
                subtr=subtr+1;
                isubtr=1;
                imem=g;
                jmem=j-1;
            else
                pts_ln=vertcat(pts_ln,[x(j-1) y(j-1) g j-1 t(j-1) -1 a1 subtr isubtr]);
                isubtr=isubtr+1;
                imem=g;
                jmem=j-1;
            end
            pts_ln=vertcat(pts_ln,[x(j) y(j) g j t(j) a1 a2 subtr isubtr]);
            imem=g;
            jmem=j;
            
            rectr=vertcat(rectr,[x(j-1) y(j-1) g j-1 t(j-1) x(j) y(j) g j t(j) subtr isubtr]);
        else
            p1=[x(j-1),y(j-1)];
            p2=[x(j),y(j)];
            a1=rad2deg(angle2Points(p1,p2));
            
            p1=[x(j),y(j)];
            p2=[x(j+1),y(j+1)];
            a2=rad2deg(angle2Points(p1,p2));
            
            pts_ln=vertcat(pts_ln,[x(j) y(j) g j t(j) a1 a2 subtr isubtr]);
            imem=g;
            jmem=j;
            
            rectr=vertcat(rectr,[x(j) y(j) g j t(j) x(j+1) y(j+1) g j+1 t(j+1) subtr isubtr]);
        end
        
        xy=vertcat(xy,[x(j), y(j)]);
        p1=[x(j-1),y(j-1)];
        p2=[x(j),y(j)];
        p3=p2;
        p4=[x(j+1),y(j+1)];
        
        angle = lines_angle(p1,p2,p3,p4); % angle difference between incoming and outgoing edge of pt
        aangle = abs(angle); %absolute angle
        
        a1=rad2deg(angle2Points(p1,p2)); %incoming angle
        a2=rad2deg(angle2Points(p3,p4)); %outgoing angle
        
        angles=vertcat(angles,[rad2deg(angle2Points(p1,p2)) angle]);
        
        if ((aangle >= angle_tol && aangle <= 180 - angle_tol) || ...
                (aangle >= 180 + angle_tol && aangle <= 360 - angle_tol)) && ((dv(j)/dv(j-1)<=1) || dv(j) <= 40 && dv(j)>=1)
            my_corner_pts=vertcat(my_corner_pts,[x(j) y(j) g j x(j-1) y(j-1) x(j+1) y(j+1) a1 a2 dt(j-1) dt(j) mcnt dv(j-1) dv(j) subtr isubtr]);
            mcnt=mcnt+1;
        end
        isubtr=isubtr+1;
        dd=dd+polylineLength(xy,'open');
        total=total+length(x);
    end
end

slk=[];
vehicles=pts_ln(:,3);
vehicles=unique(vehicles,'rows');
for g=1:length(vehicles)
    tr=find(pts_ln(:,3)==vehicles(g));
    itr=unique(pts_ln(tr,8),'rows');
    
    for k=1:length(itr)
        tx=find(pts_ln(tr,8)==itr(k));
        if length(tx)<=1
            rh=find(my_corner_pts(:,3)==pts_ln(tr(tx),3) & my_corner_pts(:,16)==pts_ln(tr(tx),8) & my_corner_pts(:,17)==pts_ln(tr(tx),9));
            if ~isempty(rh)
                my_corner_pts(rh,:)=[];
            end
        else
            xy=[pts_ln(tr(tx),1) pts_ln(tr(tx),2)];
            
            for n=1:length(tx)-1
                xy=[pts_ln(tr(tx(n)),1) pts_ln(tr(tx(n)),2)];
                xy=vertcat(xy,[pts_ln(tr(tx(n+1)),1) pts_ln(tr(tx(n+1)),2)]);
                mllk=[xy(1,1) xy(1,2) xy(2,1) xy(2,2) polylineLength(xy,'open')];
                
                tllk=mllk;
                rj=1000;
                
                if ~isempty(tllk)
                    while length(rj)~=0
                        
                        rj=find(tllk(:,5)<=10);
                        ttllk=tllk(rj,:);
                        rj=find(tllk(:,5)>10);
                        for b=1:length(rj)
                            md=midPoint([[tllk(rj(b),1) tllk(rj(b),2)] [tllk(rj(b),3) tllk(rj(b),4)]]);
                            ttllk=vertcat(ttllk,[tllk(rj(b),1), tllk(rj(b),2),md(1,1),md(1,2), distancePoints([tllk(rj(b),1) tllk(rj(b),2)],[md(1,1) md(1,2)], 2)]);
                            ttllk=vertcat(ttllk,[md(1,1),md(1,2),tllk(rj(b),3), tllk(rj(b),4), distancePoints([md(1,1) md(1,2)],[tllk(rj(b),3) tllk(rj(b),4)], 2)]);
                        end
                        tllk=ttllk;
                    end
                    mllk=tllk;
                    slk=vertcat(slk,mllk);
                end
            end
        end
    end
end

% ----------------------------------------------------------
% finding cluster type of turn, range, cluster weight
% ----------------------------------------------------------

fprintf('Finding type of turns, cluster range and weight\n');

isim=[];
idis=[];
sim_angle=40;
for g=1:length(my_corner_pts(:,1))
    x1=my_corner_pts(g,1)+dist;
    y1=my_corner_pts(g,2)+dist;
    x2=my_corner_pts(g,1)+dist;
    y2=my_corner_pts(g,2)-dist;
    x3=my_corner_pts(g,1)-dist;
    y3=my_corner_pts(g,2)+dist;
    x4=my_corner_pts(g,1)-dist;
    y4=my_corner_pts(g,2)-dist;
    
    x=[x1 x2 x3 x4];
    x=vertcat(x,[y1 y2 y3 y4]);
    mb=minBoundingBox(x);
    mb=horzcat(mb,[mb(1,1); mb(2,1)]);
    in1 = inpolygon(slk(:,1), slk(:,2),mb(1,:), mb(2,:));
    in2 = inpolygon(slk(:,3), slk(:,4),mb(1,:), mb(2,:));
    
    rs1=find(in1==1);
    rs2=find(in2==1);
    rss=union(rs1,rs2,'rows');
    
    p1=0;
    p2=0;
    p3=0;
    p4=0;
    pts=[];
    rz=[];
    rpts=[];
    
    if my_corner_pts(g,14)<my_corner_pts(g,15) && my_corner_pts(g,14)>=1 || my_corner_pts(g,14)>my_corner_pts(g,15) && my_corner_pts(g,15)<=1
        in = createLine([my_corner_pts(g,5) my_corner_pts(g,6)],[my_corner_pts(g,1) my_corner_pts(g,2)]);
        e1 =orthogonalLine(in, [my_corner_pts(g,1) my_corner_pts(g,2)]);
        ndegrees=my_corner_pts(g,9);
        
        for h=1:length(rss)
            if slk(rss(h),1)~=my_corner_pts(g,1)&slk(rss(h),2)~=my_corner_pts(g,2)&slk(rss(h),3)~=my_corner_pts(g,1)&slk(rss(h),4)~=my_corner_pts(g,2)
                e2 = createLine([slk(rss(h),1) slk(rss(h),2)], [slk(rss(h),3) slk(rss(h),4)]);
                point = intersectLines(e1, e2);
                
                if ~isnan(point(1,1)) && ~isinf(point(1,1)) && ((abs(ndegrees-rad2deg(angle2Points([slk(rss(h),1) slk(rss(h),2)], [slk(rss(h),3) slk(rss(h),4)])))<=20||abs(ndegrees-rad2deg(angle2Points([slk(rss(h),1) slk(rss(h),2)], [slk(rss(h),3) slk(rss(h),4)])))>=340)||(abs(ndegrees-rad2deg(angle2Points([slk(rss(h),1) slk(rss(h),2)], [slk(rss(h),3) slk(rss(h),4)])))>=160&&abs(ndegrees-rad2deg(angle2Points([slk(rss(h),1) slk(rss(h),2)], [slk(rss(h),3) slk(rss(h),4)])))<=200))
                    rz=vertcat(rz,rss(h));
                    pts=vertcat(pts,point);
                end
            end
        end
        
        if ~isempty(pts)
            p1=max(distancePoints([my_corner_pts(g,1) my_corner_pts(g,2)], pts,2));
        end
    else
        in = createLine([my_corner_pts(g,1) my_corner_pts(g,2)],[my_corner_pts(g,7) my_corner_pts(g,8)]);
        e1 =orthogonalLine(in, [my_corner_pts(g,1) my_corner_pts(g,2)]);
        ndegrees=my_corner_pts(g,10);
        
        for h=1:length(rss)
            if slk(rss(h),1)~=my_corner_pts(g,1)&slk(rss(h),2)~=my_corner_pts(g,2)&slk(rss(h),3)~=my_corner_pts(g,1)&slk(rss(h),4)~=my_corner_pts(g,2)
                e2 = createLine([slk(rss(h),1) slk(rss(h),2)], [slk(rss(h),3) slk(rss(h),4)]);
                point = intersectLines(e1, e2);
                if ~isnan(point(1,1)) && ~isinf(point(1,1)) && ((abs(ndegrees-rad2deg(angle2Points([slk(rss(h),1) slk(rss(h),2)], [slk(rss(h),3) slk(rss(h),4)])))<=20||abs(ndegrees-rad2deg(angle2Points([slk(rss(h),1) slk(rss(h),2)], [slk(rss(h),3) slk(rss(h),4)])))>=340)||(abs(ndegrees-rad2deg(angle2Points([slk(rss(h),1) slk(rss(h),2)], [slk(rss(h),3) slk(rss(h),4)])))>=160&&abs(ndegrees-rad2deg(angle2Points([slk(rss(h),1) slk(rss(h),2)], [slk(rss(h),3) slk(rss(h),4)])))<=200))
                    pts=vertcat(pts,point);
                    rpts=vertcat(rpts,point);
                end
            end
        end
        
        if ~isempty(rpts)
            p2=max(distancePoints([my_corner_pts(g,1) my_corner_pts(g,2)], rpts,2));
        end
    end
    
    if ~isempty(pts)
        vpts=[];
        for h=1:length(pts(:,1))
            if distancePoints([my_corner_pts(g,1) my_corner_pts(g,2)], [pts(h,1) pts(h,2)], 2)<=dist
                vpts=vertcat(vpts,[pts(h,1) pts(h,2)]);
            end
        end
        
        pts=vpts;
        
        if ~isempty(pts)
            rv=find(pts(:,1)==my_corner_pts(g,1)&pts(:,2)==my_corner_pts(g,2));
            if isempty(rv)
                pts=vertcat(pts,[my_corner_pts(g,1) my_corner_pts(g,2)]);
            end
            pts=unique([pts(:,1) pts(:,2)],'rows');
            
            if length(pts(:,1))>1
                aa=centroid(pts);
                if distancePoints([my_corner_pts(g,1) my_corner_pts(g,2)], aa, 2)<=40
                    my_corner_pts(g,1) = aa(1,1);
                    my_corner_pts(g,2) = aa(1,2);
                    p1=max(p1,p2);
                end
            end
        end
    end
    
    my_corner_pts(g,9)=rad2deg(angle2Points([my_corner_pts(g,5) my_corner_pts(g,6)],[my_corner_pts(g,1) my_corner_pts(g,2)]));
    my_corner_pts(g,10)=rad2deg(angle2Points([my_corner_pts(g,1) my_corner_pts(g,2)],[my_corner_pts(g,7) my_corner_pts(g,8)]));
    
    if p1>1
        x1=my_corner_pts(g,1)+p1;
        y1=my_corner_pts(g,2)+p1;
        x2=my_corner_pts(g,1)+p1;
        y2=my_corner_pts(g,2)-p1;
        x3=my_corner_pts(g,1)-p1;
        y3=my_corner_pts(g,2)+p1;
        x4=my_corner_pts(g,1)-p1;
        y4=my_corner_pts(g,2)-p1;
    else
        x1=my_corner_pts(g,1)+15;
        y1=my_corner_pts(g,2)+15;
        x2=my_corner_pts(g,1)+15;
        y2=my_corner_pts(g,2)-15;
        x3=my_corner_pts(g,1)-15;
        y3=my_corner_pts(g,2)+15;
        x4=my_corner_pts(g,1)-15;
        y4=my_corner_pts(g,2)-15;
    end
    
    x=[x1 x2 x3 x4];
    x=vertcat(x,[y1 y2 y3 y4]);
    mb=minBoundingBox(x);
    mb=horzcat(mb,[mb(1,1); mb(2,1)]);
    
    in1 = inpolygon(slk(:,1), slk(:,2),mb(1,:), mb(2,:));
    in2 = inpolygon(slk(:,3), slk(:,4),mb(1,:), mb(2,:));
    
    rs1=find(in1==1);
    rs2=find(in2==1);
    rs=union(rs1,rs2,'rows');
    
    itype=1;
    dir=1;
    aaa=[];
    ii=[];
    if length(rs)>1
        
        for h=1:length(rs)
            cangle=rad2deg(angle2Points([slk(rs(h),1) slk(rs(h),2)],[slk(rs(h),3) slk(rs(h),4)]));
            aaa=vertcat(aaa,cangle);
        end
        [aaa ii]=unique(aaa,'rows');
        ii=sortrows(ii);
        
        idds=[];
        idds=vertcat(idds,[my_corner_pts(g,13) my_corner_pts(g,5) my_corner_pts(g,6) my_corner_pts(g,1) my_corner_pts(g,2) my_corner_pts(g,9) 1000 -1000 0 0]);
        for h=1:length(ii)
            cangle=rad2deg(angle2Points([slk(rs(ii(h)),1) slk(rs(ii(h)),2)],[slk(rs(ii(h)),3) slk(rs(ii(h)),4)]));
            
            cnt=0;
            for d=1:length(idds(:,1))
                if cangle>idds(d,6) && abs(idds(d,8)-cangle)<=10 && idds(d,8)~=-1000
                    cnt=cnt+1;
                    idds(d,9)=idds(d,9)+1;
                    
                    if cangle > idds(d,8)
                        idds(d,8)=cangle;
                    end
                elseif cangle<idds(d,6) && abs(idds(d,7)-cangle)<=10 && idds(d,7)~=1000
                    cnt=cnt+1;
                    idds(d,9)=idds(d,9)+1;
                    if cangle < idds(d,7)
                        idds(d,7)=cangle;
                    end
                    
                elseif abs(idds(d,6)-cangle)<=sim_angle
                    cnt=cnt+1;
                    idds(d,9)=idds(d,9)+1;
                    if cangle > idds(d,8) && cangle > idds(d,6)
                        idds(d,8)=cangle;
                    elseif cangle < idds(d,7) && cangle < idds(d,6)
                        idds(d,7)=cangle;
                    end
                end
            end
            if cnt==0
                idds=vertcat(idds,[my_corner_pts(g,13) slk(rs(ii(h)),1) slk(rs(ii(h)),2) slk(rs(ii(h)),3) slk(rs(ii(h)),4) cangle 1000 -1000 0 0]);
            end
        end
        
        for h=1:length(idds(:,1))
            idis=vertcat(idis,[idds(h,1) idds(h,2) idds(h,3) idds(h,4) idds(h,5) idds(h,6) idds(h,7) idds(h,8) idds(h,9)]);
        end
        
        if ~isempty(idis)
            if length(idds)>0
                itype=length(idds(:,1));
            end
        end
    end
    
    my_corner_pts(g,18)=itype;
    my_corner_pts(g,19)=p1;
    my_corner_pts(g,20)=length(ii);
end


temp_my_corner_pts=sortrows(my_corner_pts,-19);

my_corner_pts=temp_my_corner_pts;
g=1;
idx=1;
while g<=length(my_corner_pts(:,1))
    p1=0;
    if my_corner_pts(g,19)~=0
        if my_corner_pts(g,19)>1
            p1=my_corner_pts(g,19);
        else
            p1=dist/2;
        end
    else
        p1=dist/2;
    end
    
    x1=my_corner_pts(g,1)+p1;
    y1=my_corner_pts(g,2)+p1;
    x2=my_corner_pts(g,1)+p1;
    y2=my_corner_pts(g,2)-p1;
    x3=my_corner_pts(g,1)-p1;
    y3=my_corner_pts(g,2)+p1;
    x4=my_corner_pts(g,1)-p1;
    y4=my_corner_pts(g,2)-p1;
    
    x=[x1 x2 x3 x4];
    x=vertcat(x,[y1 y2 y3 y4]);
    mb=minBoundingBox(x);
    mb=horzcat(mb,[mb(1,1); mb(2,1)]);
    
    in = inpolygon(my_corner_pts(:,1), my_corner_pts(:,2),mb(1,:), mb(2,:));
    rs=find(in==1);
    rc=find(my_corner_pts(rs,13)~=my_corner_pts(g,13));
    rs=rs(rc);
    
    sim=[];
    for h=1:length(rs)
        if abs(my_corner_pts(g,9)-my_corner_pts(rs(h),9))<=15 ||abs(my_corner_pts(g,9)+my_corner_pts(rs(h),9)-360)<=15 || abs(my_corner_pts(g,9)-my_corner_pts(rs(h),10))>=165&&abs(my_corner_pts(g,9)-my_corner_pts(rs(h),10))<=195
            sim=vertcat(sim,h);
        elseif abs(my_corner_pts(g,10)-my_corner_pts(rs(h),10))<=15 ||abs(my_corner_pts(g,10)+my_corner_pts(rs(h),10)-360)<=15 || abs(my_corner_pts(g,10)-my_corner_pts(rs(h),9))>=165&&abs(my_corner_pts(g,10)-my_corner_pts(rs(h),9))<=195
            sim=vertcat(sim,h);
        end
    end
    
    rest=setdiff(rs,rs(sim));
    ra=find(idis(:,1)==my_corner_pts(g,13));
    for h=1:length(rest)
        for l=1:length(ra)
            if abs(idis(ra(l),6)-my_corner_pts(rest(h),9))<=idis(ra(l),6)-idis(ra(l),7)
                ih=find(rs(:,1)==rest(h));
                sim=vertcat(sim,ih);
            end
        end
    end
    
    sim=unique(sim,'rows');
    
    if length(sim)~=0
        
        pts_to_inter=vertcat(pts_to_inter,[1 1 my_corner_pts(g,3) my_corner_pts(g,4) my_corner_pts(g,13) my_corner_pts(g,1) my_corner_pts(g,2) my_corner_pts(g,5) my_corner_pts(g,6) my_corner_pts(g,7) my_corner_pts(g,8) idx]);
        inter_to_cluster=vertcat(inter_to_cluster, [my_corner_pts(g,3), my_corner_pts(g,4), idx]);
        mpts=[my_corner_pts(g,5) my_corner_pts(g,6) my_corner_pts(g,1) my_corner_pts(g,2)];
        mpts=vertcat(mpts,[my_corner_pts(g,1) my_corner_pts(g,2) my_corner_pts(g,7) my_corner_pts(g,8)]);
        
        xy=[my_corner_pts(g,1) my_corner_pts(g,2)];
        for h=1:length(sim)
            xy=vertcat(xy,[my_corner_pts(rs(sim(h)),1) my_corner_pts(rs(sim(h)),2)]);
            
            mpts=vertcat(mpts,[my_corner_pts(rs(sim(h)),5) my_corner_pts(rs(sim(h)),6) my_corner_pts(rs(sim(h)),1) my_corner_pts(rs(sim(h)),2)]);
            mpts=vertcat(mpts,[my_corner_pts(rs(sim(h)),1) my_corner_pts(rs(sim(h)),2) my_corner_pts(rs(sim(h)),7) my_corner_pts(rs(sim(h)),8)]);
            
            pts_to_inter=vertcat(pts_to_inter,[1 g my_corner_pts(rs(sim(h)),3) my_corner_pts(rs(sim(h)),4) my_corner_pts(rs(sim(h)),13) my_corner_pts(rs(sim(h)),1) my_corner_pts(rs(sim(h)),2) my_corner_pts(rs(sim(h)),5) my_corner_pts(rs(sim(h)),6) my_corner_pts(rs(sim(h)),7) my_corner_pts(rs(sim(h)),8) idx]);
            inter_to_cluster=vertcat(inter_to_cluster, [my_corner_pts(rs(sim(h)),3), my_corner_pts(rs(sim(h)),4), idx]);
        end
        aa=centroid(xy);
        
        ipts=[aa(1,1) aa(1,2)];
        
        for b=1:2:length(mpts(:,1))-1
            in1 = createLine([mpts(b,1) mpts(b,2)],[mpts(b,3) mpts(b,4)]);
            d1=rad2deg(angle2Points([mpts(b,1) mpts(b,2)],[mpts(b,3) mpts(b,4)]));
            in2 = createLine([mpts(b+1,1) mpts(b+1,2)],[mpts(b+1,3) mpts(b+1,4)]);
            d2=rad2deg(angle2Points([mpts(b+1,1) mpts(b+1,2)],[mpts(b+1,3) mpts(b+1,4)]));
            ipts=vertcat(ipts,[mpts(b,3) mpts(b,4)]);
            for c=b+2:2:length(mpts(:,1))-2
                in3 = createLine([mpts(c,1) mpts(c,2)],[mpts(c,3) mpts(c,4)]);
                d3=rad2deg(angle2Points([mpts(c,1) mpts(c,2)],[mpts(c,3) mpts(c,4)]));
                in4 = createLine([mpts(c+1,1) mpts(c+1,2)],[mpts(c+1,3) mpts(c+1,4)]);
                d4=rad2deg(angle2Points([mpts(c+1,1) mpts(c+1,2)],[mpts(c+1,3) mpts(c+1,4)]));
                
                if (abs(d1-d3)>=25&&abs(d1-d3)<=320)
                    point = intersectLines(in1, in3);
                    if ~isnan(point(1,1)) && ~isinf(point(1,1)) && distancePoints([aa(1,1) aa(1,2)], [point(1,1),point(1,2)], 2)<=p1
                        ipts=vertcat(ipts,[point(1,1),point(1,2)]);
                    end
                end
                if (abs(d1-d4)>=25&&abs(d1-d4)<=320)
                    point = intersectLines(in1, in4);
                    if ~isnan(point(1,1)) && ~isinf(point(1,1)) && distancePoints([aa(1,1) aa(1,2)], [point(1,1),point(1,2)], 2)<=p1
                        ipts=vertcat(ipts,[point(1,1),point(1,2)]);
                    end
                end
                if (abs(d2-d3)>=25&&abs(d2-d3)<=320)
                    point = intersectLines(in2, in3);
                    if ~isnan(point(1,1)) && ~isinf(point(1,1)) && distancePoints([aa(1,1) aa(1,2)], [point(1,1),point(1,2)], 2)<=p1
                        ipts=vertcat(ipts,[point(1,1),point(1,2)]);
                    end
                end
                if (abs(d2-d4)>=45&&abs(d2-d4)<=320)
                    point = intersectLines(in2, in4);
                    if ~isnan(point(1,1)) && ~isinf(point(1,1)) && distancePoints([aa(1,1) aa(1,2)], [point(1,1),point(1,2)], 2)<=p1
                        ipts=vertcat(ipts,[point(1,1),point(1,2)]);
                    end
                end
            end
        end
        
        aaa=centroid(ipts);
        aa=aaa;
        turn_pts=vertcat(turn_pts,[aa(1,1), aa(1,2), idx, my_corner_pts(g,9), my_corner_pts(g,10), length(sim(:,1))+1, 1]);
        turn_seeds=vertcat(turn_seeds,[idx g aa(1,1) aa(1,2) length(sim(:,1))+1 my_corner_pts(g,9) my_corner_pts(g,10) 1]);
        
        centre_pts=vertcat(centre_pts,[aa(1,1), aa(1,2), idx, my_corner_pts(g,9), my_corner_pts(g,10), (length(sim(:,1))+1), my_corner_pts(g,18), my_corner_pts(g,19),my_corner_pts(g,20)]);
        clusterseeds=vertcat(clusterseeds,[idx i aa(1,1) aa(1,2) length(sim(:,1))+1 my_corner_pts(g,9), my_corner_pts(g,10)]);
        
        my_corner_pts(rs(sim),:)=[];
        my_corner_pts(g,:)=[];
        idx=idx+1;
    else
        rz=vertcat(rz,my_corner_pts(g,13));
        g=g+1;
    end
end

save('slk.mat','slk');
save('pts_ln.mat','pts_ln');
save('vehicles.mat','vehicles');
save('centre_pts.mat','centre_pts');
save('pts_to_inter.mat','pts_to_inter');
save('clusterseeds.mat','clusterseeds');
save('inter_to_cluster.mat','inter_to_cluster');

fprintf('Finishing at %s\n', timestamp);

