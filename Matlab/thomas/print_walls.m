function rs=plotRestartData(rs)
% rs is the restart file which will be plotted

%load the restart data file in
if (~exist('rs','var'))
    files=dir('*.restart');
    rs=load_restartdata(files.name);
end

figure(1); clf;


%plot the particles in red
ind_flow = sum(rs.Particles(:,4:6).^2,2)>0;
plot3(rs.Particles(~ind_flow,1),rs.Particles(~ind_flow,2),rs.Particles(~ind_flow,3),'k.')
hold on
plot3(rs.Particles(ind_flow,1),rs.Particles(ind_flow,2),rs.Particles(ind_flow,3),'r.')
set(gcf,'Position',[1357         672         560         420]);

% Label the axis for easier orientation
xlabel('x')
ylabel('y')
zlabel('z')

%Now plot the infinte walls
for i=1:rs.NIwalls
    
    normal=rs.Iwalls(i).normal;
    position=rs.Iwalls(i).position;
    
    x=[1;1]*[rs.xmin rs.xmax];
    y=[rs.ymin;rs.ymax]*[1 1];
    z=(position-normal(1)*x-normal(2)*y)/normal(3);
    surf(x,y,z);
end

%% Now do the finite walls



for i=1:rs.NFwalls
    
    % add domain boundaries
    Fwall = rs.Fwalls(i);   
    Fwall.planes(end+1).normal=[1 0 0];
    Fwall.planes(end).position=rs.xmax;
    Fwall.planes(end+1).normal=[-1 0 0];
    Fwall.planes(end).position=-rs.xmin;
    Fwall.planes(end+1).normal=[0 1 0];
    Fwall.planes(end).position=rs.ymax;
    Fwall.planes(end+1).normal=[0 -1 0];
    Fwall.planes(end).position=-rs.ymin;
    Fwall.planes(end+1).normal=[0 0 -1];
    Fwall.planes(end).position=-rs.zmin;
    Fwall.planes(end+1).normal=[0 0 1];
    Fwall.planes(end).position=rs.zmax;
    Fwall.num_planes = Fwall.num_planes + 6;
    
    %store intersection points in IP
    IP=[];
    no=0;
    for j1=1:Fwall.num_planes
    for j2=j1+1:Fwall.num_planes
    for j3=j2+1:Fwall.num_planes
        A=[Fwall.planes(j1).normal; Fwall.planes(j2).normal; Fwall.planes(j3).normal];
        if (abs(det(A))>1e-12) %if non-singular
            b=[Fwall.planes(j1).position; Fwall.planes(j2).position; Fwall.planes(j3).position];
            IP(end+1).point = A\b;
        else
            no=no+1;
        end
    end
    end
    end
    
    %Check all points, if they are ib the correct side of the finite wall.
    %Some of the intersections will be fake as in about the finite walls
    %has already been ended by another wall.

    A=[];
    b=[];
    % No not need to plot the domain walls, so do not do this step for them
    for j=1:Fwall.num_planes-6
        A(j,:)=Fwall.planes(j).normal;
        b(j,1)=Fwall.planes(j).position;
    end
    for j=Fwall.num_planes-5:Fwall.num_planes
        A(j,:)=-Fwall.planes(j).normal;
        b(j,1)=-Fwall.planes(j).position;
    end
    %check if A*IP>=b i.e. point is on the inside of the wall
    IPinternal=[];
    for j=1:length(IP)
        if (min(A*IP(j).point-b)>=-1e-12)
            IPinternal(end+1).point=IP(j).point;
        end
    end 
   
    %Extract the X location of the true internal intersects.
    X=[];
    for j=1:length(IPinternal)
        X(end+1,:) = IPinternal(j).point;
    end

    %Plot the surface as a convex hull, making each wall a different shade
    %of gray.
    
    K=convhulln(X);
    gs=0.5/rs.NFwalls*i;
    trisurf(K,X(:,1),X(:,2),X(:,3),'FaceAlpha',.2,'FaceColor',[gs gs gs],'EdgeColor',[gs gs gs])
    hold on

    disp('done')
    
    %keyboard;

end
axis equal

return