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
    Fwall.planes(end+1).normal=[0 0 1];
    Fwall.planes(end).position=rs.zmax;
    Fwall.planes(end+1).normal=[0 0 -1];
    Fwall.planes(end).position=-rs.zmin;
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

    A=[];
    b=[];
    for j=1:Fwall.num_planes-6
        A(j,:)=Fwall.planes(j).normal;
        b(j,1)=Fwall.planes(j).position;
    end
%     for j=Fwall.num_planes-5:Fwall.num_planes
%         A(j,:)=-Fwall.planes(j).normal;
%         b(j,1)=-Fwall.planes(j).position;
%     end
    %check if A*IP>=b 
    IPinternal=[];
    for j=1:length(IP)
        if (min(A*IP(j).point-b)>=-1e-12)
            IPinternal(end+1).point=IP(j).point;
        end
    end 
   

    X=[];
    for j=1:length(IPinternal)
        X(end+1,:) = IPinternal(j).point;
    end

    K=convhulln(X)
    gs=1.0/rs.NFwalls*i;
    trisurf(K,X(:,1),X(:,2),X(:,3),'FaceAlpha',.2,'FaceColor',[gs gs gs],'EdgeColor',[gs gs gs])
    hold on

    disp('done')

end
axis equal