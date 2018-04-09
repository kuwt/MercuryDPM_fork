%input
name = '/Users/weinhartt/Work/code/MD/DRIVERS/Thomas/run/AngledPeriodicBoundaryTest/AngledPeriodicBoundaryTestC.restart';

data=load_restartdata(name)
D=max([data.xmax-data.xmin,data.ymax-data.ymin,data.zmax-data.zmin]); %domainsize
gray=0.5*[1 1 1];
addpath('~/Work/code/MD/matlab/')
fstat=read_fstat([name(1:end-8) '.fstat']);

figure(1); clf
hold all

%% plot boundaries
if (data.Boundaries>0)
    if isfield(data,'AngledPeriodicBoundary')
        for i=1:size(data.AngledPeriodicBoundary.origin,1)
            normal_left=data.AngledPeriodicBoundary.normal_left(i,:);
            normal_right=data.AngledPeriodicBoundary.normal_right(i,:);
            origin=data.AngledPeriodicBoundary.origin(i,:);
            common_axis = cross(normal_left, normal_right);
            common_axis = common_axis/norm(common_axis);
            radialAxis_left = cross(normal_left, common_axis);
            radialAxis_right = cross(normal_right, common_axis);
            nodes = [origin;origin+D*radialAxis_left;origin+D*common_axis+D*radialAxis_left;origin+D*common_axis];
            h=fill3(nodes(:,1),nodes(:,2),nodes(:,3),gray,'EdgeColor',gray,'LineStyle',':');
            nodes = [origin;origin+D*radialAxis_right;origin+D*common_axis+D*radialAxis_right;origin+D*common_axis];
            fill3(nodes(:,1),nodes(:,2),nodes(:,3),gray,'EdgeColor',gray,'LineStyle',':');
        end
    end
end

%% plot particles (2D)
t=linspace(0,2*pi,40);
x=sin(t);
y=cos(t);
n=linspace(0,1,5);
dts=rmax/max(abs(fstat.forcet)); %scale
if (data.Particles>0)
    if isfield(data,'BP')
        vmax=sqrt(max(sum(data.BP(:,4:6).^2,2)));
        rmax=max(data.BP(:,7));
        vs=rmax/vmax; %scale
        for i=1:size(data.BP,1)
            p=data.BP(i,1:3);
            r=data.BP(i,7);
            plot(p(1)+r*x,p(2)+r*y,'b');
            v=data.BP(i,4:6);
            plot(p(1),p(2),'b.');
            plot(p(1)+[0 vs*v(1)],p(2)+[0 vs*v(2)],'b');
            a=data.BP(i,8:10);
            w=data.BP(i,11:13);
            plot(p(1)+0.8*r*sin(a(3)),p(2)+0.8*r*cos(a(3)),'b.');
            plot(p(1)+0.8*r*sin(a(3)+n*w(3)),p(2)+0.8*r*cos(a(3)+n*w(3)),'b');
            if ~isempty(data.TSPSpring{i})
                plot(p(1)+dts*[0 data.TSPSpring{i}(1)],p(2)+dts*[0 data.TSPSpring{i}(2)],'r');
            end
        end
    end
end

%% plot springs
ix=find(fstat.t==max(fstat.t));
dts=rmax/max(abs(fstat.forcet)); %scale
ds=rmax/max(abs(fstat.forcen)); %scale
for i=ix'
    c=fstat.centre(i,:);
    c=data.BP(fstat.PJ(i,:)+1,1:3);
    forcen = fstat.forcen(i);
    n=fstat.normal(i,:);
    plot(c(1)+[0 ds*forcen*n(1)],c(2)+[0 ds*forcen*n(2)],'g:')
    forcet = fstat.forcet(i);
    t=fstat.tangential(i,:);
    plot(c(1)+[0 dts*forcet*t(1)],c(2)+[0 dts*forcet*t(2)],'g:')
end

%% 
axis equal
xlim([data.xmin data.xmax+eps])
ylim([data.ymin data.ymax+eps])
zlim([data.zmin data.zmax+eps])
xlabel('x')
ylabel('y')
zlabel('z')
view(0,90)
