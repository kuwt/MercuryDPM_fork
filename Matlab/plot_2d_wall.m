function plot_2d_wall(s)
A=sscanf(s,'%*s %*i %*s %f %f %f %*s %f %*s %f %f %f %*s %f %*s %f %f %f %*s %f',12);

n1=[A(1) A(3)];
n2=[A(5) A(7)];
n3=[A(9) A(11)];
p1=A(4);
p2=A(8);
p3=A(12);

N=[n1;n2];
disp('N'),disp(N)
P=[p1;p2];
point1=N\P;
disp('point1'), disp(point1)
plot(point1(1),point1(2),'x');
N=[n2;n3];
P=[p2;p3];
point2=N\P;
disp('point2'), disp(point2)
plot(point2(1),point2(2),'x');
N=[n3;n1];
P=[p3;p1];
point3=N\P;
disp('point3'), disp(point3)
plot(point3(1),point3(2),'x');
line([point1(1), point2(1)],[point1(2),point2(2)]);
line([point2(1), point3(1)],[point2(2),point3(2)]);
line([point3(1), point1(1)],[point3(2),point1(2)]);
