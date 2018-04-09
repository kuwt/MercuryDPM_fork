function data=plotgeometry(problemname)

%open data file
datafilename = [problemname '.data'];
datafile = fopen(datafilename);

%read header line
header=sscanf(fgetl(datafile),'%f')

%read particles
[X,Y,Z]=sphere(20)
figure; hold on; xlabel('x'); ylabel('y'); zlabel('z'); view(36,18)
if length(header)==8
  dim=3;
  for i=1:header(1)
    particle=sscanf(fgetl(datafile),'%f')
    x=particle(1);
    y=particle(2);
    z=particle(3);
    r=particle(7);
    surf(X*r+x,Y*r+y,Z*r+z);
  end
else
  dim=2;
end

header=fgetl(datafile)

fclose(datafile);

return