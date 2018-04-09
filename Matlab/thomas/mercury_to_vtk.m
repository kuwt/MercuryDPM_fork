function mercury_to_vtk(filename)

%Let's read data! (ONLY FORMAT 14)
fid=fopen(filename);
if (fid==-1)
  filename_read = [filename '.data']; 
  fid=fopen(filename_read);
end

%write data
data=[];
format = 14;
while 1
  switch format
    case 14
      header=textscan(fid,'%f %f %f %f %f %f %f %f',1);
      if isempty(header{1}); break; end
      rawdata=textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f',header{1});
      %write it into structure
      data1.N = header{1};
      data1.t = header{2};
      data1.xmin = header{3};
      data1.ymin = header{4};
      data1.zmin = header{5};
      data1.xmax = header{6};
      data1.ymax = header{7};
      data1.zmax = header{8};
      data1.Position = [rawdata{1} rawdata{2} rawdata{3}];
      data1.Velocity = [rawdata{4} rawdata{5} rawdata{6}];
      data1.Radius = rawdata{7};
      data1.Angle = [rawdata{8} rawdata{9} rawdata{10}];
      data1.AngularVelocity = [rawdata{11} rawdata{12} rawdata{13}];
      data1.info = rawdata{14};
  end
  data{end+1}=data1;
  %data{end}.t
end
fclose(fid);


%Let's write data!
path=pwd;
%rmdir('VTK_Output','s');
newdir = ['VTK_Output/' filename]; %makes new subdirectory for each case
mkdir(newdir); 

for i=1:size(data, 2)
    
    data{1,i}.Position(data{1,i}.Position~=data{1,i}.Position)=0;
    data{1,i}.Velocity(data{1,i}.Velocity~=data{1,i}.Velocity)=0; %cleans NaN - not supported in vtk legacy
    data{1,i}.AngularVelocity(data{1,i}.AngularVelocity~=data{1,i}.AngularVelocity)=0;
    
    filenumber = int2str(i);
    filename_write = [filename filenumber '.vtk']; % separate files for each timestep
    
    fid = fopen(fullfile(path,newdir,filename_write), 'wt');
    fprintf(fid, '# vtk DataFile Version 2.0\n');
    fprintf(fid, 'Mercury MD Output to PARAVIEW\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, '\n');
    
    fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
    fprintf(fid, 'POINTS %d float\n', data{1,i}.N);    
    for j=1:data{1,i}.N
    fprintf(fid, '%f %f %f\n', data{1,i}.Position(j,1), data{1,i}.Position(j,2),data{1,i}.Position(j,3));
    end
    
    fprintf(fid, 'POINT_DATA %d\n', data{1,i}.N);
    
    fprintf(fid, 'SCALARS Radius float\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    for j=1:data{1,i}.N
    fprintf(fid, '%f\n', data{1,i}.Radius(j,1));    
    fprintf(fid, '\n');
    end
    
    fprintf(fid, 'VECTORS Velocity float\n');
    for j=1:data{1,i}.N
    fprintf(fid, '%f %f %f\n', data{1,i}.Velocity(j,1), data{1,i}.Velocity(j,2),data{1,i}.Velocity(j,3));    
    fprintf(fid, '\n');
    end
    
    fprintf(fid, 'VECTORS Angular_Velocity float\n');
    for j=1:data{1,i}.N
    fprintf(fid, '%f %f %f\n', data{1,i}.AngularVelocity(j,1), data{1,i}.AngularVelocity(j,2),data{1,i}.AngularVelocity(j,3));    
    fprintf(fid, '\n');
    end
    
    %easy to add fields, just follow template
    
    fclose(fid);
    
end
    
return
