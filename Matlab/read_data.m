function data = read_data(filename,format,time)

%if no filename is given, pick the first data file in the directory
if (~exist('filename','var'))
  filenames = ls('*.data -t');
  filename = strtok(filenames);
  disp(filename);
end

%if no format is given, pick format 14
if (~exist('format','var'))
  format = 14;
end

%if no format is given, pick format 14
if (~exist('time','var'))
  time = [-inf inf];
end

%write data
data=[];
fid=fopen(filename);
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
  if (data1.t>=min(time) & data1.t<=max(time))
      data{end+1}=data1;
      data{end}.t
  end
end

fclose(fid);


return