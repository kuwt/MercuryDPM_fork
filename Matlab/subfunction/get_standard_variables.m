function [data, info] = get_standard_variables(data,info)

if isfield(data,'variance'), data = rmfield(data,'variance'); end

info.variable_names = {'PVF'; ...
  'VelocityX'; 'VelocityY'; 'VelocityZ'; ...
  'StressXX'; 'StressXY'; 'StressXZ'; 'StressYY'; 'StressYZ'; 'StressZZ'; ...
  'Temperature'; ...
  'Pressure'};
newvariables = [ ...
  data.variables(:,1), ...
  data.variables(:,[3:5])./data.variables(:,[2 2 2]), ...
  data.variables(:,[15:20])+data.variables(:,[21:26]), ...
  ];
data.variables = [ newvariables ...
  sum(data.variables(:,[6 9 11])./data.variables(:,[2 2 2])-newvariables(:,[2:4]).^2 ,2)/3, ...
  sum(newvariables(:,[5 8 10]),2)/3 ...
  ];


%top-right: add in the empty spot the deviator stress ratio:
%sD1:=(S1-S3)/(2*p)  or sD2 (to be discussed)
%in one of the other empty spots:

%pls. add the angle of the eigen-direction of stress from the 
%x-direction (parallel to the inclined plane?)
data.variables(:,end+(1:4)) = 0;
%data.variables(:,end-1) = (data.variables(:,5)-data.variables(:,10))./data.variables(:,12)/2;
stress = zeros(3);
for j=1:size(data.variables,1)
  stress(:) = data.variables(j,[5 6 7 6 8 9 7 9 10]);
  %stress = stress - eye(3)*trace(stress)/3;
  if (stress==[0 0 0;0 0 0;0 0 0])
    data.variables(j,end) = NaN;
  else
    [v1,d1]=eig(stress);
    %v2=v1*diag(sign(sum(v1(1:2,:),1)));
    [v,d]=dsort_ev(v1,d1);
    data.variables(j,end-3) = (max(diag(d))-min(diag(d)))/2./data.variables(j,12);
    angle = acos([1 0 0]*v)/pi*180;
    angle = min(angle,180-angle);
    [val,ind] = min(diag(d));
    data.variables(j,end-2:end) = angle;
  end
end
info.variable_names{end+1} = 'deviator stress ratio sD1';
info.variable_names{end+1} = 'angle max stress eigendir to xaxis';
info.variable_names{end+1} = 'angle middle stress eigendir to xaxis';
info.variable_names{end+1} = 'angle min stress eigendir to xaxis';

return