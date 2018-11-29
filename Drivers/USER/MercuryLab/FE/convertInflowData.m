open pearlitol12-Scale.mat
weight = ans.weight(2:end)-ans.weight(2);
samplingRate = ans.sampling_rate
maxTime = (length(weight)-1)/samplingRate;
time = linspace(0,maxTime,length(weight))
plot(time,weight,'.-')

%% output for c++
fid = fopen('inflowData.h','w');
fprintf(fid,'const std::vector<Mdouble> weight={%s};\n',...
    num2str(weight,'%g,'));
fprintf(fid,'const Mdouble samplingRate = %g;\n',...
    samplingRate);
fclose(fid);
type inflowData.h

