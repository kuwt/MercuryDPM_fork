%% open files
clear weight samplingRate name
files = dir('Scale/*.mat');
fprintf("\n%20s\t%8s\n","File name","Throughput [kg/s]")
for i=1:length(files)
    name = [files(i).folder '/' files(i).name];
    open(name);
    d(i).weight = (ans.weight(2:end-1)-ans.weight(2))/1000;
    d(i).samplingRate = ans.sampling_rate;
    d(i).name = files(i).name(1:end-10);
    d(i).maxTime = (length(d(i).weight)-1)/d(i).samplingRate;
    d(i).time = linspace(0,d(i).maxTime,length(d(i).weight));
    d(i).meanThroughput = d(i).weight(end)/d(i).maxTime*3600; %[kg/h]
    fprintf("%20s\t%8f\n",d(i).name,d(i).meanThroughput)
end
% select specific values and rescale them as needed
% we need: MPT 0.12 0.24 6
%          APAP 8.4
%          Avicel MCC PH101 3.6
%          Pearlitol SD100 11.88, 23.76
d=d([2 3 6 9 11 14 15]);
meanThroughput=[0.12 0.24 6 8.4 3.6 .99*12 .99*24];
for i=1:length(d)
    d(i).weight = d(i).weight * meanThroughput(i)/d(i).meanThroughput;
    d(i).meanThroughput = meanThroughput(i);
end
%% plot accumulated weight
for i=1:length(d)
    plot(d(i).time,d(i).weight,'.-','Displayname',d(i).name)
    hold on
end
xlabel('time [s]')
ylabel('weight [kg]')
hold off
legend show
legend location eastoutside
%saveas(gcf,'mass.png')
%% plot throughput
for i=1:length(d)
    plot(d(i).time,smooth(gradient(d(i).weight,d(i).time)*3600,500),'.-','Displayname',[d(i).name num2str(d(i).meanThroughput,', %.3f kg/h')])
    hold on
end
xlabel('time [s]')
ylabel('throughput [kg/h] (smoothed over 5s)')
axis tight
ylim([0 25])
hold off
legend show
legend location eastoutside
set(legend,'interpreter','none')
%saveas(gcf,'throughput.png')
%% output for c++
% we need: MPT 0.12 0.24 6
%          APAP 8.4
%          Pearlitol SD100 11.88, 23.76
%          Avicel MCC PH101 3.6
fid = fopen('inflowData.h','w');
for i = 1:length(d)
fprintf(fid,'const std::vector<Mdouble> %sWeight={%s};\n', ...
    d(i).name, num2str(d(i).weight,'%g,'));
fprintf(fid,'const Mdouble %sSamplingInterval = %g;\n',...
    d(i).name, 1/d(i).samplingRate);
end
fclose(fid);
type inflowData.h

%% plot accumulated weight
d5weight = interp1(d(4).time,d(4).weight,d(5).time);
plot(d(5).time,0.3*d(5).weight,'.-','Displayname',d(4).name)
hold on
plot(d(5).time,0.7*d5weight,'.-b','Displayname',d(5).name)
plot(d(5).time,0.3*d(5).weight+0.7*d5weight,'.-','Displayname','total')
xlabel('time [s]','Interpreter','none')
ylabel('weight [kg]','Interpreter','none')
hold off
legend show
legend location eastoutside 
set(legend,'Interpreter','none')
xlim([0,5])
saveas(gcf,'nonMononotousInflow.png')
