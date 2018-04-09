function cal_partial_pressure

close all;

p=figure;
phi=figure;
sca=figure;

for i=1:10:801
    

[data, info] = load_stat_file('segregation.1.large',i);
[data, info] = get_standard_variables(data,info);


large_pressure=data.variables(:,12);
large_PVF=data.variables(:,1);

[data, info] = load_stat_file('segregation.1.small',i);
[data, info] = get_standard_variables(data,info);
small_pressure=data.variables(:,12);
small_PVF=data.variables(:,1);

[data, info] = load_stat_file('segregation.1.all',i);
[data, info] = get_standard_variables(data,info);
all_pressure=data.variables(:,12);
all_PVF=data.variables(:,1);

z=data.coordinates(:,3);

small_partial_pressure=small_pressure./all_pressure;
large_partial_pressure=large_pressure./all_pressure;

small_phi=small_PVF./all_PVF;
large_phi=large_PVF./all_PVF;

figure(p);
plot(z,[small_partial_pressure large_partial_pressure]);

figure(phi);
plot(z,[small_phi large_phi]);

figure(sca);
scatter(small_phi,small_partial_pressure,'b');
hold on
scatter(large_phi,large_partial_pressure,'r');

end

saveas(3 ,['figure1.pdf']);




end
