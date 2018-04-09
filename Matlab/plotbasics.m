% Plots some of the basic (=microscopic) fields in data
% 
function data = plotbasics(data)

if (isstr(data))
  data = loadstatistics(data);
end

for i=1:16
	subplot(4,4,i);
	switch i
		case 1
			fieldname = 'Nu';
		case 2
			fieldname = 'MomentumX';
		case 3
			fieldname = 'MomentumZ';
		case 4
			fieldname = 'MomentumFluxXX';
		case 5
			fieldname = 'MomentumFluxZZ';
		case 6
			fieldname = 'EnergyFluxX';
		case 7
			fieldname = 'EnergyFluxZ';
		case 8
			fieldname = 'NormalStressXX';
		case 9
			fieldname = 'NormalStressZZ';
		case 10
			fieldname = 'TangentialStressXX';
		case 11
			fieldname = 'TangentialStressZZ';
		case 12
			fieldname = 'FabricXZ';
		case 13
			fieldname = 'CoordinationNumber';
		otherwise
			fieldname = 'z';
	end
  plot(data.z,getfield(data,fieldname));
	axis tight
	title(fieldname);
end

return

%fieldnames(data{1})
