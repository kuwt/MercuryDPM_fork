function NablaVariables = nabla(variables, coordinates)

% stores the grid dimensions
for i=1:3
  n(i) = size(coordinates,1) / sum(coordinates(1,i)==coordinates(:,i));
end

% stores coordinates in grid
for i=1:3
  X{i} = zeros(n(end:-1:1));
  X{i}(:) = coordinates(:,i);
end

dX{1} = max(1e-60,(X{1}(:,:,[2:end end])-X{1}(:,:,[1 1:end-1])));
dX{2} = max(1e-60,(X{2}(:,[2:end end],:)-X{2}(:,[1 1:end-1],:)));
dX{3} = max(1e-60,(X{3}([2:end end],:,:)-X{3}([1 1:end-1],:,:)));

NablaVariables{1} = zeros(size(variables));
NablaVariables{2} = zeros(size(variables));
NablaVariables{3} = zeros(size(variables));
for j=1:size(variables,2)
  variablesGrid{j} = zeros(n(end:-1:1));
  variablesGrid{j}(:) = variables(:,j);

  dummy = (variablesGrid{j}(:,:,[2:end end])-variablesGrid{j}(:,:,[1 1:end-1]))./dX{1};
  NablaVariables{1}(:,j) = dummy(:);
  dummy= (variablesGrid{j}(:,[2:end end],:)-variablesGrid{j}(:,[1 1:end-1],:))./dX{2};
  NablaVariables{2}(:,j) = dummy(:);
  dummy= (variablesGrid{j}([2:end end],:,:)-variablesGrid{j}([1 1:end-1],:,:))./dX{3};
  NablaVariables{3}(:,j) = dummy(:);
end

return