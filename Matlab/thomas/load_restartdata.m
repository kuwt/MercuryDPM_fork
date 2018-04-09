function data=load_restartdata(name)

global rawdata;

%% First load the data all the data in rawdata
rawdata=textread(name,'%s');

%now read in values into a structured data
L=length(rawdata)
data=struct();
i=1;
while (i<=L)
    j=i+1;
    if strcmp(rawdata{i},'name') %exception: name is a string
        data.(rawdata{i})=rawdata{j};
        j=j+1;
    else %standard format in restart: first variable name, then array of values (or no value)
        if strcmp(rawdata{i},'TSP')
            j=j+1; %skip the first value (BP)
        end
        %read in the length and values of the variable
        values = []; 
        while (j<=L) 
            val = str2double(rawdata{j});
            if (isfinite(val))
                values(end+1) = val;
            else
                break;
            end
            j=j+1;
        end
        if strcmp(rawdata{i},'TSP')
            if ~isfield(data,'BP'), data.BP=[]; end
            data.TSPSpring{size(data.BP+1,1)+1}=values(18:end);
            values = values(1:17);
            i=i+1; %name BP
        end
        %write into data member variable (create matrix if the name is repeated)
        if isfield(data,rawdata{i})
            data.(rawdata{i})(end+1,:)=values;
        else
            data.(rawdata{i})=values;
        end
    end
    i=j;
end

%handle exceptions
if isfield(data,'AngledPeriodicBoundary')
    data.AngledPeriodicBoundary=struct();
    data.AngledPeriodicBoundary.normal_left=data.normal_left;
    data.AngledPeriodicBoundary.normal_right=data.normal_right;
    data.AngledPeriodicBoundary.origin=data.origin;
    rmfield(data,{'normal_left','normal_right','origin'})
end
return
