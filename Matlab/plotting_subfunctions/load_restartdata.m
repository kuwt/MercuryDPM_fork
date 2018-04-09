function data=load_restartdata(name)

global rawdata;



%% First load the data all the data in rawdata
rawdata=textread(name,'%s');







%% Now deal with particles

% First get the number of particeles
[data.Nparticles,index]=extract_info('Particles',1,'double');


% Now load all the data of the particles
data.Particles=zeros(data.Nparticles,15);

for i=1:data.Nparticles
    data.Particles(i,:)=str2double(rawdata(index-13+i*15:index+i*15+1));
end


%% Domain

data.xmax=extract_info('xmax',1,'double');
data.xmin=extract_info('xmin',1,'double');

data.ymax=extract_info('ymax',1,'double');
data.ymin=extract_info('ymin',1,'double');

data.zmax=extract_info('zmax',1,'double');
data.zmin=extract_info('zmin',1,'double');

%% Walls

%First extract how many walls (total)
data.Nwalls=extract_info('Walls',1,'double');

%Now get the number of finite walls for each wall.
[Nfinite,index]=extract_info('numFiniteWalls',1,'double');

%Set the number of finite and infinite walls both to zero
data.NIwalls=0;
data.NFwalls=0;
    

%Loop over each wall reading in the information for it.
for i=1:data.Nwalls
    
 
    if Nfinite(i)==0
        
        %In this case the wall is infinte so have a position and a single
        %normal
        
        
        data.NIwalls=data.NIwalls+1;   
        data.Iwalls(data.NIwalls).position=str2double(rawdata(index(i)+7));
        data.Iwalls(data.NIwalls).normal(1:3)=str2double(rawdata(index(i)+(3:5)));
        
    else
        
        % In this case the wall is finite and considers of a number of
        % dim-1 planes
        
        data.NFwalls=data.NFwalls+1;
        data.Fwalls(data.NFwalls).num_planes=Nfinite(i);
            
        
        %Now loop over each plane makes up the Wall
        for j=1:Nfinite(i)
                        
            data.Fwalls(data.NFwalls).planes(j).position=str2double(rawdata(index(i)+1+j*6));
            data.Fwalls(data.NFwalls).planes(j).normal(1:3)=str2double(rawdata(index(i)+(3:5)+j*6-6));
            
        end% loop over planes of the finite wall
        
    end% loop over all the walls

end% procedure load data




%% Finally add the rawdata in the data object in case there is extra data
%% the user might want later

data.rawdata=rawdata;




end

function [info,index]=extract_info(name,no_of_lines,type)

global rawdata

index=find(strcmp(name,rawdata));

% This can be replace with a line like the one which extract particles, but
% I do not get this syntex ask Thomas, as he wrote the orginally line.
for i=1:length(index)
    info(i)=rawdata(index(i)+1:index(i)+no_of_lines);
end

if strcmp(type,'double')
    info=str2double(info);
end

end