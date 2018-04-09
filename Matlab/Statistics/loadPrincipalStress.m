function data=loadPrincipalStress(data)

if iscell(data)
  % this is used if data is a cell array (of structs)
  for i=1:length(data)
    data{i}=loadPrincipalStress(data{i});
  end
else
  %here, data is a singel data ste (a struct)
  
  data.PrincipalEVec=ones(3,3,length(data.x))*nan;
  data.PrincipalEVal=ones(3,length(data.x))*nan;
  data.PrincipalAngle=ones(3,length(data.x))*nan;
  data.PrincipalEVecKin=ones(3,3,length(data.x))*nan;
  data.PrincipalEValKin=ones(3,length(data.x))*nan;
  data.PrincipalAngleKin=ones(3,length(data.x))*nan;
  data.PrincipalEVecNor=ones(3,3,length(data.x))*nan;
  data.PrincipalEValNor=ones(3,length(data.x))*nan;
  data.PrincipalAngleNor=ones(3,length(data.x))*nan;
  data.PrincipalEVecTan=ones(3,3,length(data.x))*nan;
  data.PrincipalEValTan=ones(3,length(data.x))*nan;
  data.PrincipalAngleTan=ones(3,length(data.x))*nan;
  data.FabricPrincipalEVec=ones(3,3,length(data.x))*nan;
  data.FabricPrincipalEVal=ones(3,length(data.x))*nan;
  data.FabricPrincipalAngle=ones(3,length(data.x))*nan;

  %go through all points
  for j=1:length(data.x)
    %put Stress in a matrix
    Stress = [ ...
      data.StressXX(j) data.StressXY(j) data.StressXZ(j); ...
      data.StressXY(j) data.StressYY(j) data.StressYZ(j); ...
      data.StressXZ(j) data.StressYZ(j) data.StressZZ(j) ];
    %check if stress is non-zero and finite, i.e. can be split into ev's
    if ~(Stress==zeros(3)| isnan(max(Stress(:))))
      [v1,d1]=eig(Stress);
      [v,d]=dsort_ev(v1,d1);
      d = diag(d);
      p=sum(d)/3;
      data.PrincipalEVec(:,:,j)=v;
      data.PrincipalEVal(:,j)=d';
      %angle with z-axis
      %P=acosd([0 0 1]*v); P(P>90)=180-P(P>90); 
      P=atand(v(1,:)./v(3,:)); 
      data.PrincipalAngle(:,j)=P;
    end
    %put KineticStress in a matrix
    KineticStress = [ ...
      data.KineticStressXX(j) data.KineticStressXY(j) data.KineticStressXZ(j); ...
      data.KineticStressXY(j) data.KineticStressYY(j) data.KineticStressYZ(j); ...
      data.KineticStressXZ(j) data.KineticStressYZ(j) data.KineticStressZZ(j) ];
    %check if KineticStress is non-zero and finite, i.e. can be split into ev's
    if ~(KineticStress==zeros(3)| isnan(max(KineticStress(:))))
      [v1,d1]=eig(KineticStress);
      [v,d]=dsort_ev(v1,d1);
      d = diag(d);
      p=sum(d)/3;
      data.PrincipalEVecKin(:,:,j)=v;
      data.PrincipalEValKin(:,j)=d';
      %angle with z-axis
      %P=acosd([0 0 1]*v); P(P>90)=180-P(P>90); 
      P=atand(v(1,:)./v(3,:)); 
      data.PrincipalAngleKin(:,j)=P;
    end
    %put NormalStress in a matrix
    NormalStress = [ ...
      data.NormalStressXX(j) data.NormalStressXY(j) data.NormalStressXZ(j); ...
      data.NormalStressXY(j) data.NormalStressYY(j) data.NormalStressYZ(j); ...
      data.NormalStressXZ(j) data.NormalStressYZ(j) data.NormalStressZZ(j) ];
    %check if NormalStress is non-zero and finite, i.e. can be split into ev's
    if ~(NormalStress==zeros(3)| isnan(max(NormalStress(:))))
      [v1,d1]=eig(NormalStress);
      [v,d]=dsort_ev(v1,d1);
      d = diag(d);
      p=sum(d)/3;
      data.PrincipalEVecNor(:,:,j)=v;
      data.PrincipalEValNor(:,j)=d';
      %angle with z-axis
      %P=acosd([0 0 1]*v); P(P>90)=180-P(P>90); 
      P=atand(v(1,:)./v(3,:)); 
      data.PrincipalAngleNor(:,j)=P;
    end
    %put TangentialStress in a matrix
    TangentialStress = [ ...
      data.TangentialStressXX(j) data.TangentialStressXY(j) data.TangentialStressXZ(j); ...
      data.TangentialStressXY(j) data.TangentialStressYY(j) data.TangentialStressYZ(j); ...
      data.TangentialStressXZ(j) data.TangentialStressYZ(j) data.TangentialStressZZ(j) ];
    %check if TangentialStress is non-zero and finite, i.e. can be split into ev's
    if ~(TangentialStress==zeros(3)| isnan(max(TangentialStress(:))))
      [v1,d1]=eig(TangentialStress);
      [v,d]=dsort_ev(v1,d1);
      d = diag(d);
      p=sum(d)/3;
      data.PrincipalEVecTan(:,:,j)=v;
      data.PrincipalEValTan(:,j)=d';
      %angle with z-axis
      %P=acosd([0 0 1]*v); P(P>90)=180-P(P>90); 
      P=atand(v(1,:)./v(3,:)); 
      data.PrincipalAngleTan(:,j)=P;
    end
    %put Fabric in a matrix
    Fabric = [ ...
      data.FabricXX(j) data.FabricXY(j) data.FabricXZ(j); ...
      data.FabricXY(j) data.FabricYY(j) data.FabricYZ(j); ...
      data.FabricXZ(j) data.FabricYZ(j) data.FabricZZ(j) ];
    %check if Fabric is non-zero and finite, i.e. can be split into ev's
    if ~(Fabric==zeros(3)| isnan(max(Fabric(:))))
      [v1,d1]=eig(Fabric);
      [v,d]=dsort_ev(v1,d1);
      d = diag(d);
      p=sum(d)/3;
      data.FabricPrincipalEVec(:,:,j)=v;
      data.FabricPrincipalEVal(:,j)=d';
      %angle with z-axis
      %P=acosd([0 0 1]*v); P(P>90)=180-P(P>90); 
      P=atand(v(1,:)./v(3,:)); 
      data.FabricPrincipalAngle(:,j)=P;
    end
  end

  %now do the same with the flow value
  Stress=data.FlowStress;
  if ~(Stress==zeros(3)| isnan(max(Stress(:))))
    [v1,d1]=eig(Stress);
    [v,d]=dsort_ev(v1,d1);
    d = diag(d);
    p=sum(d)/3;
    data.FlowPrincipalEVec=v;
    data.FlowPrincipalEVal=d';
    %angle with z-axis
    %P=acosd([0 0 1]*v); P(P>90)=180-P(P>90); 
    P=atand(v(1,:)./v(3,:)); 
    data.FlowPrincipalAngle=P;
  end
end

return