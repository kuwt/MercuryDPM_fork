function data=loadPrincipalStress(data)

if isstruct(data)
  for i=1:length(data)
    data{i}=loadPrincipalStress(data{i});
  end
else
  for j=1:length(data.x)
    Stress = [ ...
      data.StressXX(j) data.StressXY(j) data.StressXZ(j); ...
      data.StressXY(j) data.StressYY(j) data.StressYZ(j); ...
      data.StressXZ(j) data.StressYZ(j) data.StressZZ(j) ];
    if (Stress==zeros(3)| isnan(max(Stress(:))))

    else
      [v1,d1]=eig(Stress);
      [v,d]=dsort_ev(v1,d1);
      d = diag(d);
      p=sum(d)/3;
      data.PrincipalEVec(:,:,j)=v;
      data.PrincipalEVal(j,:)=d;
      %angle with z-axis
      %data.PrincipalAngle(j,:)=atan(v(3,:)./v(1,:))/pi*180;
      P=acosd([0 0 1]*v);
      P(P>90)=180-P(P>90)
      data.PrincipalAngle(j,:)=P;
    end
  end
end

return