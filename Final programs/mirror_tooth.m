function [xtooth2,ytooth2,xtooth2k,ytooth2k] = mirror_tooth(xtooth1,ytooth1,ron,sos,k)

d=sqrt(xtooth1.^2+ytooth1.^2);
dmax=10000;
for i=1:length(d)
    if abs(d(i)-ron)<dmax
        dmax=abs(d(i)-ron);
        p=i;
    end
end
th=atan(ytooth1(p)/xtooth1(p));

%% syntetagmenes shmeiwn xwris metafora

th2=sos/ron-2*(atan(ytooth1./xtooth1)-th);                     
xtooth2=xtooth1.*cos(th2)-ytooth1.*sin(th2);
ytooth2=xtooth1.*sin(th2)+ytooth1.*cos(th2);


%% syntetagmenes shmeiwn katatomhs me metafora ws pros ton arxiko kyklo
 
xtooth2k=xtooth1.*cos(th2+k)-ytooth1.*sin(th2+k);
ytooth2k=xtooth1.*sin(th2+k)+ytooth1.*cos(th2+k);


end
