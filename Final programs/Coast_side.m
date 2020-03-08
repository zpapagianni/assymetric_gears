function [xtooth2,ytooth2,xtooth2k,ytooth2k,p,k] = Coast_side(Z,Cs,Cc,mn,Ck,Cf,xdrive,ydrive,b0,tos);

%% dialog box

prompt = {'Profile angle for coast side ',};
dlg_title = 'COAST SIDE';
num_lines = 1;
def = {'25'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
aon= str2num(answer{1});


% aon=25; 
aon=degtorad(aon);       
Cm=0;
x=Cm;
%% metatroph ms se mn
ms=mn/cos(b0);
aos1=atan(tan(aon)/cos(b0));

%% Ypologismos aktinwn               
ron=Z*ms/2;
rg=ron*cos(aos1);        % aktina vasikou kuklou
rk=ron+(Ck+x)*mn;            % aktina kefalis
rf=ron-(Cf-x)*mn;            % aktina podos
rc=Cc*mn;
hf=Cf*mn;    
sos1=Cs*pi*ms+2*ms*Cm*tan(aon)
l0=tos-sos1; 

z1=(2*(Cf-x-Cc*(1-sin(aon))))/sin(aon)^2;
fprintf('Minimum number of Teeth for coast side is: %f\n\n',z1);

%% Ypologsismos 
ph1=((pi/2)*Cs*mn+(hf+rc*((1-sin(aon))/sin(aon)))*tan(aon))/((Z/2)*mn); 
phif=((2*cos(b0)^2)*(hf-x*mn+rc*(1+sin(aon)*tan(b0)^2)))/(Z*mn*tan(aon));

S=sqrt(cos(b0)^2+tan(aon)^2);
rifo=sqrt((Z*mn/(2*S))^2+((tan(aon)/cos(b0))*(Z*mn/(2*S))-(S/tan(aon))*(hf-x*mn-rc*(1-sin(aon))))^2); 

%fillet
j=1;
for ph=0:phif/1000:phif 
    for ln=pi/100000:(pi-2*pi/100000)/1000:pi-pi/100000
        phi=((2*(cos(b0))^2)/(Z*mn*tan(ln)))*(hf-x*mn-rc*(1+sin(ln)*tan(b0)^2)); 
        if abs(ph-phi)<0.0001                      
            l=atan(tan(ln)/cos(b0)); 
            ksfo(j)=Z*mn*sin(ph1+phi)/(2*cos(b0))-((hf-x*mn-rc*(1-sin(ln)))*cos(ph1+phi-l))/sin(l); 
            hfo(j)=Z*mn*cos(ph1+phi)/(2*cos(b0))+((hf-x*mn-rc*(1-sin(ln)))*sin(ph1+phi-l))/sin(l); 
            j=j+1;
        end
    end
end

%ekseligmeni
j=1;
clear ph
for r=(Z*mn*cos(aos1)/(2*cos(b0))):(1.1*rk-(Z*mn*cos(aos1)/(2*cos(b0))))/100000:1.1*rk            
    ph=((2*cos(b0)/(Z*mn))*((Z*mn)/(2*cos(b0))-hf+x*mn-rc*((1-sin(aon))/sin(aon)))-(sqrt(((r*2*cos(b0)/(Z*mn))^2)-cos(aos1)^2)/sin(aos1)))*tan(aos1);
    ksio(j)=(Z./2).*(mn./cos(b0)).*sin(ph1+ph)-((hf-x.*mn+rc.*((1-sin(aon))./sin(aon))*tan(aos1)+(Z.*mn.*ph./(2.*cos(b0))))).*cos(aos1).*cos(ph1+ph-aos1);
    hio(j)=(Z./2).*(mn./cos(b0)).*cos(ph1+ph)+((hf-x.*mn+rc.*((1-sin(aon))./sin(aon))*tan(aos1)+(Z.*mn.*ph./(2.*cos(b0))))).*cos(aos1).*sin(ph1+ph-aos1);
    j=j+1;
end

drf=sqrt(ksfo.^2+hfo.^2);
dri=sqrt(ksio.^2+hio.^2);



%% kataxwrhsh timwn tou fillet 
% j=length(drf);
% p=1;
% if (Z>z1)
%   for i=1:j
%       if (drf(i)>=rf) && (drf(i)<=rifo)
%           xfillet(p)=ksfo(i);
%           yfillet(p)=hfo(i);
%           p=p+1;
%       end
%    end
% elseif (Z<=z1) 
%     for i=1:j
%         if (drf(i)>=rf) && (drf(i)<=rifo)
%             xfillet(p)=-ksfo(i);
%             yfillet(p)=-hfo(i);
%             p=p+1;
%         end
%      end
% end
 
j=length(drf);
p=1;
for i=1:j
      if (drf(i)>=rf) && (drf(i)<=rifo)
          xfillet(p)=ksfo(i);
          yfillet(p)=hfo(i);
          p=p+1;
      end
end
%% kataxwrhsh timwn ths ekseligmenhs
j=length(dri);
p=1;
for i=1:j
    if (dri(i)>=rifo) && (dri(i)<=rk)
        xinv(p)=ksio(i);
        yinv(p)=hio(i);
        p=p+1;
    end
end

%% peristrofh fillet-involute wste na tairiaksoun ston rs

ainv=atan(yinv(1)/xinv(1));
y=length(xfillet);
afillet=atan(yfillet(y)/xfillet(y));
rot=(abs(ainv)-afillet)
if x<0
    rot=rot+pi;
end

xtooth=[xfillet.*cos(rot)-yfillet.*sin(rot) xinv];   % peristrofh tou pinaka fillet
ytooth=[xfillet.*sin(rot)+yfillet.*cos(rot) yinv]; 




%% find ron-involute intersection
d=sqrt(xtooth.^2+ytooth.^2);
dmax=10000;
for i=1:length(d)
    if abs(d(i)-ron)<dmax
        dmax=abs(d(i)-ron);
        p=i;
    end
end
p2=p;
th3=atan(ytooth(p)/xtooth(p));

%% metafora shmeiwn symmetrika ston aksona syntetagmenwn
rot_angle=abs(-pi/2+abs(th3)+(sos1/ron)/2);
xtooth1= xtooth.*cos(rot_angle)+ytooth.*sin(rot_angle);
ytooth1= -xtooth.*sin(rot_angle)+ytooth.*cos(rot_angle);


%% ypologismos gonias kai metasxhmatismos wste coast side na metaferthei sto shmeio tomhs
xcoast=xtooth1(p2);
ycoast=ytooth1(p2);
drive=atan(ydrive/xdrive)
coast=atan(ycoast/xcoast)
k=(atan(ydrive/xdrive))-(atan(ycoast/xcoast));

%% mirror tooth

[xtooth2,ytooth2,xtooth2k,ytooth2k] = mirror_tooth(xtooth1,ytooth1,ron,sos1,k);

%% plot dontiou
% plot(xtooth1,ytooth1,'b')
% plot(xtooth2,ytooth2,'r') 
plot(xtooth2k,ytooth2k,'k') 
% %% plot multiple times
% for th=0:2*pi/Z:2*pi
%     x1k=(xtooth2k).*cos(th)-(ytooth2k).*sin(th);
%     y1k=(xtooth2k).*sin(th)+(ytooth2k).*cos(th);
%     plot(x1k,y1k,'b');
% 
% end

% 

end

