% ********************Asymmetric helical Gears ******************

clear all
clc
hold on






%% dialog box

prompt = {'Profile angle','Helix angle','Number of tooth'};
dlg_title = 'DRIVE SIDE';
num_lines = 3;
def = {'14.5','45','20'};
answer = inputdlg(prompt,dlg_title,num_lines,def)
aon= str2num(answer{1}); % einai i hmigwnia odontws tis kathetis katatomis
b0=str2num(answer{2});% einai i klisi twn odontwn ws pros aksona toy troxou ston arxiko kulindro
Z=str2num(answer{3}); % arithmos dontiwn !!!!!!!!!!!!!!!!

fprintf('');
aon,b0,Z





%% arxikes synthikes
% aon=14.5;                % einai i hmigwnia odontws tis kathetis katatomis
% b0=20;                   % einai i klisi twn odontwn ws pros aksona toy troxou ston arxiko kulindro
% mn=4;                    % einai to module tis ka8eths tomis ston arxiko kulindro
% Z=20;                    % arithmos dontiwn !!!!!!!!!!!!!!!!

Cs=0.5;                    % pairnei times [0.475,0.5]-syntelestis paxous
mn=20;                   
Cc=0.25;                 % sintelestis gia tin aktina kampulotitas
Ck=1;
Cf=1.25;                 %synhthws [1,25-1,30]
Cm=0;
x=Cm


%% metatroph ms se mn
aon=degtorad(aon);       % metatropi tis gwnias a0n apo moires se rad
b0=degtorad(b0);         % metatropi tis gwnias b0 apo moires se radms=mn/cos(b0);                  % module ths metopikhs tomhs

ms=mn/cos(b0);   
aos=atan(tan(aon)/cos(b0));

sos=Cs*pi*ms+2*ms*Cm*tan(aon);  % paxos tou kanonikou odontos ston arxiko kyklo tis metwpikis tomis
tos=pi*ms;                      % vima ston arxiko kyklo
ton=tos*cos(b0);                
fprintf('Tooth thickness = %f\n',sos);
fprintf('Width of space tos = %f\n',tos)


l0=tos-sos;                     %diakeno ston arxiko kyklo

fprintf('Backlash l0 = %f\n',l0)
%% elegxos timwn Cf,Ck,Cs
% Cf
cfmax=((1-Cs)*pi)/(2*tan(aon));
if Cf>cfmax
    error('please change cf value\nMax cf value is: %f\n',cfmax);
end
%Cc
ccmax=((1-Cs)*pi/2-Cf*tan(aon))/(cos(aon)-(1-sin(aon))*tan(aon));
if Cc>ccmax
    error('please change cc value\nMax cc value is: %f\n',ccmax);
end

%% ******* upologismos akmax kai Ck *********

invakmax=tan(aon)-aon+(pi/Z)*Cs;
dth=pi/100;
akmax=0;
j=1;
fi=0;
while j>0.0001
    if fi>invakmax
        akmax=akmax-dth;
        dth=dth/10;
    end
    
    akmax=akmax+dth;
    fi=tan(akmax)-akmax;
    j=abs(invakmax-fi);
end
ckmax=Z/2*(((cos(aon))/(cos(akmax)))-1);
if Ck<0
    error('please reduce number of teeth\n\n');
end
if Ck>ckmax
    error('please change ck value\nMax ck value is: %f\n',ckmax);
end

%ypologismos elexistou arithmoy odontwn
zmin=(2*(Cf-x-Cc*(1-sin(aon))))/sin(aon)^2;
fprintf('Minimum number of Teeth for drive side is: %f\n\n',zmin);


%% ypologismos akitinon
               
ron=Z*ms/2;
rg=ron*cos(aos);             % aktina vasikou kuklou
rk=ron+(Ck+x)*mn;            % aktina kyklou kefalis
rf=ron-(Cf-x)*mn;            % aktina kyklou podos
rc=Cc*mn;
hf=Cf*mn;               %ypsos kefalis tou odontos tou kanona
fprintf('Pitch circle ron = %f\n',ron)
%% ypologismos shmeion

ph1=((pi/2)*Cs*mn+(hf+rc*((1-sin(aon))/sin(aon)))*tan(aon))/((Z/2)*mn); %i gwnia pou sximatizei mia eutheia Op1 me ton aksona ?
phif=((2*cos(b0)^2)*(hf-x*mn+rc*(1+sin(aon)*tan(b0)^2)))/(Z*mn*tan(aon));% dinei tis sintetagmenes tou simeiou A, opou sinantatai i ekseiligmeni me to fillet

S=sqrt(cos(b0)^2+tan(aon)^2);
rifo=sqrt((Z*mn/(2*S))^2+((tan(aon)/cos(b0))*(Z*mn/(2*S))-(S/tan(aon))*(hf-x*mn-rc*(1-sin(aon))))^2); %otan dn yparxoun ypokopes ston odonta tou koinou simeiou A
%fillet
j=1;
for ph=0:phif/1000:phif %perigrafi tou fillet me xrhsh tis gwnias ? metavallomenis metaksi, otan dn yparxoun ypokopes
    for ln=pi/100000:(pi-2*pi/100000)/1000:pi-pi/100000
        phi=((2*(cos(b0))^2)/(Z*mn*tan(ln)))*(hf-x*mn-rc*(1+sin(ln)*tan(b0)^2)); %sxesi 19 apo to vivlio kwstopoulou
        if abs(ph-phi)<0.0001                       % +
                 %sth thesi Fn toy fillet i katheti epi ti katatomh sxhmatizei ti gwnia ?n (sti katheti tomh)
            l=atan(tan(ln)/cos(b0)); %sxesi tis gwnias ? metaksi metwpikis kai kathetis tomis
            a1=Z*mn*sin(ph1+phi)/(2*cos(b0));
            a2=((hf-x*mn-rc*(1-sin(ln)))*cos(ph1+phi-l))/sin(l);
            b1=Z*mn*cos(ph1+phi)/(2*cos(b0));
            b2=((hf-x*mn-rc*(1-sin(ln)))*sin(ph1+phi-l))/sin(l);
            ksfo(j)=a1-a2; %einai i tetmimeni sto simeio F tou fillet
            hfo(j)=b1+b2;  %einai i tetagmimeni sto simeio F tou fillet
            j=j+1;
        end
    end
end
%  plot(ksfo,hfo,'r')
 %ekseleigmeni
j=1;
clear ph
for r=(Z*mn*cos(aos)/(2*cos(b0))):(1.1*rk-(Z*mn*cos(aos)/(2*cos(b0))))/100000:1.1*rk            % 1.1*rk giati den ftanei to rk
    ph=((2*cos(b0)/(Z*mn))*((Z*mn)/(2*cos(b0))-hf+x*mn-rc*((1-sin(aon))/sin(aon)))-(sqrt(((r*2*cos(b0)/(Z*mn))^2)-cos(aos)^2)/sin(aos)))*tan(aos);

    ksio(j)=(Z./2).*(mn./cos(b0)).*sin(ph1+ph)-((hf-x.*mn+rc.*((1-sin(aon))./sin(aon))*tan(aos)+(Z.*mn.*ph./(2.*cos(b0))))).*cos(aos).*cos(ph1+ph-aos);
    hio(j)=(Z./2).*(mn./cos(b0)).*cos(ph1+ph)+((hf-x.*mn+rc.*((1-sin(aon))./sin(aon))*tan(aos)+(Z.*mn.*ph./(2.*cos(b0))))).*cos(aos).*sin(ph1+ph-aos);
    j=j+1;
end
%  plot(ksio,hio,'b')
drf=sqrt(ksfo.^2+hfo.^2);   %apostasi apo to (0,0)
dri=sqrt(ksio.^2+hio.^2);

%%  kataxwrhsh timwn tou fillet 
j=length(drf);
p=1;
for i=1:j
    if (drf(i)>=rf) && (drf(i)<=rifo)
          xfillet(p)=ksfo(i);
          yfillet(p)=hfo(i);
          p=p+1;
    end
end
% plot(xfillet,yfillet,'r')
%% kataxwrhsh timwn ths ekseligmenhs
j=length(dri);
p=1;
for i=1:j
    if (dri(i)>=rifo) && (dri(i)<rk)
        xinv(p)=ksio(i);
        yinv(p)=hio(i);
        p=p+1;
    end
end

% plot(xinv,yinv,'g')
%% peristrofh fillet-involute wste na tairiaksoun ston rs

ainv=atan(yinv(1)/xinv(1));
y=length(xfillet);
afillet=atan(yfillet(y)/xfillet(y));
% plot(xfillet(y),yfillet(y),'g')
% if ainv<0 
%    rot=(-ainv-afillet)
% else
%    rot=(-afillet+ainv);
% end
rot=(abs(ainv)-afillet)
if x<0
    rot=rot+pi;
end
xtooth=[xfillet.*cos(rot)-yfillet.*sin(rot) xinv];   % peristrofh tou pinaka fillet
ytooth=[xfillet.*sin(rot)+yfillet.*cos(rot) yinv]; 
% plot(xtooth, ytooth)

%% plot kyklon
 
vhma=2*pi/360;
a=-pi:vhma:pi;
%rk
x=rk*cos(a);
y=rk*sin(a);
plot(x,y,'k--');

%rf
x=rf*cos(a);
y=rf*sin(a);
plot(x,y,'k');


%ron
x=ron*cos(a);
y=ron*sin(a);
plot(x,y,'k-.');

%rg
x=rg*cos(a);
y=rg*sin(a);
plot(x,y,'r-.');

legend ('rk','rf', 'ron', 'rg')
%% axes
xL = xlim;
yL = ylim;
line([0 0], yL);  %x-axis
line(xL, [0 0]);  %y-axis
% legend('rk','rf','ron','rg','final tooth','coast side','coast side xwris metasxhmatismo','Location','southwest')


%% find ron-involute intersection
% 
d=sqrt(xtooth.^2+ytooth.^2);
dmax=10000;
for i=1:length(d)
    if abs(d(i)-ron)<dmax
        dmax=abs(d(i)-ron);
        p=i;
    end
end
th=atan(ytooth(p)/xtooth(p));
p1=p;
 

%% metafora shmeiwn symmetrika ston aksona syntetagmenwn
rot_angle=abs(-pi/2+abs(th)+((sos/ron)/2));    %(sos/ron)/2=gonia
xtooth1= xtooth.*cos(rot_angle)+ytooth.*sin(rot_angle);
ytooth1= -xtooth.*sin(rot_angle)+ytooth.*cos(rot_angle); 

xdrive=xtooth1(p1);
ydrive=ytooth1(p1);

fprintf('Point on pitch circle')
xdrive
ydrive



% % plot dontiou
plot(xtooth1,ytooth1,'k')
% % % plot multiple times
% % % 
% % % for th=0:2*pi/Z:2*pi
% % %     x1k=(xtooth1).*cos(th)-(ytooth1).*sin(th);
% % %     y1k=(xtooth1).*sin(th)+(ytooth1).*cos(th);
% % %     plot(x1k,y1k,'b');
% % % 
% % % end


%% ypologismos coast side
[xtooth2,ytooth2,xtooth2k,ytooth2k,p2,k] = Coast_side(Z,Cs,Cc,mn,Ck,Cf,xdrive,ydrive,b0,tos);

% %% find rk-involute intersection
% 
% w1=atan(ytooth1(length(xtooth1))/xtooth1(length(xtooth1)));
% w2=atan(ytooth2k(length(xtooth2k)))/xtooth2k(length(xtooth2k));
% radtodeg(w1)
% radtodeg(w2)
% %%plot toxou
% 
% vhma=2*pi/360;
% 
% a2=90-w2
% a=w1:vhma:a2;
% x=rk*cos(a);
% y=rk*sin(a);
% plot(x,y);

%% ypologismos paxous kai gonias

fprintf('Rotation angle k= %f\n',k)
dist=abs(xtooth2k(p2)-xtooth1(p1));           %apostasi shmeiwn arxikou kyklou 
fprintf('Tooth thickness with rotation dist = %f\n',dist)
dist1=abs(xtooth2(p2)-xtooth1(p1));           
fprintf('Tooth thickness without rotation dist1 = %f\n',dist1)
% plot(xtooth2,ytooth2,'g')

angledrive=atan(ytooth1(p1)/xtooth1(p1));
anglecoast=atan(ytooth2k(p2)/xtooth2k(p2));
angle=pi-abs(anglecoast)-abs(angledrive); 
angleInDegrees = radtodeg(angle);
fprintf('Spots angle on drive side angledrive =%f\n',angledrive)
fprintf('Spots angle on coast side anglecoast = %f\n',anglecoast)
fprintf('Real thickness angle = %f\n',angle)
fprintf('Real thickness angle = %f\n',angleInDegrees)
angletrue=sos/ron;
fprintf('Theoretical thickness angle angletrue = %f\n',angletrue)

 
% %% axes
% xL = xlim;
% yL = ylim;
% line([0 0], yL);  %x-axis
% line(xL, [0 0]);  %y-axis
% % legend('rk','rf','ron','rg','final tooth','coast side','coast side xwris metasxhmatismo','Location','southwest')

%% ******************** eksagwgh se txt ***********************
fileID=fopen('helical_gear.txt','w'); % ??????? ??????? ??? ???????
size=length(xtooth2)+length(xtooth1);
xm=length(xtooth1);
diaf1=ytooth1(xm)-rk
for i=1:5:xm  
    fprintf(fileID,'%f \t %f \t %f \t \r\n',[xtooth1(i) ytooth1(i) 0]);
%     fprintf(fileID,'%f \t %f \t %f \t \r\n',[xtooth2(i) ytooth2(i) 0]);
end
fclose(fileID);

fileID=fopen('helical_gear_coast.txt','w'); % ??????? ??????? ??? ???????
size=length(xtooth2k)+length(xtooth2k);
xm=length(xtooth2k);
diaf2=ytooth2k(xm)-rk
for i=1:5:xm  
%     fprintf(fileID,'%f \t %f \t %f \t \r\n',[xtooth(i) ytooth(i) 0]);
    fprintf(fileID,'%f \t %f \t %f \t \r\n',[xtooth2k(i) ytooth2k(i) 0]);
end
fclose(fileID);
