cochain=spinw;
beta=105.579;
%beta=90;
db=0.180125;
%db=0.1;
cochain.genlattice('lat_const',[5.4701 9.297 10.256],'angled',[90 beta 90 ],'spgr',0);
cochain.addatom('r',[0 0 0],'S',1/2,'label','MCo2','color','blue');
cochain.addatom('r',[1/2 db 0],'S',1/2,'label','MCo2','color','red');
cochain.addatom('r',[0 1/2 1/2],'S',1/2,'label','MCo2','color','blue');
cochain.addatom('r',[1/2 db+1/2 1/2],'S',1/2,'label','MCo2','color','red');
disp('Magnetic lattice:')
cochain.table('matom')
plot(cochain,'range',[5 2 2])
cochain.gencoupling('maxDistance',7);
cochain.table('bond',[])

%Jxy=-0.640281;
%Delta=0.0;
%Jpp=0.154127;
%Jzp=0.114217;
%Jab=0.921438;

%Jxy=-0.653215;
%Delta=0.5;
%Jpp=-0.173512;
%Jzp=0.117498;
%Jab=0.733949;

%Jxy=-0.648642;
%Delta=1.0;
%Jpp=-0.495152;
%Jzp=0.116331;
%Jab=0.581395;

%Jxy=-0.629213;
%Delta=1.5;
%Jpp=-0.789017;
%Jzp=0.111458;
%Jab=0.47085;

%Jxy=-0.602486;
%Delta=2;
%Jpp=-1.0493;
%Jzp=0.104974;
%Jab=0.394348;

%Delta=1.6;%fixed without Jab2
%Jxy=-0.639462 - 0.0502642*Delta + 0.0482202*Delta^2 - 0.00684786*Delta^3;
%Jpp=0.155883 - 0.660586*Delta - 0.00872262*Delta^2 + 0.0192051*Delta^3;
%Jzp=0.113987 + 0.0129309*Delta - 0.0125754*Delta^2 + 0.00190864*Delta^3;
%Jab=.923463 - 0.419368*Delta + 0.0783831*Delta^2 - 0.000204204*Delta^3;
%Jab=0.0;
%Jab2=0.2;

%Delta=2.0;%fixed e2=2.7,e1=1.9
%Jxy=-1.19484 + 0.604703*Delta - 0.251514*Delta^2 + 0.0476954*Delta^3;
%Jpp=-0.297308 - 0.456785*Delta + 0.0692861*Delta^2 - 0.00678025*Delta^3;
%Jzp=0.348908 - 0.294853*Delta + 0.133483*Delta^2 - 0.0259955*Delta^3;
%Jab=0.103413;
%Jab2=0.672124;

Delta=1.95;%fixed emin=1,e2=1.9
Jxy=-202.015 - 139.476/Delta + 304.675/Delta^0.5 + 37.6273*Delta - 5.64084*Delta^2 + 0.377229*Delta^3;
Jpp=-0.882279;
Jzp=0.0648605;
Jab=0.103413;
Jab2=0.672124;

% Jxy=-0.586085;
% Delta=1.70954;
% Jpp=-0.67518;
% Jzp=0.0609267;
% Jab=0.1;
% Jab2=1.0;

% Jxy=-0.565325;
% Delta=0.78607;
% Jpp=-0.15596;
% Jzp=0.0650547;
% Jab=0.5;
% Jab2=1.6-2*Jab;

% Jxy=-0.564428;
% Delta=0.971207;
% Jpp=-0.260012;
% Jzp=0.0649052;
% Jab=0.6;
% Jab2=1.5-2*Jab;

Jxy=-0.56195;
Delta=1.17218;
Jpp=-0.371267;
Jzp=0.0644937;
Jab=0.4;
Jab2=1.4-2*Jab;

% Jxy=-0.557565;
% Delta=1.39453;
% Jpp=-0.49141;
% Jzp=0.06377;
% Jab=0.6;
% Jab2=1.3-2*Jab;

% Jxy=-0.550872;
% Delta=1.64611;
% Jpp=-0.622706;
% Jzp=0.0626773;
% Jab=0.6;
% Jab2=1.2-2*Jab;

Jxy=-0.826666;
Delta=1.0;
Jpp=-0.18666666;
Jzp=0.527973;
Jab=0.74;
Jab2=1.06;

%R=[0 cos((180-105.579)/180*pi) 1*cos(15.579/180*pi);0 -1 0; 1 0 0];
%R=[0 0 1;0 -1 0; 1 0 0];
%spinangle=0.115377;
%Rmat=[-sin(spinangle) 0 cos(spinangle);0 -1 0; cos(spinangle) 0 sin(spinangle)];
Jxxmat=[Jxy+Jpp -sqrt(3)*Jpp Jzp/2;-sqrt(3)*Jpp Jxy-Jpp -Jzp*sqrt(3)/2;Jzp/2 -Jzp*sqrt(3)/2 Delta*Jxy];
Jyymat=[Jxy+Jpp sqrt(3)*Jpp Jzp/2;sqrt(3)*Jpp Jxy-Jpp Jzp*sqrt(3)/2;Jzp/2 Jzp*sqrt(3)/2 Delta*Jxy];
jxx1=-0.440615;
jyy1=-0.409365;
%jzz1=-0.725054;
%Jxxmat=Rmat*[jxx1 0 0;0 jyy1 0; 0 0 jzz1]*inv(Rmat);
%Jab=0.5;
%Jab2=0.5;
%Jyymat=Jxxmat;
cochain.addmatrix('label','Jxx','value',Jxxmat,'color','r');
cochain.addmatrix('label','Jyy','value',Jyymat,'color','g');
cochain.addmatrix('label','Jperp','value',Jab*eye(3),'color','b');
cochain.addmatrix('label','Jperp2','value',Jab2*eye(3),'color','m');
cochain.addcoupling('mat','Jxx','bond',1,'subidx',1);
cochain.addcoupling('mat','Jxx','bond',1,'subidx',3);
cochain.addcoupling('mat','Jyy','bond',1,'subidx',2);
cochain.addcoupling('mat','Jyy','bond',1,'subidx',4);
cochain.addcoupling('mat','Jperp','bond',4);
cochain.addcoupling('mat','Jperp2','bond',3,'subidx',1);
cochain.addcoupling('mat','Jperp2','bond',3,'subidx',2);
cochain.addcoupling('mat','Jperp2','bond',3,'subidx',3);
cochain.addcoupling('mat','Jperp2','bond',3,'subidx',4);
plot(cochain,'range',[5 2 2],'atomMode','mag','cellMode','inside','atomLegend',false,'cellcolor','gray','bondMode','line','bondLinewidth0',2)
B=0.0*4.9/2.0;
n = [0.0 0.0 1.0];
cochain.field(n/norm(n)*B);
cochain.genmagstr('mode','direct', 'k',[0 0 0],'S',[0.1 0.2 0.23 0.12; 0 0 0 0;1 1 -1 -1]);
cochain.optmagsteep('random',true,'nRun',5e3);
cochain.mag_str.F
cochain.energy;
plot(cochain,'range',[5 2 2],'bondMode','line','bondLineWidth0',3)

clear('Qp');

Qp{1} = [ 0;   0.0; 0];
Qp{2} = [  1;   0.0; 0];
Qp{3} = [1;0.5;0.5];
Qp{4} = [1;1;1];
Qp{5} = [0;0;0];
Qp{6} = 101;
colormapeditor
spec1 = cochain.spinwave(Qp,'hermit',true,'formfact',false,'saveH',true,'saveT',true);
figure
sw_plotspec(spec1,'mode',1,'axLim',[0 4.5],'dashed', true)

spec1 = sw_neutron(spec1);
spec1 = sw_egrid(spec1,'component','Sperp','evect',linspace(0,4.5,2001),'imagChk',false);

f=figure('units','normalized','position',[.1 .1 .84 .5])
set(gca,'XTick',[], 'YTick', [])
set(gca,'XColor', 'none','YColor','none')
sw_plotspec(spec1,'mode','color','dE',0.1,'axLim',[0 1],'legend',false);

figure
E = linspace(0,4.5,200);
Q = linspace(0.0,3,600);
powSpec = cochain.powspec(Q,'Evect',E,'nRand',1e4,'formfact',true);
set(gca,'XTick',[], 'YTick', [])
set(gca,'XColor', 'none','YColor','none')
swpref.setpref('colormap',@jet)
sw_plotspec(powSpec,'dE',0.2,'axLim',[0 0.6]);
%csvwrite('/home/pavel/Documents/TFIM/powder_fit/fit_XYZ/Jab_05_Jab2_06.csv',powSpec.swConv);
