close all
clc

forcing=load('calibrate_forcing.txt');
forcing_PrecCl=load('PrecipChloride_calibrate.txt');
forcing_QCl=load('DischargeChloride_calibrate.txt');
forcing_Temperature=load('Temperature_calibrate.txt');

%    Si, Su, Sf
Sin= [0   40  0 ];

%          Si,     Su,    Sf,    Ss   Sm
Cl_in = [800000 1000000 800000 700000  1];

ExtraPar.forcing=forcing(:,2:5);
ExtraPar.Sin=Sin;
ExtraPar.Cl_in=Cl_in;

ExtraPar.Cl=forcing_PrecCl(:,2);

ExtraPar.PrecCl=forcing_PrecCl(:,2);
ExtraPar.PrecClDate=forcing_PrecCl(:,1);
ExtraPar.QCl=forcing_QCl(:,2);
ExtraPar.QClDate=forcing_QCl(:,1);
ExtraPar.Temp=forcing_Temperature(:,2);

Qo=ExtraPar.forcing(:,2);
Prec=ExtraPar.forcing(:,1);
tmax=length(Prec);
hour=1:tmax;

A=load('MC.txt');

ParetoNS = 1-A(:,13);
ParetoLNS= 1-A(:,14);

figure(1)
plot(ParetoNS,ParetoLNS,'o')
xlabel('1-Nash-Sutcliffe');
ylabel('1-Log Nash-Sutcliffe');

%%  Determines date indices to use for creating a vector of QClmod to be compare to QClO
QClO=ExtraPar.QClDate;          %mg/mm
QClODate=ExtraPar.QClDate;
QClModDate=ExtraPar.PrecClDate;

indices=zeros(length(QClODate),1);
for i=1:length(QClODate)
    indices(i) = find(QClModDate==QClODate(i));
end

endyeardate=QClModDate(365);
endyearrow=find(QClODate>endyeardate,1);

indMod=indices([endyearrow:end],:);
indObs=transpose(endyearrow:length(QClODate));

ExtraPar.indMod=indMod;
ExtraPar.indObs=indObs;

%%
[Opt,ind]=max(A(:,15));
OptPar=A(ind(1),1:12);

[Qm,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,SiClConcF,SuClConcF,SfClConcF,SsClConcF] = HBVMod_calibrate_snow(OptPar,ExtraPar);
%[Qm,ObjNS,ObjLNS,ObjClNS,Base,Fast,SumRu,SumRf,SumRp,SumRs,RC,QtotClConc,SiF,SuF,SfF,SstotF,SiClConcF,SuClConcF,SfClConcF,SsClConcF,Cm,RC_wet,RC_dry]

QClObs=forcing_QCl(:,2);
DateQClObs=forcing_QCl(:,1);

load('QCl.mat');                    %modeled QCl called QtotClConc
DateQClmod=forcing_PrecCl(:,1);

figure(6)
plot(hour,Qo,'r');
hold on
plot(hour,Qm,'g');
xlabel('time [days]');
ylabel('Q [mm/day]');
yyaxis right
plot(DateQClObs-DateQClmod(1),QClObs,'or:');
ylabel('Q Cl [mg/mm]');
legend('Qobs','Qmod','QClObs');

%%
figure(3)
subplot(6,2,1)
plot(A(:,1),A(:,13),'.');
xlabel('I_{max}');
ylabel('NS Q');

subplot(6,2,2)
plot(A(:,2),A(:,13),'.');
xlabel('S_{u,max}');
ylabel('NS Q');

subplot(6,2,3)
plot(A(:,3),A(:,13),'.');
xlabel('\beta');
ylabel('NS Q');

subplot(6,2,4)
plot(A(:,4),A(:,13),'.');
xlabel('P_{max}');
ylabel('NS Q');

subplot(6,2,5)
plot(A(:,5),A(:,13),'.');
xlabel('K_{f}');
ylabel('NS Q');

subplot(6,2,6)
plot(A(:,6),A(:,13),'.');
xlabel('Q_{L}');
ylabel('NS Q');

subplot(6,2,7)
plot(A(:,7),A(:,13),'.');
xlabel('Cp');
ylabel('NS Q');

subplot(6,2,8)
plot(A(:,8),A(:,13),'.');
xlabel('Lp');
ylabel('NS Q');

subplot(6,2,9)
plot(A(:,9),A(:,13),'.');
xlabel('K_{S}');
ylabel('NS Q');

subplot(6,2,10)
plot(A(:,10),A(:,13),'.');
xlabel('MixVol');
ylabel('NS Cl');

subplot(6,2,11)
plot(A(:,11),A(:,13),'.');
xlabel('Ttresh');
ylabel('NS Cl');

subplot(6,2,12)
plot(A(:,12),A(:,13),'.');
xlabel('M');
ylabel('NS Cl');


%%
QClObs=forcing_QCl(:,2);
DateQClObs=forcing_QCl(:,1);

load('QCl.mat');                    %modeled QCl called QtotClConc
DateQClmod=forcing_PrecCl(:,1);

figure(7)
plot(DateQClmod-DateQClmod(1),QtotClConc,'g');
hold on
plot(DateQClObs-DateQClmod(1),QClObs,'or:');
legend('QClmod','QClobs');
ylabel('mg/mm');
%ylim([0 15000000]);

%%
figure(8)
plot(DateQClmod-DateQClmod(1)+1,SuClConcF,DateQClmod-DateQClmod(1)+1,SfClConcF,DateQClmod-DateQClmod(1)+1,SsClConcF,DateQClmod-DateQClmod(1)+1,QtotClConc);
legend('SuCl','SfCl','SsCl','QtotClConc');
%ylim([0 15000000]);
