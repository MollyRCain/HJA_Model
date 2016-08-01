%
clear all
clc

forcing=load('calibrate_forcing.txt');
forcing_PrecCl=load('PrecipChloride_calibrate.txt');
forcing_QCl=load('DischargeChloride_calibrate.txt');
forcing_Temperature=load('Temperature_calibrate.txt');

          %     Imax Sumax    beta  Pmax   Kf    QL    Cp  Lp  Ks  SsMixVol  Tthresh  M
ParRange.minn = [0     300      0     .001  .03  0.001  0   0  0.01    100   -2      0.1];
ParRange.maxn = [3     600     100    4     1    0.7   1   1  0.03   5000     2      5];
%    Si, Su, Sf
Sin= [0   40  0 ];

%          Si,     Su,    Sf,    Ss   Sm
Cl_in = [800000 1000000 800000 700000  1];


ExtraPar.forcing=forcing(:,2:5);

ExtraPar.Sin=Sin;

ExtraPar.Cl_in=Cl_in;

ExtraPar.PrecCl=forcing_PrecCl(:,2);
ExtraPar.PrecClDate=forcing_PrecCl(:,1);

ExtraPar.QCl=forcing_QCl(:,2);
ExtraPar.QClDate=forcing_QCl(:,1);

ExtraPar.Temp=forcing_Temperature(:,2);

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
Prec=ExtraPar.forcing(:,1);     %mm/day
Qo=ExtraPar.forcing(:,2);       %mm/day
month=ExtraPar.forcing(:,4);
iwet=find(month==1 | month==2 | month==3 | month==4 | month==11 | month==12);
idry=find(month==5 | month==6 | month==7 | month==8 | month==9 | month==10);

RC_obs = sum(Qo)/sum(Prec);
RCobs_wet = sum(Qo(iwet))/sum(Prec(iwet));
RCobs_dry = sum(Qo(idry))/sum(Prec(idry));


%%
% GLUE
A=[];
nmax=100000;
h = waitbar(0,'Please wait...');
for n=1:nmax
Rnum=rand(1,12);
Par=ParRange.minn+(ParRange.maxn-ParRange.minn).*Rnum;
[~, ObjNS, ObjLNS,ObjClNS,Base,Fast,Ru,Rf,Rp,Rs,RC,~,~,~,~,~,~,~,~,~,~,~,~,RC_wet,RC_dry] = HBVMod_calibrate_snow(Par,ExtraPar);
if ObjNS>0.5 && ObjLNS>0.5  && ObjClNS>-10 %&& RC_dry<0.4 && RC_dry<0.7                                         
    A=[A;[Par ObjNS ObjLNS ObjClNS Base Fast Ru Rf Rp Rs RC]];
end
waitbar(n/nmax)

end
close(h)
save('MC.txt','A','-ascii');

[Opt,ind]=max(A(:,15));
OptPar=A(ind,1:12);

[Qm,Obj]=HBVMod_calibrate_snow(OptPar,ExtraPar);
[~,ObjNS,ObjLNS,ObjClNS,Base,~,~,~,~,~,RC,QtotClConc,SiF,SuF,SfF,SmF,SstotF,SiClConc,SuClConc,SfClConc,SsClConc,SmClConc,Cm,RC_wet,RC_dry]=HBVMod_calibrate_snow(OptPar,ExtraPar);
FinalCond = [SiF SuF SfF SmF SstotF];
FinalCondCl = [SiClConc SuClConc SfClConc SsClConc SmClConc];
save('QCl.mat','QtotClConc');