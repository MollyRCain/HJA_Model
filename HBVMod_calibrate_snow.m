function [Qm,ObjNS,ObjLNS,ObjClNS,Base,Fast,SumRu,SumRf,SumRp,SumRs,RC,QtotClConc,SiF,SuF,SfF,SmF,SstotF,SiClConc,SuClConc,SfClConc,SsClConc,SmClConc,Cm,RC_wet,RC_dry] = HBVMod_calibrate( Par, ExtraPar )
%HBVpareto Calculates values of 3 objective functions for HBV model



Imax=Par(1);        %Threshold of IR
Sumax=Par(2);       %Max UR storage?
beta=Par(3);        %Shape parameter of runoff generation
Pmax=Par(4);        %Max percolation rate
Kf=Par(5);          %FR timescale
QL=Par(6);          %Constant deep iNS loss
Cp=Par(7);          %Pref recharge coeff
Lp=Par(8);          %Transpiration threshold
Ks=Par(9);          %SR timescale
SsMixVol=Par(10);   %GW mixing volume
Tthresh=Par(11);     %Treshold temperature for snow [deg]
M=Par(12);          %Melting parameter [mm/deg*day]

Prec=ExtraPar.forcing(:,1);     %mm/day
Qo=ExtraPar.forcing(:,2);       %mm/day
Etp=ExtraPar.forcing(:,3);      %mm/day
month=ExtraPar.forcing(:,4);
PConc = ExtraPar.PrecCl;        %mg/mm
QClO=ExtraPar.QCl;          %mg/mm
QClODate=ExtraPar.QClDate;
QClmodDate=ExtraPar.PrecClDate;
Temp=ExtraPar.Temp;         %deg C

indMod=ExtraPar.indMod;
indObs=ExtraPar.indObs;

tmax=length(Prec);
Si=zeros(tmax,1);       %interception reservoir storage
Su=zeros(tmax,1);       %unsaturated soil reservoir storage
Sf=zeros(tmax,1);       %fast reacting reservoir storage
Sstot=zeros(tmax,1);    %slow reacting reservoir storage
Sm=zeros(tmax,1);       %snow reservoir storage
Eidt=zeros(tmax,1);
Eudt=zeros(tmax,1);
Qtotdt=zeros(tmax,1);
Qsdt=zeros(tmax,1);
Qfdt=zeros(tmax,1);
Rs=zeros(tmax,1);
Rf=zeros(tmax,1);
Rp=zeros(tmax,1);
Ru=zeros(tmax,1);
SiCl=zeros(tmax,1);
SiClConc=zeros(tmax,1);
SuCl=ones(tmax,1);
SuClConc=ones(tmax,1);
SfCl=zeros(tmax,1);
SfClConc=zeros(tmax,1);
SsCl=zeros(tmax,1);
SsClConc=zeros(tmax,1);
SmCl=zeros(tmax,1);
SmClConc=zeros(tmax,1);
QfCl=zeros(tmax,1);
QsCl=zeros(tmax,1);
QtotClConc=zeros(tmax,1);
QtotCl=zeros(tmax,1);
QsClConc=zeros(tmax,1);
QfClConc=zeros(tmax,1);
Cm=zeros(tmax,1);
SiClConc_R=zeros(tmax,1);
SmClConc_R=zeros(tmax,1);
Pedt=zeros(tmax,1);
PeCl=zeros(tmax,1);
PeClConc=zeros(tmax,1);
Ridt=zeros(tmax,1);
Rmdt=zeros(tmax,1);

Si(1)=ExtraPar.Sin(1);
Su(1)=ExtraPar.Sin(2);
Sf(1)=ExtraPar.Sin(3);
Sstot(1)=0;
Sm(1)=0;

SiClConc(1)=ExtraPar.Cl_in(1);
SiCl(1)=SiClConc(1)*Si(1);
SiCl_input=SiCl(1);
SuClConc(1)=ExtraPar.Cl_in(2);
SuCl(1)=SuClConc(1)*Su(1);
SuCl_input=SuCl(1);
SfClConc(1)=ExtraPar.Cl_in(3);
SfCl(1)=SfClConc(1)*Sf(1);
SfCl_input=SfCl(1);
SsClConc(1)=ExtraPar.Cl_in(4);
SsCl(1)=SsClConc(1)*(Sstot(1)+SsMixVol);
SsCl_input=SsCl(1);
SmClConc(1)=ExtraPar.Cl_in(5);
SmCl(1)=SmClConc(1)*Sm(1);
SmCl_input=SmCl(1);

dt=1;       %hour

%%
% 
for i=1:tmax
    
    Pdt=Prec(i)*dt;                     %dt is 1 so just Precip
    Epdt=Etp(i)*dt;                     %dt is 1 so just Potential Evap
    
    if Temp(i)>Tthresh
        % Snow Reservoir
        SmClConc_R(i)=SmClConc(i);
        Rmdt(i)=min(Sm(i),M*(Temp(i)-Tthresh)*dt);                 %Melt water
        Sm(i)=Sm(i)-Rmdt(i);
        SmCl(i)=SmCl(i)-Rmdt(i)*SmClConc(i);
%         if Sm(i)>0
%             SmClConc(i)=SmCl(i)/Sm(i);
%         else
%             SmClConc(i)=0;
%         end
        
        % Interception Reservoir
        if Pdt>0                            %If it rains
            Si(i)=Si(i)+Pdt;                %SI increases by precip amount
            SiCl(i)=SiCl(i)+Pdt*PConc(i);
            if Si(i)>0
                SiClConc(i)=SiCl(i)/Si(i);
            else
                SiClConc(i)=0;
            end

            SiClConc_R(i)=SiClConc(i);

            Ridt(i)=max(0,Si(i)-Imax);         %Effective rainfall (spills from Si when exceeds Imax; percolation; exludes runoff)
            RidtCl=Ridt(i)*SiClConc(i);
            Si(i)=Si(i)-Ridt(i);               %New Si, either Si before or with overflow subtracted to = Imax
            SiCl(i)=SiCl(i)-Ridt(i)*SiClConc(i);
            if Si(i)>0
                SiClConc(i)=SiCl(i)/Si(i);
            else
                SiClConc(i)=0;
            end


            Eidt(i)=0;                      %Assume no actual evap when rainfall                    
        else                                %Evaporation only when there is no raiNSall
            Ridt(i)=0;                         %Effective rainfall is 0
            Eidt(i)=min(Epdt,Si(i));        %Evap from Si is potential evap, unless exceeds Si
            Si(i)=Si(i)-Eidt(i);            %New Si is old Si minus evap
            if Si(i)>0
                SiClConc(i)=SiCl(i)/Si(i);
            else
                SiClConc(i)=0;
            end

            %SiClConc_R(i)=SiClConc(i);
        end
        
    else
        %Snow Reservoir
        Sm(i)=Sm(i)+Pdt;
        SmCl(i)=SmCl(i)+Pdt*PConc(i);
        if Sm(i)>0
            SmClConc(i)=SmClConc(i)/Sm(i);
        else
            SmClConc(i)=0;
        end
        
        SmClConc_R(i)=SmClConc(i);
        Rmdt(i)=0;
        
        % Interception Reservoir
        Ridt(i)=0;                         %Effective rainfall is 0
        Eidt(i)=min(Epdt,Si(i));        %Evap from Si is potential evap, unless exceeds Si
        Si(i)=Si(i)-Eidt(i);            %New Si is old Si minus evap
        if Si(i)>0
            SiClConc(i)=SiCl(i)/Si(i);
        else
            SiClConc(i)=0;
        end
        
        SiClConc_R(i)=SiClConc(i);
        
    end
        
    if i<tmax                           
        Si(i+1)=Si(i);       
        SiCl(i+1)=SiCl(i);
        SiClConc(i+1)=SiClConc(i);
        
        Sm(i+1)=Sm(i);       
        SmCl(i+1)=SmCl(i);
        SmClConc(i+1)=SmClConc(i);
    end
    
    
    if i>1
        ClBSiSm=Prec(i)*PConc(i)-Ridt(i)*SiClConc_R(i)-Rmdt(i)*SmClConc_R(i)+SiCl(i-1)-SiCl(i)+SmCl(i-1)-SmCl(i); 
        if ClBSiSm>0.001
            fprintf('ClBSiSm = %d\n',round(ClBSiSm));
        end
    end
    
    if i>1
        WBSiSm=Prec(i)-Ridt(i)-Rmdt(i)-Eidt(i)+Si(i-1)-Si(i)+Sm(i-1)-Sm(i);
        if WBSiSm>0.001
            fprintf('WBSiSm = %d\n',round(WBSiSm));
        end
    end
    
    
    % Unsaturated Reservoir
    %inflows and outflows
    Pedt(i)=Ridt(i)+Rmdt(i);
    PeCl(i)=Ridt(i)*SiClConc_R(i)+Rmdt(i)*SmClConc_R(i);
    
    if Pedt(i)>0
        PeClConc(i)=PeCl(i)/Pedt(i);
    else
        PeClConc(i)=0;
    end
    
    
    CR= 1/(1+exp(((-Su(i)/Sumax)+0.5)/beta));      %determines how Pe is partitioned to Su and Sf
    Ru(i) = (1-CR)*Pedt(i);                           %portion of Pe that discharges to Su
    Rf(i)= CR*(1-Cp)*Pedt(i);                         %recharge of fast reservoir
    Rp(i)=CR*Cp*Pedt(i);                              %preferential recharge of slow reservoir
    Su(i)=Su(i)+ Ru(i);                            %new Su is Su plus inflow from Pe
    
    Cm(i)=1-(Su(i)/Sumax);                           %mixing coefficient
%     if Cm(i)<0
%         Cm(i)=0;
%     end
    
    RuClB_in=Ru(i)*PeClConc(i);
    RpClB_in=Cm(i)*Rp(i)*PeClConc(i);
    RfClB_in=Cm(i)*Rf(i)*PeClConc(i);
    
    SuCl(i)=SuCl(i)+RuClB_in;
    SuCl(i)=SuCl(i)+RpClB_in+RfClB_in;
    
   
    % Transpiration
    Eudt(i)=(Epdt-Eidt(i))*min(1,(Su(i)/(Sumax.*Lp)));     %Evap from unsat res
    Eudt(i)=min(Eudt(i),Su(i));                         %If actual evap exceeds Su, then = Su
    Su(i) = Su(i) - Eudt(i);                            %New Su is Su minus transpiration
    
    if Su(i)>0
        SuClConc(i)=SuCl(i)/(Su(i)+(Cm(i)*(Rp(i)+Rf(i))));
    else
        SuClConc(i)=0;
    end
        
    Rs(i)=Pmax*(Su(i)/Sumax)*dt;           %recharge of slow reservoir
    Su(i)=Su(i)-Rs(i);
    
    RsClB_out=Rs(i)*SuClConc(i);
    RpClB_out=Cm(i)*Rp(i)*SuClConc(i);
    RfClB_out=Cm(i)*Rf(i)*SuClConc(i);
    
    SuCl(i)=SuCl(i)-RsClB_out-RpClB_out-RfClB_out;
    
    
    if Su(i)>0
        SuClConc(i)=SuCl(i)/Su(i);
    else
        SuClConc(i)=0;
    end
    
   
    if i<tmax                           
        Su(i+1)=Su(i);
        SuCl(i+1)=SuCl(i);
        SuClConc(i+1)=SuClConc(i);
    end
    
    if i>1
        ClBSu= RuClB_in+RpClB_in+RfClB_in-RsClB_out-RpClB_out-RfClB_out+SuCl(i-1)-SuCl(i); 
%         fprintf('ClBSu = %d\n',round(ClBSu));
    end
    
     
% Fast Reservoir
    Sf(i)=Sf(i)+Rf(i)*dt;                                %Sf old plus portion of Pe that discharges to 
    
    RfClB_in=RfClB_out+((1-Cm(i))*Rf(i)*PeClConc(i));
    SfCl(i)=SfCl(i)+RfClB_in;
    
    if Sf(i)>0
        SfClConc(i)=SfCl(i)/Sf(i);
    else
        SfClConc(i)=0;
    end
    
    Qfdt(i)= min(Sf(i),(Sf(i)*(1-exp(-Kf*dt)))/dt);      %Fast discharge. i=t
    
    Sf(i)=Sf(i)-Qfdt(i);                                 %New Sf subtracting out discharge
    
    QfCl(i) = Qfdt(i)*SfClConc(i);
    SfCl(i)=SfCl(i)-Qfdt(i)*SfClConc(i);
    
    if Sf(i)>0
        SfClConc(i)=SfCl(i)/Sf(i);
    else
        SfClConc(i)=0;
    end
    
    if i<tmax
        Sf(i+1)=Sf(i);
        SfCl(i+1)=SfCl(i);
        SfClConc(i+1)=SfClConc(i);
    end
    
    if i>1
        ClBSf= RfClB_in-QfCl(i)+SfCl(i-1)-SfCl(i); 
%        fprintf('ClBSf = %d\n',round(ClBSf));
    end
    
    
    % Slow Reservoir
    Sstotin=Sstot(i)+Rs(i)*dt+Rp(i)*dt;                                                                         %Sstot after inflows of pref recharge of slow res and recharge of slow res
   
    RsRpClB_in=RsClB_out+RpClB_out+((1-Cm(i))*Rp(i)*PeClConc(i));
    SsCl(i)=SsCl(i)+RsRpClB_in;
    
    SsClConc(i)=SsCl(i)/(Sstotin+SsMixVol);
    
    if Sstotin>0
        Sstotout = (Sstotin*exp(-Ks*dt))-((QL/Ks)*(1-exp(-Ks*dt)));                                             %Sstot after inflows and outflows
    else Sstotout = Sstotin-QL;
    end
    
    Qstot=(Sstotin-Sstotout)/dt;                                                                                %Total discharge from slow res
    Qsdt(i)= max(0,(Qstot-QL));                                                                                 %runoff from slow reservoir (subtracting deep percolation)
    
    Ssa = max(0,Sstotout);
    Ssp = min(0,Sstotout);
    
    QsCl(i)=Qsdt(i)*SsClConc(i);
    SsCl(i)=SsCl(i)-QL*SsClConc(i)-Qsdt(i)*SsClConc(i);
    QLClB(i)=QL*SsClConc(i);
    SsClConc(i)=SsCl(i)/(Sstotout+SsMixVol);
    
    Sstot(i)= Ssa+Ssp;
    
    if i<tmax
        Sstot(i+1)=Sstot(i);
        SsCl(i+1)=SsCl(i);
        SsClConc(i+1)=SsClConc(i);
    end
    
    if i>1
        ClBSs=RsRpClB_in-QLClB(i)-QsCl(i)+SsCl(i-1)-SsCl(i);
%        fprintf('ClBSs = %d\n',round(ClBSs));
    end
    

    if i==1
        ClB1=Prec(i)*PConc(i)-QfCl(i)-QsCl(i)-QLClB(i)+SiCl_input+SuCl_input+SfCl_input+SsCl_input+SmCl_input-SiCl(i)-SuCl(i)-SfCl(i)-SsCl(i)-SmCl(i);
    else
        ClB1=Prec(i)*PConc(i)-QfCl(i)-QsCl(i)-QLClB(i)+SiCl(i-1)+SuCl(i-1)+SfCl(i-1)+SsCl(i-1)+SmCl(i-1)-SiCl(i)-SuCl(i)-SfCl(i)-SsCl(i)-SmCl(i);
%        fprintf('ClB = %d\n',round(ClB));
    end
    
%     if ClB1>0.01
%         fprintf('ClB1 = %d\n',round(ClB1));
%     end
    
    
    %%
    Cl_input=SuCl_input+SsCl_input+SfCl_input+SiCl_input+SmCl_input;
    Cl_output=SiCl(i)+SuCl(i)+SfCl(i)+SsCl(i)+SmCl(i);
    ClB2=sum(Prec(1:i).*PConc(1:i))-sum(QfCl(1:i))-sum(QsCl(1:i))-sum(QLClB(1:i))+Cl_input-Cl_output;

%     if ClB2>0.01
%         fprintf('ClB2 = %d\n',round(ClB2));
%     end
    

    
    Qtotdt(i)=Qsdt(i)+Qfdt(i);                        %discharge to stream from fast and slow reservoirs
    QtotCl(i)=QfCl(i)+QsCl(i);
    QsClConc(i)=max(0,QsCl(i)/Qsdt(i));
    QfClConc(i)=QfCl(i)/Qfdt(i);
    QtotClConc(i)=max(0,QtotCl(i)/Qtotdt(i));
    
end 

% Check Water Balance
Sen=Si(tmax)+Su(tmax)+Sf(tmax)+Sstot(tmax)+Sm(tmax);                         %Final total storage
SiF=Si;
SuF=Su;
SfF=Sf;
SmF=Sm;
SstotF=Sstot;

Sin=sum(ExtraPar.Sin);                                              %Initial total storage
WB=sum(Prec)-sum(Eidt)-sum(Eudt)-sum(Qtotdt)-(QL*tmax)+Sin-Sen;     %Precip - Evap - Transp - Discharge - deep percolation + initial S - final S
if WB>0.01
    fprintf('WB = %d\n',round(WB));                                     %Should = 0
end

%%Check Chloride Balance

Cl_input=SuCl_input+SsCl_input+SfCl_input+SiCl_input+SmCl_input;
%Cl_input=SiCl(1)+SuCl(1)+SfCl(1)+SsCl(1);
Cl_output=SiCl(end)+SuCl(end)+SfCl(end)+SsCl(end)+SmCl(end);
ClB=sum(Prec.*PConc)-sum(QfCl)-sum(QsCl)-sum(QLClB(1:i))+Cl_input-Cl_output;

if ClB>0.01
    fprintf('ClB = %d\n',round(ClB));
end

%%
Base=sum(Qsdt)/sum(Qtotdt);
Fast=sum(Qfdt)/sum(Qtotdt);
SumRu=sum(Ru);
SumRf=sum(Rf);
SumRp=sum(Rp);
SumRs=sum(Rs);

RC=mean(Qtotdt)/mean(Prec);
iwet=find(month==1 | month==2 | month==3 | month==4 | month==11 | month==12);
idry=find(month==5 | month==6 | month==7 | month==8 | month==9 | month==10);
RC_wet=mean(Qtotdt(iwet))/mean(Prec(iwet));
RC_dry=mean(Qtotdt(idry))/mean(Prec(idry));

% % Offset Q
% Weigths=Weigfun(Tlag);
% Qm = conv(Qtotdt,Weigths);
% Qm=Qm(1:tmax);
Qm = Qtotdt;

% Calculate objective
ind1=find(Qo>=0);
ind2=find(ind1>365);                            %only use data after first year to calculate Obj
ind3=ind1(ind2);
QoAv=mean(Qo(ind3));

ErrUp=sum((Qm(ind3)-Qo(ind3)).^2);              %Nash-Sutcliff Coeff
ErrDo=sum((Qo(ind3)-QoAv).^2);
ObjNS=1-ErrUp/ErrDo;

ErrUp=sum((log(Qm(ind3))-log(Qo(ind3))).^2);      %Log Nash-Sutcliff Coeff
ErrDo=sum((log(Qo(ind3))-log(QoAv)).^2);
ObjLNS=1-ErrUp/ErrDo;

%% Calculate objective for Chloride
QClOAv=mean(QClO(indObs));

% ErrDo=sum((QCLO(indObs)-QtotClConc(indMod)).^2);              %Coeff of determination
% ErrUp=sum((QClO(indObs)-QClOAv).^2);
% ObjCl=ErrUp/ErrDo;

ErrUp=sum((QtotClConc(indMod)-QClO(indObs)).^2);              %Nash-Sutcliff Coeff
ErrDo=sum((QClO(indObs)-QClOAv).^2);
ObjClNS=1-ErrUp/ErrDo;

% ErrUp=sum((log(QtotClConc(indMod))-log(QClO(indObs))).^2);      %Log Nash-Sutcliff Coeff
% ErrDo=sum((log(QClO(indObs))-log(QClOAv)).^2);
% ObjClLNS=1-ErrUp/ErrDo;


% %% Plot
% hour=1:tmax;
% plot(hour,Qo,'g');
% hold on
% plot(hour,Qm,'r');
% legend('Qobs','Qmod');
