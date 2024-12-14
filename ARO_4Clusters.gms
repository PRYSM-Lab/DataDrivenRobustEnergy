
$set outputname monolithic_REW_uncert

option threads=36;
option optcr= 0.05;
option MIP= gurobi;
*gurobi;
option rmip=gurobi;
option limrow=100000;
option reslim=172800;
option lp=gurobi;
sets
g                'Regions'
l                'Transportation modes'
p                'Prodution tehnologies'
r                'Reservoir'
s                'Storage tehnologies'
t                'years'
d                'Diameter size'
c                'clusters'
h                'hours  /1*8760/' 
e                'renewables'
j                '24h'
;

$onecho>gurobi.opt
         threads 36
         mipstart=1
$offecho

*$onecho>gurobi.opt
*numericfocus 1
*scaleflag 1
*$offecho

sets
sc(s)            'storage caverns'
sv (s)           'storage vessels'

alias (g,g1)
alias (h,h1)
;

set
N(g,g1)          'set of neighbouring regions g and g1'
Npipe(g,g1)      'set of pipeline connections between region g and g1'
GR(g,r)          'set of collection points g and reservoir r connections'
GS(g,s)          'set of regions g in which storage caverns sc are located'
Gimp(g)          'set of regions g in which international import can take place'
;

Parameters

beta                     'ratio of stored amount (%)'
y_c(p,t)                 'coefficient of CO2 capture for plant type p in time period t(tn CO2 / MWh H2)'
y_e(p,t)                 'coefficient of CO2 emission for plant type p in time period t(tn CO2 / MWh H2)'
deltaH                   'Ratio of hydrogen regional pipeline operating costs to capital costs (%)'
deltaC_onshore           'ratio of onshore CO2 pipeline operating costs to capital costs'
deltaC_offshore          'ratio of offshore CO2 pipeline operating costs to capital costs'
diaH(d)                  'diameter of a regional hydrogen pipeline of diameter size d (m)'
diaC_onshore(d)          'diameter of an onshore CO2 pipeline of diameter size d (m)'
diaC_offshore(d)         'diameter of an offshore CO2 hydrogen pipeline of diameter size d (m)'
iota                     'maximum percentage of international hydrogen imports over the total demand (%)'
dur                      'duration of time periods (y)'
LTonshore                'useful life of onshore CO2 pipelines (y)'
LToffshore               'useful life of offshore CO2 pipelines (y)'
LTpipe                   'useful life of hydrogen pipelines (y)'
LTp(p)                   'useful life of hydrogen production plants (y)'
LTs(s)                   'useful life of hydrogen storage facilities (y)'
LTt(l)                   'useful life of hydrogen road transportation modes {Trailer,Tanker} (y)'
a                        'days in a year (days)'
AV(c,h,g,e)              'availability of renewable e in region g, cluster c and hour h (%)'
ayHR0(d,g,g1)            'initial availability of a regional hydrogen pipeline of diameter size d between regions g and g1 (0-1)'
ayC0(d,g,g1)             'initial availability of a onshore CO2 pipeline of diameter size d between regions g and g1 (0-1)'
aeC0(r)                  'initial availability of a offshore CO2 pipeline of diameter size d between collection point in regions g and reservoir r (0-1)'
BA(g,t)                  'biomass availability in region g and time period t'
br(g)                    'biomass regional availability'
cbio(t)                  'biomass cost in time period t (?/MWh)'
cccH(d)                  'capital costs of a regional hydrogen pipeline of diameter size q d (?k km-1)'
cccC_onshore(d)          'capital costs of an onshore CO2 pipeline of diameter size d (?k km-1)'
cccC_offshore(d)         'capital costs of an offshore CO2 pipeline of diameter size d (?k km-1)'
cgas(t)                  'natural gas cost in time period t (?/MWh)'
crf                      'capital recovery factor'
Cstart(p)                'cost for starting up for each technology type (?)'
Cshut(p)                 'cost for shutting down for each technology type (?)'
ct(t)                    'carbon tax in time period t (? kg-1 CO2)'
dc(t)                    'demand coefficient at time period t'
dem(g,t,c,h)             'total hydrogen demand in region g in time period t (MW)'
dfc(t)                   'discount factor for capital costs in time period t'
dfo(t)                   'discount factor for operating costs in time period t'
dw(l)                    'driver wage of road transportation mode l (? h-1)'
DistPipe(g,g1)           'delivery distance of a onshore CO2 pipeline between regions g and g1 (km)'
DistRes(g,r)             'Distance from CO2 collection point in region g to reservoir r (km)'
Dist(g,g1)           'regional delivery distance of hydrogen transportation mode l in region g (km)'
DistSt(g,s)              'distance between region g and underground storage type s'
DT(p)                    'min down time (h)'
ec(t)                    'cost of electricity back to grid (?/MWe)'
eta(p,t)                 'efficiency of WE in time period t (%)'
emtarget(t)              'emissions target in time period t (kgCO2)'
feR(l)                   'fuel economy of road transportation mode l transporting product type i within a region (km l-1)'
fp(l)                    'fuel price of road transportation mode l (? l-1)'
GasDem(c,h,t,g)            'hydrogen demand for each region g each cluster c and hour h (MWh)'
ge(l)                    'general expenses of road transportation mode l transporting product type i (? d-1)'
ir                       'discount rate (%)'
landAV(e,g)              'land availability of renewable e in region g (MW)'
lut(l)                   'load and unload time of road transportation mode l  (h)'
me(l)                    'maintenance expenses of road transportation mode l  (? km-1)'
nel                      'economic life cycle of capital investments (y)'
np0(p,g)                 'initial number of hydrogen produCtion plants of teChnology p and size j producing product type i in region g'
ns0(s,g)                 'initial number of hydrogen storage facilities of type s and size j storing product type i in region g'
pcap_max(p)              'maximum capacity of a hydrogen production plant of type p and size j (MW)'
pcap_min(p)              'manimum capacity of a hydrogen production plant of type p and size j (MW)'
pccost(p,t)              'capital cost of a production plant of type p (?k/kW)'
pimp                     'Price of hydrogen import (?/MWh)'
pocostF(p,t)             'operating production cost  in a production plant of type p  (?/MW/y)'
pocostV(p,t)             'operating production cost  in a production plant of type p  (?/MW)'
qHmax(d)                 'maximum flowrate in a hydrogen pipeline of diameter size d (kg H2 d-1)'
qCmax(d)                 'maximum flowrate in a CO2 pipeline of dimetere size d (kg H2 d-1)'
QImax(s)                 'maximum injection rate for each storage type s'
QRmax(s)                 'maximum retrieval rate for each storage type s'
rcap(r)                  'total capacity of reservoir r (kg CO2-eq)'
ri0                      'initial CO2 inventory in reservoir r (kg CO2)'
RD(p)                    'Comit Ramp down '
rccost(e,t)              'renewable e cacital cost in time period t(?/MW)'
rocost(e,t)              'renewable e operating cost in time period t(?/MW)'
RU(p)                    'Comit Ramp up '
scap_max(s)              'maximum capacity of a storage facility of type s (MWh H2)'
scap_min(s)              'minimum capacity of a storage facility of type s (MWh H2)'
sccost(s)                'capital cost of a storage facility of type s (?/MW)'
socostF(s)               'Fixed operating storage cost  in a production plant of type p (?/MW/y)'
socostV(s)               'Variable operating storage cost  in a production plant of type p (£/kWh stored)'
spR(l)                   'regional average speed of road transportation mode l transporting product type i within a region(km h-1)'
st0(s,g)                 'storage at time 0'
tcap(l)                  'capacity of road transportation mode l transporting product type i (MWh H2 unit-1)'
tmc(l)                   'capital cost of establishing a road transportation unit of transportation mode l (? unit-1)'
tmaR(l)                  'regional availability of road transportation mode l (h d-1)'
PCap(p)                  'unit capacity for production type p (MW)'
SCap(s)                  'unit capacity for storage type s (MW)'
uInit(p,g,t)             'Initial operating units type p in region g at time period t'
UT(p)                    'min up time (h)'
Vbio_max(t)              'maximum biomass consumption in year t'
WF(c)                    'weight of clusters (d)'
demH(g,t,c,h)            'demand in hour intervals (MW)'
AVH(c,h,g,e)             'availability of renewable e in region g, cluster c and hour h (%)'
BigQ_neg(c,h,j,t,g)
BigQ_pos(c,h,j,t,g)
Gasdemave(c,h,t,g) 
;

$call gdxxrw newhydro_clusters4.xlsx skipempty=0 trace=3 index=Index!A1:F200
$gdxin newhydro_clusters4.gdx
$load p,h,g,l,r,c,s,sc,sv,t,d,e,N,GR,GS,Gimp,j
$load beta ,y_c, y_e,dc,deltaH ,deltaC_onshore,deltaC_offshore ,diaH,diaC_onshore
$load diaC_offshore, dur , ec,emtarget, DT,LTonshore ,LToffshore ,LTpipe ,LTp, LTs, LTt
$load a,AV,ayHR0, ayC0,aeC0,br,cbio,cccH,cccC_onshore,cccC_offshore,ct,eta
$load feR,fp,GasDem,ge,landAV,lut,me,nel,np0,ns0,pcap_max,dw,DistPipe,DistRes,Dist,DistSt,PCap,pimp,Gasdemave
$load pcap_min,pccost,pocostF,pocostV,qHmax,qCmax,QImax,QRmax,rcap,rccost,rocost,ri0,RD,RU,scap_max,SCap,scap_min,sccost
$load socostF,socostV,spR,tcap,tmc,tmaR,UT,uInit,iota,Vbio_max
$load WF,BigQ_neg,BigQ_pos
$load Cstart,Cshut
;

sets
TT(t)   /3*6/
CC(c)   /1*4/
HH(h)   /1*24/;

scalar
y1       /3/
y2       /6/
theta    /1/;

ec(t)=3*ec(t);

cgas('1')=26.62;
cgas('2')=28.33;
cgas('3')=28.33;
cgas('4')=28.33;
cgas('5')=29.04;
cgas('6')=29.78;

ir=0.06;

*regional biomass availability
parameter bp;
bp=0.5;
display bp;
BA(g,t)=bp*br(g)*Vbio_max(t)*1000000;

*discount factor of capital costs
dfc(t)=round(1/(1+ir)**(dur*ord(t)-dur),2);
*discount factor of operating costs
dfo(t)=round(1/(1+ir)**(dur*ord(t)-5) + 1/(1+ir)**(5*ord(t)-4) + 1/(1+ir)**(5*ord(t)-3)+
1/(1+ir)**(5*ord(t)-2) + 1/(1+ir)**(5*ord(t)-1),2);
*capital recovery factor
nel=30;
crf=round((ir*(1+ir)**nel)/((1+ir)**nel-1),2);

*SET Npipe DEFINITION
Npipe(g,g1)$(DistPipe(g,g1) gt 0)=yes;
display N,Npipe;

*Set GS
st0(s,g)=0;
GS(g,sv)=yes;

*demand
*dem(g,t,c,h)= dc(t)*GasDem(c,h,g);
dem(g,t,c,h)= GasDem(c,h,t,g);

display dem;
display g,l,p,r,s,c,h,t,d,e,N,GR,GS,Gimp,beta ,y_c,y_e,Cstart, Cshut, deltaH ,deltaC_onshore,deltaC_offshore ,diaH,diaC_onshore,
diaC_offshore, dc,dur,LTonshore ,LToffshore ,LTpipe ,LTp, LTs, LTt,
ayHR0, AV,ayC0,aeC0,BA,br,cccH,cccC_onshore,cccC_offshore,cbio,cgas,crf,ct,dem,dfc,
dfo,dw,DistPipe,DistRes,Dist,DistSt,ec,emtarget,eta,feR,fp,GasDem,ge,ir,landAV,lut,me,nel,np0,ns0,pcap_max,
pcap_min,pccost,pocostF,pocostV,pimp,qHmax,qCmax,rcap,rccost,rocost,ri0,scap_max,scap_min,sccost,
socostF,socostV,spR,tcap,tmc,tmaR,pcap,scap,Vbio_max,WF;

set it outer iterations /it1*it8/;
set iter(it);
set iter_fesi(it);
set Stage2_it /Stage2_it1*Stage2_it30/;

iter(it) =no;
iter_fesi(it)=no;


integer Variables
InvP(p,g,t)      'Investment of new plants of type p producing in region g in time period t'
InvS(s,g,t)      'Investment of new storage facilities of type s in region g in time period t'
ITU(l,g,g1,t)    'Number of new transportation units of type l for regional transportation byroad in region g to region g acquired in time period t'
;

positive Variables
NP(p,g,t)        'Number of plants of type j and size p in region g in time period t'
NS(s,g,t)        'Number of storage facilities of type s and size p in region g in time period t'
NTU(l,g,g1,t)    'Number of transportation units of type l for regional transportation by road in region g in time period t'
;

positive Variables
AY(d,g,g1,t)     'availability of hydrogen pipelines of diameter size d for regional distribution in region g in time period t'
AYon(d,g,g1,t)   'availability of onshore CO2 pipelines of diameter size d for local distribution in region g in time period t'
AYoff(d,g,r,t)   'availability of offshore CO2 pipelines of diameter size d for local distribution in region g in time period t'
AYst(d,g,sc,t)    'availability of hydrogen pipelines of diameter size d for distribution in region g in time period t'
;

integer Variables
Yh(d,g,g1,t)     'establishment of hydrogen pipelines of diameter size d for regioanal distribution in region g in time period t'
Yon(d,g,g1,t)    'establishment of onshore CO2pipelines of diameter size d in region g in time period t'
Yoff(d,g,r,t)    'establishment of offshore CO2 pipelines of diameter size d in region g in time period t'
Yst(d,g,sc,t)     'establishment of hydrogen pipelines of diameter size d in region g in time period d to storage type s'
;

positive Variables
BC(it)               'Biomass Cost (?)'
CL(it,g,t,c,h)      'curtailment in region g, time period t, cluster c, hour h (MW)'
InvR(e,g,t)      'invested capacity of renewable e in region g and time period t (MW)'
FCR(it)             'Fuel cost for regional transport (?)'
GC(it)               'Gas Cost (?)'
GCR(it)              'General Cost for regional transport (?)'
IIC(it)            'International import cost (?)'
IMP(it,g,t,c,h)     'flow rate of international import in region g in time period t (MW)'
LCR(it)              'Labour cost for regional transport (?)'
MCR(it)              'Maintenance cost for regional transport (?)'
NR(e,g,t)        'capacity of renewable e in region g and time period t (MW)'
Pr(it,p,g,t,c,h)    'Production rate of product type i produced by a plant of type j and size p in region g in time period t (MW)'
Pre(it,e,g,t,c,h)   'electricity production from renewable e in region g, time period t, cluster c, hour h (MW)'
PipeCC           'Pipeline capital cost (?k)'
PipeOC           'Pipeline operating cost (?k)'
PCC              'Production capital cost (?k)'
POC1              'Production operating cost (?)'
POC2(it)
Q(it,l,g,g1,t,c,h)  'regional flowrate of H2 via transportation mode l in region g in time period t (MWh)'
Qi(it,g,s,t,c,h)    'flowrate of H2 via pipeline from region g to storage type s in time period t(MWh)'
Qr(it,s,g,t,c,h)    'flowrate of H2 via pipeline from region g to storage type s in time period t(MWh)'
Qon(it,g,g1,t,c,h)  'regional flowrate of CO2 via onshore pipelines between regions g and g?? in time period t (kg CO2/d)'
Qoff(it,g,r,t,c,h)  'flowrate of CO2 via oofshore pipelines from a collection point in region g to a reservoir r in time period t (kg CO2/d)'
RCC              'Road transportation capital cost (?)'
Rdown(p,g,t,c,h) 'upward reserve contribution (MWh)'
ReC              'renewables capital cost (?)'
RI(r,t)          'Inventory of CO2 in reservoir r in time period t (kg CO2-eq)'
ROC(it)              'Road transportation cost (?)'
Rup(p,g,t,c,h)   'downward reserve contribution (MWh)'
SCC              'Storage capical cost (?)'
SOC
POC
SOC1             'Storage operating cost (?)'
SOC2(it)
St(it,s,g,t,c,h)    'Average inventory of product type i stored in a storage facility of type s and size p in region g in time period t (kW)'
TCC              'Transportation capital cost (?)'
TOC              'Transportation operating cost (?)'
Vbio(t)          'Biomass consumption in time period t, cluster c and hour h (kg)'
Vgas(t)          'Gas consumpiton in timer period t, cluster c and hour h ('

;
variable
Stage1_TC
Stage1_eta;

variables
Stage2_TC 
TC               'Total cost (?)'
*em(t)            'total emissions (tCO2)'
CEC(it)              'Carbon emissions cost (?)'
;

Parameter
Stage1_dem(it,g,t,c,h);

RU(p)$(ord(p)<=3)=0.1;
RD(p)$(ord(p)<=3)=0.1;
display  pcap_min,RU,RD;


Equations
Stage1_Obj,Stage2_Obj,POCeq, SOCeq, PCCeq,POCeq1,POCeq2(it),SCCeq,SOCeq1,SOCeq2(it),TCCeq,RCCeq,PipeCCeq,TOCeq,PipeOCeq,
IICeq(it),ReCeq,co2MassBalance(it,g,t,c,h),
PCapacity1(it,p,g,t,c,h),PCapacity2(it,p,g,t,c,h),
PAvailability(p,g,t),SAvailability(s,g,t),
MaxInj(it,t,g,s,c,h),MaxRetr(it,t,g,s,c,h),
SCapacityU(s,g,t,c,h),SCapacity1(it,s,g,t,c,h),SCapacity2(it,s,g,t,c,h),
TAvailability(g,g1,t),
H2PipeMax(it,g,g1,t,c,h),OnshorePipeMax(it,g,g1,t,c,h),OffshorePipeMax(it,g,r,t,c,h),PipeStAvailability(d,g,s,t),
H2PAvailability(d,g,g1,t),OnPAvailability(d,g,g1,t),OffPAvailability(d,g,r,t),
H2Pipe(g,g1,t),OnPipe(g,g1,t),OffPipe(g,r,t),StPipe(s,g,t),
ResInventory(it,r,t),
ImpLimit(it,t,c,h),
ElecProd(it,g,t,c,h),RenewAv(it,e,g,t,c,h),RenewCap(e,g,t),
CurtPercentage(it,c,h),
Emissions(t),
GasCons(t),
BioCons(t),
LandAvailability(e,g,t),
BiomassAvailability(it,g,t),
Stage1_auxislack(it),
POCeq2(it),ROCeq(it),FCReq(it),GCReq(it),LCReq(it),MCReq(it),CECeq(it),GasCost(it),BioCost(it)
h2MassBalance(it,g,t,c,h),
EmTargeteq(it,t),SInventory2(it,s,g,t,c,h),
SFinal(it,s,g,t,c),
RampUp(it,p,g,c,h,t),RampDown(it,p,g,c,h,t)
;

*%%%%%%% The initial objective function with all continous variables
Stage1_Obj..
Stage1_TC  =e= 1000*PCC + SCC + 1000*PipeCC + SOC1 + POC1 + 1000*PipeOC  + ReC + Stage1_eta;

*facilities capital cost
PCCeq..
PCC =e= sum((p,g,t)$TT(t),dfc(t)*pccost(p,t)*PCap(p)*InvP(p,g,t));

SCCeq..
SCC =e= sum((s,g,t)$(TT(t) and GS(g,s)),dfc(t)*sccost(s)*SCap(s)*InvS(s,g,t));

*pipeline transportation capital cost
PipeCCeq..
PipeCC=e= sum((t,g,g1,d)$(Npipe(g,g1) and TT(t) and (ord(g)<ord(g1)) and (ord(d) = 3)),    dfc(t)*cccH(d)* DistPipe(g,g1)*Yh(d,g,g1,t))+
          sum((t,g,g1,d)$(N(g,g1)     and TT(t) and (ord(g)<ord(g1)) and (ord(d) =2)), dfc(t)*cccC_onshore(d)* Dist(g,g1)*    Yon(d,g,g1,t))+
          sum((t,g,r,d)$(GR(g,r)      and TT(t) and (ord(d) =2)),                      dfc(t)*cccC_offshore(d)*DistRes(g,r)*  Yoff(d,g,r,t))+
          sum((t,g,sc,d)$(GS(g,sc)and TT(t) and (ord(d) = 3)),                                    dfc(t)*cccH(d)*         DistSt(g,sc)*  Yst(d,g,sc,t));

*pipline operating cost
PipeOCeq..
PipeOC=e=sum((t,g,g1,d)$(Npipe(g,g1) and TT(t) and (ord(g)<ord(g1)) and (ord(d) = 3)),    dfo(t)*deltaH* crf*cccH(d)* DistPipe(g,g1)*AY(d,g,g1,t))
      +  sum((t,g,g1,d)$(N(g,g1)     and TT(t) and (ord(g)<ord(g1)) and (ord(d) = 2)),    dfo(t)*deltaC_onshore* crf*cccC_onshore(d)* Dist(g,g1)*    AYon(d,g,g1,t))
      +  sum((t,g,d,r)$(GR(g,r)      and TT(t) and (ord(d) = 2)) ,                        dfo(t)*deltaC_offshore*crf*cccC_offshore(d)*DistRes(g,r)*  AYoff(d,g,r,t))
      +  sum((t,g,sc,d)$(GS(g,sc) and TT(t) and (ord(d) = 3))  ,              dfo(t)*deltaH*         crf*cccH(d)*         DistSt(g,sc)*  AYst(d,g,sc,t));


*facilities operation cost
POCeq1..
POC1  =e= sum((p,g,t)$TT(t),dfo(t)*pocostF(p,t)*PCap(p)*NP(p,g,t));

SOCeq1..
SOC1  =e= sum((s,g,t)$(TT(t) and GS(g,s)),dfo(t)*socostF(s)*SCap(s)*NS(s,g,t));

*renewables total cost
ReCeq..
ReC =E= sum((e,g,t)$TT(t), (dfc(t)*rccost(e,t)*InvR(e,g,t)+ dfo(t)*rocost(e,t)*NR(e,g,t)));

Stage1_auxislack(iter_fesi)..
Stage1_eta =g= POC2(iter_fesi) + SOC2(iter_fesi)  + CEC(iter_fesi) + IIC(iter_fesi) + GC(iter_fesi) + BC(iter_fesi);

POCeq2(iter_fesi)..
POC2(iter_fesi)  =e= sum((p,g,t,c,h)$(TT(t) and CC(c) and HH(h)),dfo(t)*WF(c)*pocostV(p,t)*theta*Pr(iter_fesi,p,g,t,c,h));

SOCeq2(iter_fesi)..
SOC2(iter_fesi)  =e= sum((s,g,t,c,h)$(TT(t) and GS(g,s) and CC(c) and HH(h)),dfo(t)*WF(c)*socostV(s)*theta*Qi(iter_fesi,g,s,t,c,h) );

*carbon emissions cost
CECeq(iter_fesi)..
CEC(iter_fesi) =e= sum((p,g,t,c,h)$(TT(t) and CC(c) and HH(h)),WF(c)*dfo(t)*ct(t)*y_e(p,t)*theta*Pr(iter_fesi,p,g,t,c,h));

*international import cost
IICeq(iter_fesi)..
IIC(iter_fesi) =e= sum((t,g,c,h)$(gimp(g)and TT(t) and CC(c) and HH(h)),WF(c)*dfo(t)*pimp*theta*IMP(iter_fesi,g,t,c,h));

*Fuels cost
GasCost(iter_fesi)..
GC(iter_fesi) =e= sum((g,p,c,h,t)$(TT(t) and ord(p)<=2 and CC(c) and HH(h)), dfo(t)*cgas(t)*WF(c)*theta*Pr(iter_fesi,p,g,t,c,h)/eta(p,t));

BioCost(iter_fesi)..
BC(iter_fesi) =e= sum((g,c,h,t)$(TT(t) and CC(c) and HH(h)), dfo(t)*cbio(t)*WF(c)*theta*Pr(iter_fesi,'BECCS',g,t,c,h)/eta('BECCS',t));

*==============Constraints===========================================

*=Hygrogen MASS BALANCES=
h2MassBalance(iter,g,t,c,h)$(TT(t) and CC(c) and HH(h))..
sum((p),Pr(iter,p,g,t,c,h)) + sum(g1$Npipe(g1,g), Q(iter,'Pipe',g1,g,t,c,h))+ IMP(iter,g,t,c,h)$(gimp(g))+ sum(s$(GS(g,s)),Qr(iter,s,g,t,c,h))
=g=
sum(g1$Npipe(g,g1), Q(iter,'Pipe',g,g1,t,c,h))+ sum(s$(GS(g,s)),Qi(iter,g,s,t,c,h))+ Stage1_dem(iter,g,t,c,h);

*=co2 MASS BALANCES=
co2MassBalance(iter,g,t,c,h)$(TT(t) and CC(c) and HH(h))..
sum(g1$N(g1,g),Qon(iter,g1,g,t,c,h))+sum(p,y_c(p,t)*Pr(iter,p,g,t,c,h)) =e=
sum(g1$N(g,g1),Qon(iter,g,g1,t,c,h))+sum(r$GR(g,r),Qoff(iter,g,r,t,c,h));

*Biomass availability
BiomassAvailability(iter,g,t)$(TT(t))..
sum((c,h)$(CC(c) and HH(h)),WF(c)*theta*Pr(iter,'BECCS',g,t,c,h)/eta('BECCS',t)) =l= BA(g,t);

*======================================RAMP UP/DOWN========================================================*
*Ramp Up
RampUp(iter,p,g,c,h,t)$(TT(t) and CC(c) and HH(h) and ord(h)>1)..
Pr(iter,p,g,t,c,h) - Pr(iter,p,g,t,c,h-1)  =l= theta*RU(p)*PCap(p)*NP(p,g,t);

*Ramp Down
RampDown(iter,p,g,c,h,t)$(TT(t) and CC(c) and HH(h) and ord(h)>1)..
Pr(iter,p,g,t,c,h-1)-Pr(iter,p,g,t,c,h)  =l= theta*RD(p)*PCap(p)*NP(p,g,t);

*======================================PRODUCTION CONSTRAINTS==============================================*
PCapacity1(iter,p,g,t,c,h)$(TT(t) and CC(c) and HH(h))..
Pr(iter,p,g,t,c,h)=g=PCap(p)*pcap_min(p)*NP(p,g,t)  ;

PCapacity2(iter,p,g,t,c,h)$(TT(t) and CC(c) and HH(h))..
Pr(iter,p,g,t,c,h)=l=PCap(p)*pcap_max(p)*NP(p,g,t)  ;

*availability in production plants
PAvailability(p,g,t)$(TT(t))..
InvP(p,g,t)  =e= NP(p,g,t) - NP(p,g,t-1)$(ord(t)>y1) - np0(p,g)$(ord(t)=y1);

SInventory2(iter,s,g,t,c,h)$(GS(g,s)and TT(t) and CC(c) and HH(h))..
St(iter,s,g,t,c,h) =e= St(iter,s,g,t,c,h-1)$(ord(h)>1) + st0(s,g)$(ord(h)=1) +theta*(Qi(iter,g,s,t,c,h) - Qr(iter,s,g,t,c,h)) ;

*maximum injection and retrieval rate
MaxInj(iter,t,g,s,c,h)$(GS(g,s) and TT(t) and CC(c) and HH(h))..
Qi(iter,g,s,t,c,h)=l=QImax(s)*NS(s,g,t);

MaxRetr(iter,t,g,s,c,h)$(GS(g,s) and TT(t) and CC(c) and HH(h))..
Qr(iter,s,g,t,c,h)=l=QRmax(s)*NS(s,g,t);

*storage capacity constraints for underground storage
SCapacityU(sc,g,t,c,h)$(GS(g,sc) and TT(t) and CC(c) and HH(h))..
InvS(sc,g,t)=l=sum(d,Yst(d,g,sc,t));

*storage capacity constraints
SCapacity1(iter,s,g,t,c,h)$(GS(g,s) and TT(t) and CC(c) and HH(h))..
St(iter,s,g,t,c,h)=g=SCap(s)*scap_min(s)*NS(s,g,t)  ;

SCapacity2(iter,s,g,t,c,h)$(GS(g,s) and TT(t) and CC(c) and HH(h))..
St(iter,s,g,t,c,h)=l=SCap(s)*scap_max(s)*NS(s,g,t) ;

*availability of storage facility
SAvailability(s,g,t)$(GS(g,s) and TT(t))..
NS(s,g,t) =e=  NS(s,g,t-1)$(ord(t)>y1) +  ns0(s,g)$(ord(t)=y1) +InvS(s,g,t);

SFinal(iter,s,g,t,c)$(GS(g,s) and TT(t) and CC(c))..
St(iter,s,g,t,c,'24') =e= 0;

*======================================RENEWABLES CONSTRAINTS==============================================*
*Electricity production for electrolysis
ElecProd(iter,g,t,c,h)$(TT(t) and CC(c) and HH(h))..
Pr(iter,'WE',g,t,c,h) =e= eta('WE',t)*(sum(e,Pre(iter,e,g,t,c,h)) - CL(iter,g,t,c,h));

*Renewables Availability
RenewAv(iter,e,g,t,c,h)$(TT(t) and CC(c) and HH(h))..
Pre(iter,e,g,t,c,h) =e= 0.7*AV(c,h,g,e)*NR(e,g,t);
******************HA1*************************************************
*Pre(e,g,t,c,h) =e= 0.7*AVH(c,h,g,e)*NR(e,g,t);

*Renewables capacity expansion
RenewCap(e,g,t)$(TT(t))..
NR(e,g,t) =e= NR(e,g,t-1)$(ord(t)>y1) + InvR(e,g,t);

*Land availability
LandAvailability(e,g,t)$(TT(t))..
NR(e,g,t) =l= landAV(e,g);

*Curtailment limit
CurtPercentage(iter,c,h)$(CC(c) and HH(h))..
sum((g,t),CL(iter,g,t,c,h)) =l= 0.1*sum((e,g,t), Pre(iter,e,g,t,c,h)) ;

*=====================================PIPELINE CONSTRAINTS=====================================================*
H2PipeMax(iter,g,g1,t,c,h)$(Npipe(g,g1) and TT(t) and CC(c) and HH(h))..
Q(iter,'Pipe',g,g1,t,c,h)=l= sum(d$(ord(d) =3), qHmax(d)*(AY(d,g,g1,t)$(ord(g)<ord(g1))+AY(d,g1,g,t)$(ord(g1)<ord(g))));

OnshorePipeMax(iter,g,g1,t,c,h)$(N(g,g1) and TT(t) and CC(c) and HH(h))..
*Qon(g,g1,t,c,h)=l= sum(d,10*qCmax(d)*(AYon(d,g,g1,t)$(ord(g)<ord(g1))+AYon(d,g1,g,t)$(ord(g1)<ord(g))));
Qon(iter,g,g1,t,c,h)=l= sum(d$(ord(d)=2),qCmax(d)*(AYon(d,g,g1,t)$(ord(g)<ord(g1))+AYon(d,g1,g,t)$(ord(g1)<ord(g))));

OffshorePipeMax(iter,g,r,t,c,h)$(GR(g,r) and TT(t) and CC(c) and HH(h))..
*Qoff(g,r,t,c,h)=l= sum(d,10*qCmax(d)*AYoff(d,g,r,t));
Qoff(iter,g,r,t,c,h)=l= sum(d$(ord(d)=2),qCmax(d)*AYoff(d,g,r,t));

*availability of pipelines
H2PAvailability(d,g,g1,t)$(Npipe(g,g1) and (ord(g)<ord(g1)) and TT(t)  and ord(d) =3)..
*AY(d,g,g1,t)   =e= AY(d,g,g1,t-1)$(ord(t)>y1)  +  ayHR0(d,g,g1)$(ord(t)=y1) + Yh(d,g,g1,t) - Yh(d,g,g1,t-(LTpipe/dur))$(ord(t)-(LTpipe/dur) gt y1)  ;
AY(d,g,g1,t)   =e= AY(d,g,g1,t-1)$(ord(t)>y1)  +  ayHR0(d,g,g1)$(ord(t)=y1) + Yh(d,g,g1,t) ;

OnPAvailability(d,g,g1,t)$(N(g,g1) and (ord(g)<ord(g1)) and TT(t) and ord(d)=2)..
*AYon(d,g,g1,t) =e= AYon(d,g,g1,t-1)$(ord(t)>y1) +  ayC0(d,g,g1)$(ord(t)=y1) + Yon(d,g,g1,t) - Yon(d,g,g1,t-(LTonshore/dur))$((ord(t)-(LTonshore/dur)) gt y1)  ;
AYon(d,g,g1,t) =e= AYon(d,g,g1,t-1)$(ord(t)>y1) +  ayC0(d,g,g1)$(ord(t)=y1) + Yon(d,g,g1,t) ;

OffPAvailability(d,g,r,t)$(GR(g,r) and TT(t) and ord(d)=2)..
*AYoff(d,g,r,t) =e= AYoff(d,g,r,t-1)$(ord(t)>y1) +  aeC0(r)$(ord(t)=1) + Yoff(d,g,r,t) - Yoff(d,g,r,t-(LToffshore/dur))$((ord(t)-(LToffshore/dur)) gt y1)   ;
AYoff(d,g,r,t) =e= AYoff(d,g,r,t-1)$(ord(t)>y1) +  aeC0(r)$(ord(t)=y1) + Yoff(d,g,r,t) ;

PipeStAvailability(d,g,sc,t)$(GS(g,sc) and TT(t)  and ord(d) =3)..
*AYst(d,g,sc,t) =e= AYst(d,g,sc,t-1)$(ord(t)>y1) +  0$(ord(t)=1) + Yst(d,g,sc,t) - Yst(d,g,sc,t-(LTpipe/dur))$((ord(t)-(LTpipe/dur)) gt y1)   ;
AYst(d,g,sc,t) =e= AYst(d,g,sc,t-1)$(ord(t)>y1) + Yst(d,g,sc,t) ;

*======================================RESERVOIRS==============================================================*
*inventory
ResInventory(iter,r,t)$TT(t)..
RI(r,t)=e=RI(r,t-1)$(ord(t)>y1) +ri0(r)$(ord(t)=y1)/1000 + dur*sum((g,c,h)$GR(g,r),WF(c)*theta*Qoff(iter,g,r,t,c,h))$(ord(t)>=y1)/1000;

*inventory level maximum
*ResCapacity(r,t)$TT(t)..                      RI(r,t)=l=sum((d,g)$GR(g,r),AYoff(d,g,r,t)*rcap(r));

*======================================IMPORT=================================================================*
*Import limit
ImpLimit(iter,t,c,h)$(TT(t) and CC(c) and HH(h))..
sum(g$Gimp(g),Imp(iter,g,t,c,h))=l=iota*sum(g,Gasdemave(c,h,t,g));

*EMISSIONS==
EmTargeteq(iter,t)$TT(t)..
sum((p,g,c,h)$(CC(c) and HH(h)),WF(c)*y_e(p,t)*theta*Pr(iter,p,g,t,c,h)) =l= emtarget(t);


Qon.up(it,g,g1,t,c,h)$(TT(t) and CC(c) and HH(h))=1.17E+04 ;
Qoff.up(it,g,r,t,c,h)$(TT(t) and CC(c) and HH(h))=1.17E+04 ;


*ITU.up(l,g,g1,t)$TT(t)=25;
Q.up(it,'Pipe',g,g1,t,c,h)$(TT(t) and CC(c) and HH(h) and Npipe(g,g1))=15343;


model Stage1_Problem /
Stage1_Obj,PCCeq,SCCeq,PipeCCeq,PipeOCeq, POCeq1, SOCeq1, ReCeq, Stage1_auxislack,
POCeq2, SOCeq2,CECeq,IICeq,GasCost,BioCost,h2MassBalance,co2MassBalance,BiomassAvailability,RampUp,RampDown,
PCapacity2,PAvailability,
SInventory2,MaxInj,MaxRetr,SCapacityU,SCapacity1,SCapacity2,SAvailability,
SFinal,ElecProd,RenewAv,RenewCap,LandAvailability,
*CurtPercentage,
H2PipeMax,OnshorePipeMax,OffshorePipeMax
PipeStAvailability,
H2PAvailability,OnPAvailability,OffPAvailability,ResInventory,ImpLimit,EmTargeteq
/;

model Stage1_Problem_HA2 /
Stage1_Obj,PCCeq,SCCeq,PipeCCeq,PipeOCeq, POCeq1, SOCeq1, ReCeq, Stage1_auxislack,
POCeq2, SOCeq2,CECeq,IICeq,GasCost,BioCost,h2MassBalance,co2MassBalance,BiomassAvailability,RampUp,RampDown,
PCapacity2,PAvailability,
SInventory2,MaxInj,MaxRetr,SCapacityU,SCapacity1,SCapacity2,SAvailability,
SFinal,ElecProd,RenewAv,RenewCap,LandAvailability,
*CurtPercentage,
*H2PipeMax,OnshorePipeMax,OffshorePipeMax
PipeStAvailability,
H2PAvailability,OnPAvailability,OffPAvailability,ResInventory,ImpLimit,EmTargeteq
/;


Stage1_Problem.optfile=1;
Stage1_Problem_HA2.optfile=1;

SCCeq.scale=1000;
PCCeq.scale=1000;
POCeq1.scale=1000;
SOCeq1.scale=1000;
POCeq2.scale(it)=1000;
SOCeq2.scale(it)=1000;
EmTargeteq.scale(it,t)=10000;

Stage1_Problem.scaleopt=1;
Stage1_Problem_HA2.scaleopt=1;

*============solve maxmin problem====================
Equations
Stage2_middle_Obj,
Stage2_lower_Obj,POCeq2_Stage2,SOCeq2_Stage2,IICeq_Stage2,CECeq_Stage2,GasCost_Stage2,BioCost_Stage2,
h2MassBalance_Stage2(g,t,c,h),co2MassBalance_Stage2(g,t,c,h),BiomassAvailability_Stage2(g,t),
RampUp_Stage2(p,g,c,h,t),RampDown_Stage2(p,g,c,h,t),
PCapacity2_Stage2(p,g,t,c,h),SInventory2_Stage2(s,g,t,c,h),
MaxInj_Stage2(t,g,s,c,h),MaxRetr_Stage2(t,g,s,c,h),SCapacity1_Stage2(s,g,t,c,h),SCapacity2_Stage2(s,g,t,c,h),SFinal_Stage2(s,g,t,c),
ElecProd_Stage2(g,t,c,h),RenewAv_Stage2(e,g,t,c,h),
CurtPercentage_Stage2(c,h),
H2PipeMax_Stage2(g,g1,t,c,h),
OnshorePipeMax_Stage2(g,g1,t,c,h),OffshorePipeMax_Stage2(g,r,t,c,h),
ResInventory_Stage2(r,t),ImpLimit_Stage2(t,c,h),Demand_Stage2_cons(g,t,c,h),EmTargeteq_Stage2(t)
;

positive Variables
Pr_Stage2(p,g,t,c,h),Q_Stage2(l,g,g1,t,c,h),
Qi_Stage2(g,s,t,c,h),Qr_Stage2(s,g,t,c,h),
Qon_Stage2(g1,g,t,c,h),Qoff_Stage2(g,r,t,c,h),
IMP_Stage2(g,t,c,h),St_Stage2(s,g,t,c,h),CL_Stage2(g,t,c,h),Pre_Stage2(e,g,t,c,h),
POC2_Stage2,SOC2_Stage2,IIC_Stage2,GC_Stage2,BC_Stage2
;

variables
CEC_Stage2,
Stage2_dem(g,t,c,h),
Stage2_middle_dem(g,t,c,h)
Stage2_lower_TC,
Stage2_middle_TC;

Parameters
Stage2_fixed_dem(g,t,c,h),
Lem_stage2(g,t,c,h);

Stage2_lower_Obj..
Stage2_lower_TC =e= POC2_Stage2 + SOC2_Stage2 + CEC_Stage2 + IIC_Stage2 + GC_Stage2 + BC_Stage2;

POCeq2_Stage2..
POC2_Stage2  =e= sum((p,g,t,c,h)$(TT(t) and CC(c) and HH(h)),dfo(t)*WF(c)*pocostV(p,t)*theta*Pr_Stage2(p,g,t,c,h));

SOCeq2_Stage2..
SOC2_Stage2  =e= sum((s,g,t,c,h)$(TT(t) and GS(g,s) and CC(c) and HH(h)),dfo(t)*WF(c)*socostV(s)*theta*Qi_Stage2(g,s,t,c,h) );

*international import cost
IICeq_Stage2..
IIC_Stage2 =e= sum((t,g,c,h)$(gimp(g)and TT(t) and CC(c) and HH(h)),WF(c)*dfo(t)*pimp*theta*IMP_Stage2(g,t,c,h));

*carbon emissions cost
CECeq_Stage2..
CEC_Stage2 =e= sum((p,g,t,c,h)$(TT(t) and CC(c) and HH(h)),WF(c)*dfo(t)*ct(t)*y_e(p,t)*theta*Pr_Stage2(p,g,t,c,h));

*Fuels cost
GasCost_Stage2..
GC_Stage2 =e= sum((g,p,c,h,t)$(TT(t) and ord(p)<=2 and CC(c) and HH(h)), dfo(t)*cgas(t)*WF(c)*theta*Pr_Stage2(p,g,t,c,h)/eta(p,t));

BioCost_Stage2..
BC_Stage2 =e= sum((g,c,h,t)$(TT(t) and CC(c) and HH(h)), dfo(t)*cbio(t)*WF(c)*theta*Pr_Stage2('BECCS',g,t,c,h)/eta('BECCS',t));

h2MassBalance_Stage2(g,t,c,h)$(TT(t) and CC(c) and HH(h))..
sum((p),Pr_Stage2(p,g,t,c,h)) + sum(g1$Npipe(g1,g), Q_Stage2('Pipe',g1,g,t,c,h))+ IMP_Stage2(g,t,c,h)$(gimp(g))+ sum(s$(GS(g,s)),Qr_Stage2(s,g,t,c,h))
=g=
sum(g1$Npipe(g,g1), Q_Stage2('Pipe',g,g1,t,c,h))+ sum(s$(GS(g,s)),Qi_Stage2(g,s,t,c,h))+ Stage2_dem(g,t,c,h);

*=co2 MASS BALANCES=
co2MassBalance_Stage2(g,t,c,h)$(TT(t) and CC(c) and HH(h))..
sum(g1$N(g1,g),Qon_Stage2(g1,g,t,c,h))+sum(p,y_c(p,t)*Pr_Stage2(p,g,t,c,h)) =e=
sum(g1$N(g,g1),Qon_Stage2(g,g1,t,c,h))+sum(r$GR(g,r),Qoff_Stage2(g,r,t,c,h));

*Biomass availability
BiomassAvailability_Stage2(g,t)$(TT(t))..
sum((c,h)$(CC(c) and HH(h)),WF(c)*theta*Pr_Stage2('BECCS',g,t,c,h)/eta('BECCS',t)) =l= BA(g,t);

*======================================RAMP UP/DOWN========================================================*
*Ramp Up
RampUp_Stage2(p,g,c,h,t)$(TT(t) and CC(c) and HH(h) and ord(h)>1)..
Pr_Stage2(p,g,t,c,h) - Pr_Stage2(p,g,t,c,h-1)  =l= theta*RU(p)*PCap(p)*NP.l(p,g,t);

*Ramp Down
RampDown_Stage2(p,g,c,h,t)$(TT(t) and CC(c) and HH(h) and ord(h)>1)..
Pr_Stage2(p,g,t,c,h-1)-Pr_Stage2(p,g,t,c,h)  =l= theta*RD(p)*PCap(p)*NP.l(p,g,t);

*======================================PRODUCTION CONSTRAINTS==============================================*

PCapacity2_Stage2(p,g,t,c,h)$(TT(t) and CC(c) and HH(h))..
Pr_Stage2(p,g,t,c,h)=l=PCap(p)*pcap_max(p)*NP.l(p,g,t) ;

SInventory2_Stage2(s,g,t,c,h)$(GS(g,s)and TT(t) and CC(c) and HH(h))..
St_Stage2(s,g,t,c,h) =e= St_Stage2(s,g,t,c,h-1)$(ord(h)>1) + st0(s,g)$(ord(h)=1) +theta*(Qi_Stage2(g,s,t,c,h) - Qr_Stage2(s,g,t,c,h)) ;

*maximum injection and retrieval rate
MaxInj_Stage2(t,g,s,c,h)$(GS(g,s) and TT(t) and CC(c) and HH(h))..
Qi_Stage2(g,s,t,c,h)=l=QImax(s)*NS.l(s,g,t);

MaxRetr_Stage2(t,g,s,c,h)$(GS(g,s) and TT(t) and CC(c) and HH(h))..
Qr_Stage2(s,g,t,c,h)=l=QRmax(s)*NS.l(s,g,t);

*storage capacity constraints
SCapacity1_Stage2(s,g,t,c,h)$(GS(g,s) and TT(t) and CC(c) and HH(h))..
St_Stage2(s,g,t,c,h)=g=SCap(s)*scap_min(s)*NS.l(s,g,t)  ;

SCapacity2_Stage2(s,g,t,c,h)$(GS(g,s) and TT(t) and CC(c) and HH(h))..
St_Stage2(s,g,t,c,h)=l=SCap(s)*scap_max(s)*NS.l(s,g,t) ;

SFinal_Stage2(s,g,t,c)$(GS(g,s) and TT(t) and CC(c))..
St_Stage2(s,g,t,c,'24') =e= 0;

*======================================RENEWABLES CONSTRAINTS==============================================*
*Electricity production for electrolysis
ElecProd_Stage2(g,t,c,h)$(TT(t) and CC(c) and HH(h))..
Pr_Stage2('WE',g,t,c,h) =e= eta('WE',t)*(sum(e,Pre_Stage2(e,g,t,c,h)) - CL_Stage2(g,t,c,h));

*Renewables Availability
RenewAv_Stage2(e,g,t,c,h)$(TT(t) and CC(c) and HH(h))..
Pre_Stage2(e,g,t,c,h) =e= 0.7*AV(c,h,g,e)*NR.l(e,g,t);
******************HA1*************************************************
*Pre(e,g,t,c,h) =e= 0.7*AVH(c,h,g,e)*NR(e,g,t);

*Curtailment limit
CurtPercentage_Stage2(c,h)$(CC(c) and HH(h))..
sum((g,t),CL_Stage2(g,t,c,h)) =l= 0.1*sum((e,g,t), Pre_Stage2(e,g,t,c,h));

*=====================================PIPELINE CONSTRAINTS=====================================================*
H2PipeMax_Stage2(g,g1,t,c,h)$(Npipe(g,g1) and TT(t) and CC(c) and HH(h))..
Q_Stage2('Pipe',g,g1,t,c,h)=l= sum(d$(ord(d) =3), qHmax(d)*(AY.l(d,g,g1,t)$(ord(g)<ord(g1))+AY.l(d,g1,g,t)$(ord(g1)<ord(g))));

OnshorePipeMax_Stage2(g,g1,t,c,h)$(N(g,g1) and TT(t) and CC(c) and HH(h))..
*Qon(g,g1,t,c,h)=l= sum(d,10*qCmax(d)*(AYon(d,g,g1,t)$(ord(g)<ord(g1))+AYon(d,g1,g,t)$(ord(g1)<ord(g))));
Qon_Stage2(g,g1,t,c,h)=l= sum(d$(ord(d)=2),qCmax(d)*(AYon.l(d,g,g1,t)$(ord(g)<ord(g1))+AYon.l(d,g1,g,t)$(ord(g1)<ord(g))));

OffshorePipeMax_Stage2(g,r,t,c,h)$(GR(g,r) and TT(t) and CC(c) and HH(h))..
*Qoff(g,r,t,c,h)=l= sum(d,10*qCmax(d)*AYoff(d,g,r,t));
Qoff_Stage2(g,r,t,c,h)=l= sum(d$(ord(d)=2),qCmax(d)*AYoff.l(d,g,r,t));

ResInventory_Stage2(r,t)$TT(t)..
RI.l(r,t)=e=RI.l(r,t-1)$(ord(t)>y1) +ri0(r)$(ord(t)=y1)/1000 + dur*sum((g,c,h)$GR(g,r),WF(c)*theta*Qoff_Stage2(g,r,t,c,h))$(ord(t)>=y1)/1000;

*======================================IMPORT=================================================================*
*Import limit
ImpLimit_Stage2(t,c,h)$(TT(t) and CC(c) and HH(h))..
sum(g$Gimp(g),Imp_Stage2(g,t,c,h))=l=iota*sum(g,Gasdemave(c,h,t,g));

EmTargeteq_Stage2(t)$TT(t)..
sum((p,g,c,h)$(CC(c) and HH(h)),WF(c)*y_e(p,t)*theta*Pr_Stage2(p,g,t,c,h)) =l= emtarget(t) ;

Demand_Stage2_cons(g,t,c,h)$(TT(t) and CC(c) and HH(h))..
Stage2_dem(g,t,c,h) =e= Stage2_fixed_dem(g,t,c,h);

model Stage2_lower_problem/
Stage2_lower_Obj,POCeq2_Stage2,SOCeq2_Stage2,IICeq_Stage2,CECeq_Stage2,GasCost_Stage2,BioCost_Stage2,
h2MassBalance_Stage2,co2MassBalance_Stage2,BiomassAvailability_Stage2,
RampUp_Stage2,RampDown_Stage2,PCapacity2_Stage2,SInventory2_Stage2,
MaxInj_Stage2,MaxRetr_Stage2,SCapacity1_Stage2,SCapacity2_Stage2,SFinal_Stage2,
ElecProd_Stage2,RenewAv_Stage2,
*CurtPercentage_Stage2,
H2PipeMax_Stage2,OnshorePipeMax_Stage2,OffshorePipeMax_Stage2,
ResInventory_Stage2,ImpLimit_Stage2,Demand_Stage2_cons,EmTargeteq_Stage2
/;


Parameters
Stage2_phi          'uncertainty budget'
;
Stage2_phi = 2;

*Binary
positive variables
Z_neg(t,g,c,j)  'backward deviation vector'
Z_pos(t,g,c,j)  'forward deviation vector'
;

Equations
Stage2_Middle_demand_const1(g,t,c,h)
Stage2_Middle_demand_const2(g,t,c,h)
Constr_Z1(t,g,c,j)
Constr_Z2(t,g,c)
Constr_Z3(t,g,c,j)
Constr_Z4(t,g,c,j)
;

Stage2_middle_Obj..
Stage2_middle_TC =e= Stage2_lower_TC.l
+ sum((g,t,c,h)$(TT(t) and CC(c) and HH(h)), Lem_stage2(g,t,c,h)*(Stage2_middle_dem(g,t,c,h) -Stage2_fixed_dem(g,t,c,h)));

Stage2_Middle_demand_const1(g,t,c,h)$(TT(t) and ord(c)=1 and HH(h))..
Stage2_middle_dem(g,t,c,h) =e= dem(g,t,c,h);

Stage2_Middle_demand_const2(g,t,c,h)$(TT(t) and CC(c) and ord(c)>1 and HH(h))..
Stage2_middle_dem(g,t,c,h) =e= Gasdemave(c,h,t,g) + sum(j, BigQ_neg(c,h,j,t,g)*Z_neg(t,g,c,j) + BigQ_pos(c,h,j,t,g)*Z_pos(t,g,c,j));
*dem(g,t,c,h)=round(dc(t)*GasDem(c,h,g),2);

*==================== Constr Z================

Constr_Z1(t,g,c,j)$(CC(c) and ord(c)>1)..
Z_neg(t,g,c,j) + Z_pos(t,g,c,j) =l= 1;

Constr_Z2(t,g,c)$(CC(c) and ord(c)>1)..
sum(j, Z_neg(t,g,c,j) + Z_pos(t,g,c,j)) =l= Stage2_phi;

Constr_Z3(t,g,c,j)$(CC(c) and ord(c)>1)..
Z_neg(t,g,c,j) =l= 1;

Constr_Z4(t,g,c,j)$(CC(c) and ord(c)>1)..
Z_pos(t,g,c,j) =l= 1;

model Stage2_middle_problem/
Stage2_middle_Obj,Stage2_Middle_demand_const1,Stage2_Middle_demand_const2,Constr_Z1,Constr_Z2,Constr_Z3,Constr_Z4
/;

Qon_stage2.up(g,g1,t,c,h)$(TT(t) and CC(c) and HH(h))=1.17E+04 ;
Qoff_stage2.up(g,r,t,c,h)$(TT(t) and CC(c) and HH(h))=1.17E+04 ;

*ITU.up(l,g,g1,t)$TT(t)=25;
Q_stage2.up('Pipe',g,g1,t,c,h)$(TT(t) and CC(c) and HH(h) and Npipe(g,g1))=15343;

scalar converged /0/;
parameters log(it,*),log_Stage2(Stage2_it,*),log_Stage2_fesi(Stage2_it,*),lastiteration,
stage2_converged,Stage2_bound,Stage2_fesi_bound,iteration_Stage2,fesi_iteration_Stage2,aa_index,aa_error(it,*),aa_e;
scalar ub/+inf/, lb/-inf/;
scalar  outer_upper /+inf/, outer_lower/-inf/;
Stage2_fixed_dem(g,t,c,h)$TT(t) = Gasdemave(c,h,t,g)$TT(t);

aa_index =1;

loop(it$(not converged),

iter(it) = yes;
iter_fesi(it) = yes;
*iter_fesi(it)$(aa_index = 0) = no;
aa_index = 1;

Stage1_dem(it,g,t,c,h) = Stage2_fixed_dem(g,t,c,h);
display Stage1_dem,iter;

NS.up(s,g,t)$(sc(s) and GS(g,s))=1;
InvP.up('SMRCCS',g,t)$TT(t)=10;
InvP.up('ATRCCS',g,t)$TT(t)=10;
InvP.up('BECCS',g,t)$TT(t)=10;
InvP.up('WE',g,t)$TT(t)=50;
InvS.up('MPSV',g,t)$TT(t)=80;
InvS.up('HPSV',g,t)$TT(t)=80;
InvR.up(e,g,t)$TT(t)=10000;
RI.up(r,t)$TT(t)=rcap(r)/1000;

solve Stage1_Problem_HA2 minimizing Stage1_TC using mip;
display  InvP.l,InvS.l,NP.l,NS.l
AY.l,AYon.l,AYoff.l,Yh.l,Yon.l,Yoff.l, InvR.l,NR.l, CL.l,
RI.l,St.l,IMP.l,Pr.l,Pre.l,Q.l,Qi.l,Qr.l,Qon.l,Qoff.l,
IIC.l,PipeCC.l,PipeOC.l,PCC.l,POC1.l,POC2.l,SCC.l,SOC1.l,SOC2.l,BC.l,GC.l,Stage1_TC.l;

InvP.fx(p,g,t)$TT(t)=round(InvP.l(p,g,t));
InvS.fx(s,g,t)$TT(t)=round(InvS.l(s,g,t));

solve Stage1_Problem minimizing Stage1_TC using mip;
display  InvP.l,InvS.l,NP.l,NS.l
AY.l,AYon.l,AYoff.l,Yh.l,Yon.l,Yoff.l, InvR.l,NR.l, CL.l,
RI.l,St.l,IMP.l,Pr.l,Pre.l,Q.l,Qi.l,Qr.l,Qon.l,Qoff.l,
IIC.l,PipeCC.l,PipeOC.l,PCC.l,POC1.l,POC2.l,SCC.l,SOC1.l,SOC2.l,BC.l,GC.l,Stage1_TC.l;

Stage2_bound = +inf;
stage2_converged = 0;

loop(Stage2_it$(not stage2_converged),
      
    solve Stage2_lower_problem  minimizing Stage2_lower_TC using lp;
    
     if(Stage2_lower_problem.modelstat<3, 
     
      log_Stage2(Stage2_it,'Sta2_error') = abs((Stage2_bound - Stage2_lower_TC.l)/Stage2_bound);
      iteration_Stage2 = ord(Stage2_it);
      
      stage2_converged$(abs((Stage2_bound - Stage2_lower_TC.l)/Stage2_bound)< 1E-8 ) = 1;
      Stage2_bound = Stage2_lower_TC.l;
      
      Lem_stage2(g,t,c,h) = Demand_Stage2_cons.m(g,t,c,h);
      display Stage2_lower_TC.l,Stage2_bound,Stage2_dem.l, Q_Stage2.l, Pr_Stage2.l,Stage2_fixed_dem, Lem_stage2;
      );
      
      if(Stage2_lower_problem.modelstat>2,
       aa_index = 0;
       break;
       );

     solve Stage2_middle_problem maximizing Stage2_middle_TC using lp;
     display Stage2_middle_TC.l, Z_neg.l, Z_pos.l, Stage2_middle_dem.l;
     Stage2_fixed_dem(g,t,c,h)  =  Stage2_middle_dem.l(g,t,c,h);
     
      
     display log_Stage2,iteration_Stage2,Stage2_middle_dem.l;

    );
 display IMP.l,Pr.l,Pre.l,Q.l,Qi.l,Qr.l,Qon.l,Qoff.l;
 
 log_Stage2(Stage2_it,'Stage2_error') = no;
 
outer_lower =  Stage1_TC.l;
outer_upper = Stage2_middle_TC.l + 1000*PCC.l + SCC.l + 1000*PipeCC.l + 1000*PipeOC.l  + ReC.l + POC1.l + SOC1.l;
*               + 1000*sum((p,g,t)$TT(t),dfc(t)*pccost(p,t)*PCap(p)*InvP.l(p,g,t))
*               + sum((s,g,t)$(TT(t) and GS(g,s)),dfc(t)*sccost(s)*SCap(s)*InvS.l(s,g,t)) + 1000*PipeCC.l + 1000*PipeOC.l  + ReC.l
*               + sum((p,g,t)$TT(t),dfo(t)*pocostF(p,t)*PCap(p)*NP.l(p,g,t))
*               + sum((s,g,t)$(TT(t) and GS(g,s)),dfo(t)*socostF(s)*SCap(s)*NS.l(s,g,t));
*1000*PCC + SCC + 1000*PipeCC + POC + SOC + 1000*PipeOC           
             

converged$(abs((outer_upper - outer_lower)/outer_upper) < 1E-3 and ord(it) >1) = 1;
*and ord(it)>3
log(it,'lb') = outer_lower;
log(it,'ub') = outer_upper;
lastiteration = ord(it);
display log;

);

**************************Results
parameters
Pr_final(p,g,t,c,h), St_final(s,g,t,c,h), Q_final(l,g,g1,t,c,h),Qon_final(g,g1,t,c,h),Qoff_final(g,r,t,c,h),
IMP_final(g,t,c,h),SOC2_final,IIC_final,h2MassB_final(g,t,c,h),
POC2_final,CEC_final, BC_final, GC_final;

Pr_final(p,g,t,c,h) = sum(it, Pr.l(it,p,g,t,c,h)$(ord(it) = lastiteration));
St_final(s,g,t,c,h) =  sum(it, St.l(it,s,g,t,c,h)$(ord(it) = lastiteration));
Q_final(l,g,g1,t,c,h) = sum(it, Q.l(it,l,g,g1,t,c,h)$(ord(it) = lastiteration));
Qon_final(g,g1,t,c,h) = sum(it, Qon.l(it,g,g1,t,c,h)$(ord(it) = lastiteration));
Qoff_final(g,r,t,c,h) = sum(it, Qoff.l(it,g,r,t,c,h)$(ord(it) = lastiteration));
IMP_final(g,t,c,h) = sum(it, IMP.l(it,g,t,c,h)$(ord(it) = lastiteration));
POC2_final = sum(it,POC2.l(it)$(ord(it) = lastiteration));
SOC2_final = sum(it,SOC2.l(it)$(ord(it) = lastiteration));
CEC_final =  sum(it,CEC.l(it)$(ord(it) = lastiteration));
BC_final = sum(it, BC.l(it)$(ord(it) = lastiteration));
GC_final = sum(it, GC.l(it)$(ord(it) = lastiteration));
IIC_final = sum(it, IIC.l(it)$(ord(it) = lastiteration));
h2MassB_final(g,t,c,h) = sum(it, h2MassBalance.m(it,g,t,c,h)$(ord(it) = lastiteration));


parameters prodcapacity(p,g,TT),prodcapacityperregion(g,TT),prodcapacitypertech(p,TT), storcapacity(s,g,TT),storcapacityperregion(g,TT);
prodcapacity(p,g,TT)=NP.l(p,g,TT)*Pcap(p) + eps$[ NOT NP.l(p,g,TT)];
prodcapacityperregion(g,TT)=sum(p,prodcapacity(p,g,TT));
prodcapacitypertech(p,TT)= sum(g,prodcapacity(p,g,TT));
storcapacity(s,g,TT)=NS.l(s,g,TT)*Scap(s) + eps$[ NOT NS.l(s,g,TT)];
storcapacityperregion(g,TT)=sum(s,storcapacity(s,g,TT));

parameters prodprofil(p,g,TT,c,h),storprofil(s,g,TT,c,h);
prodprofil(p,g,TT,CC,HH) = Pr_final(p,g,TT,CC,HH)+ eps$[ NOT Pr_final(p,g,TT,CC,HH)];
storprofil(s,g,TT,CC,HH)$GS(g,s) = St_final(s,g,TT,CC,HH)+ eps$[ NOT St_final(s,g,TT,CC,HH)];

parameters H2flow(l,g,g1,TT,c,h),Co2onflow(g,g1,TT,c,h),Co2offflow(g,r,TT,c,h);
H2flow(l,g,g1,TT,CC,HH)=  Q_final(l,g,g1,TT,CC,HH) + eps$[ NOT Q_final(l,g,g1,TT,CC,HH)];
Co2onflow(g,g1,TT,CC,HH)= Qon_final(g,g1,TT,CC,HH) + eps$[ NOT Qon_final(g,g1,TT,CC,HH)];
Co2offflow(g,r,TT,CC,HH)= Qoff_final(g,r,TT,CC,HH) + eps$[ NOT Qoff_final(g,r,TT,CC,HH)];

display prodcapacity,prodcapacityperregion,storcapacity,storcapacityperregion,prodprofil,storprofil;

parameter H2cost(TT),levelised;
H2cost(TT) = sum((g,CC,HH),WF(CC)*h2MassB_final(g,TT,CC,HH))/(card(g)*card(CC)*card(HH));
levelised = Stage1_TC.l/(sum((g,TT,CC,HH),dfo(TT)*WF(CC)*Stage2_fixed_dem(g,TT,CC,HH)));
display H2cost,levelised,h2MassBalance.m;

parameter totalimport(TT),per_import(TT);
totalimport(TT)=sum((g,CC,HH),WF(CC)*IMP_final(g,TT,CC,HH));
per_import(TT)=totalimport(TT)/sum((g,CC,HH),WF(CC)*Stage2_fixed_dem(g,TT,CC,HH)) ;
display totalimport,per_import;

parameter totdem(TT,c,h),totprod(TT,c,h),totstor(TT,c,h);
totdem(TT,CC,HH)=sum((g),Stage2_fixed_dem(g,TT,CC,HH))/1000;
totprod(TT,CC,HH)=sum((p,g),Pr_final(p,g,TT,CC,HH))/1000 ;
totstor(TT,CC,HH)=sum((s,g)$GS(g,s),St_final(s,g,TT,CC,HH))/1000 ;
display totdem,totprod,totstor;

parameter Vbioregional(g,TT);
*Vbioregional(g,TT)= sum((c,h)$(CC(c) and HH(h)),WF(c)*theta*Pr.l('BECCS',g,TT,c,h)/eta('BECCS',TT));
*display Vbioregional;

parameter techload(p,TT), techloadregion(p,g,TT),techloadcluster(p,c,TT);
techload(p,TT)= sum((g,CC,HH),WF(CC)*Pr_final(p,g,TT,CC,HH))/(PCap(p)*sum(g,NP.l(p,g,TT))*8760);
techloadregion(p,g,TT)= sum((CC,HH),WF(CC)*Pr_final(p,g,TT,CC,HH))/(PCap(p)*NP.l(p,g,TT)*8760);
techloadcluster(p,c,TT)$CC(c)= sum((g,HH),Pr_final(p,g,TT,c,HH))/(24*PCap(p)*sum(g,NP.l(p,g,TT)));
parameter techload, techloadregion,techloadcluster;

parameter
param(*)
em(t)
;
em(t)$TT(t) = sum((p,g,c,h)$(CC(c) and HH(h)),WF(c)*y_e(p,t)*theta*Pr_final(p,g,t,c,h));

param("Production Capital Cost") = 1000*PCC.l/1E9;
param("Production Operating Cost") = (POC1.l + POC2_final)/1E9;
param("Storage Capital Cost") = SCC.l/1E9;
param("Storage Operating Cost") = (SOC1.l + SOC2_final)/1E9;
param("Pipeline Cost") = 1000*(PipeCC.l + PipeOC.l)/1E9;
*param("Road Transportation Cost") = RCC.l + ROC.l;
param("Carbon Emissions Cost") = CEC_final/1E9;
param("Imports Cost") = IIC_final/1E9;
param("Renewables Cost")=ReC.l/1E9;
param("Biomass Cost")=BC_final/1E9;
param("Gas Cost")=GC_final/1E9;
param("Total Cost")=Stage1_TC.l/1E9;
display param;

*EXTRACT TO EXCEL
execute_unload  "%outputname%.gdx",
prodcapacity,
prodcapacityperregion,
prodcapacitypertech,
storcapacity,
storcapacityperregion,
prodprofil,
storprofil,
H2flow,
Co2onflow,
Co2offflow,
AY.l,
AYon.l,
AYoff.l,
Q.l,Qon.l,Qoff.l
em,
emtarget,
param,
H2cost,
levelised,
Stage2_fixed_dem,
totalimport,
per_import
totdem,totprod,totstor,
Vbioregional,BA,
techload, techloadregion,techloadcluster;

$onecho >  exceloutput.txt

epsout=0
par=prodcapacity rng=prodcapacity!                       Rdim=2 Cdim=1
par=prodcapacityperregion rng=prodcapacityperregion!     Rdim=1 Cdim=1
par=prodcapacitypertech rng=prodcapacitypertech!         Rdim=1 Cdim=1
par=storcapacity rng=storcapacity!                       Rdim=2 Cdim=1
par=storcapacityperregion rng=storcapacityperregion!     Rdim=1 Cdim=1
par=prodprofil rng=prodprofil!                           Rdim=3 Cdim=2
par=storprofil rng=storprofil!                           Rdim=3 Cdim=2
var=AY.L rng=Pipelines!A1                                Rdim=3 Cdim=1
var=AYon.L rng=Pipelines!I1                              Rdim=3 Cdim=1
var=AYoff.L rng=Pipelines!Q1                             Rdim=3 Cdim=1
par=H2flow rng=H2flows!                                  Rdim=5 cdim=1
par=Co2onflow rng=CO2onflows!                            Rdim=4 cdim=1
par=Co2offflow rng=CO2offflows!                          Rdim=4 cdim=1
par=param rng=Costs!                                     Rdim=1
par=levelised rng=H2levelised!
par=H2cost  rng=H2levelised!A4                           Rdim=1
par=em    rng=Emissions!A1                               Cdim=1
par=emtarget    rng=Emissions!F1                         Cdim=1
par=Stage2_fixed_dem     rng=Demand!A1                   Rdim=3 cdim=1
par=totalimport rng=imports!A1                           Cdim=1
par=per_import rng=imports!A5                            Cdim=1
par=totdem       rng=totprofiles!A1                      rdim=2 cdim=1
par=totprod      rng=totprofiles!A40                     rdim=2 cdim=1
par=Vbioregional rng=Vbioregional!a1                     rdim=1 cdim=1
par=totstor      rng=totprofiles!a80                     rdim=2 cdim=1
par=Vbioregional rng=Vbioregional!a1                     rdim=1 cdim=1
par=BA           rng=Vbioregional!a18                    rdim=1 cdim=1
par=techload     rng=techload!A1                         rdim=1 cdim=1
par=techloadregion   rng=techload!A10                    rdim=1 cdim=2
par=techloadcluster      rng=techload!A20                rdim=1 cdim=2
$offecho

Execute 'gdxxrw %outputname%.gdx output=%outputname%.xlsx @exceloutput';