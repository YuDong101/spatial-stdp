function x = lif_make_inputs(p,coords)
%LIF_MAKE_INPUTS Separate random streams; W(target,source).
n=p.N;off=~eye(n);assert(size(coords,1)==n && size(coords,2)==3);
d=sqrt(sum((reshape(coords,n,1,3)-reshape(coords,1,n,3)).^2,3));
dref=median(d(off & d>0));assert(isfinite(dref) && dref>0);
stream=RandStream('mt19937ar','Seed',p.weightSeed);
w=10.^(-6+(log10(p.gmax)+6)*rand(stream,n));
if p.initialWeightFloor~=1e-6
 w=10.^(log10(p.initialWeightFloor)+(log10(p.gmax)-log10(p.initialWeightFloor))*rand(RandStream('mt19937ar','Seed',p.weightSeed),n));
end
w(~off)=0;
phase=rand(RandStream('mt19937ar','Seed',p.phaseSeed),1,n);
vInf=p.EL+p.Iext/p.gL;assert(vInf>p.Vthreshold,'Initialization requires tonic isolated firing.');
tau=p.C/p.gL;period=p.refractoryMs+tau*log((vInf-p.Vreset)/(vInf-p.Vthreshold));
age=phase*period;v=vInf+(p.Vreset-vInf)*exp(-max(0,age-p.refractoryMs)/tau);
v(age<p.refractoryMs)=p.Vreset;until=max(0,p.refractoryMs-age);
if strcmp(p.transport,'spatial')
 delayMs=p.tau0Ms+dref*log1p(d/dref)/p.conductionSpeed;
elseif strcmp(p.transport,'source_gaussian')
 a=p.delayMeanMs+p.delayStdMs*randn(RandStream('mt19937ar','Seed',p.delaySeed),1,n);
 a=min(p.delayLimitMs,max(0,a));delayMs=repmat(a,n,1);
else,error('Unknown transport');end
delaySteps=max(1,round(delayMs/p.dtMs));delaySteps(~off)=1;
k=(1+d/dref).^(-p.alpha);k(~off)=0;
switch p.arm
 case 'S'
 case 'U',k(off)=mean(k(off));
 case 'P'
  u=find(triu(off,1));vals=k(u);s=RandStream('mt19937ar','Seed',p.shuffleSeed);
  vals=vals(randperm(s,numel(vals)));k=zeros(n);k(u)=vals;k=k+k';
 case 'N',k(:)=0;
 otherwise,error('Unknown arm');
end
bins=ones(n);bins(d/dref>0.5)=2;bins(d/dref>1)=3;bins(d/dref>2)=4;
x=struct('W',w,'V',v,'refractoryUntil',until,'delaySteps',delaySteps,'delayRequestedMs',delayMs,'kernel',k,'distance',d,'distanceReference',dref,'distanceBins',bins,'coordinates',coords,'initialPhaseFraction',phase,'isolatedPeriodMs',period);
end
