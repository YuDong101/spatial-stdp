function result = lif_simulate(p,x,startState)
%LIF_SIMULATE Edge-delay conductance LIF with restartable full state.
% Event order: arrivals/LTD at t; membrane step; post spikes/LTP at t+dt.
n=p.N;dt=p.dtMs;assert(dt>0 && dt<p.synTauMs && p.synJump>0);
assert(all(x.delaySteps(:)>=1) && all(x.delaySteps(:)==round(x.delaySteps(:))));
assert(all(x.W(:)>=0 & x.W(:)<=p.gmax) && all(diag(x.W)==0));
assert(ismember(p.stdpClock,{'somatic','arrival'}));
lastStep=round(p.endMs/dt);assert(abs(lastStep*dt-p.endMs)<1e-8);
weightStride=round(p.weightEveryMs/dt);traceStride=round(p.traceEveryMs/dt);
assert(weightStride>=1 && traceStride>=1);
if nargin<3 || isempty(startState)
 st=struct('step',0,'dtMs',dt,'V',x.V,'W',x.W,'s',zeros(n), ...
 'refractoryUntil',x.refractoryUntil,'lastSpike',-inf(1,n),'lastArrival',-inf(n), ...
 'buffer',false(max(x.delaySteps,[],'all')+1,n),'updateCounts',zeros(4,6));
else
 st=startState;assert(st.dtMs==dt && st.step<lastStep && size(st.buffer,2)==n);
 assert(size(st.buffer,1)>max(x.delaySteps,[],'all'));
end
firstStep=st.step;nSteps=lastStep-firstStep;
maxSpikes=ceil((p.endMs-firstStep*dt)/p.refractoryMs+2)*n;
times=zeros(1,maxSpikes);ids=zeros(1,maxSpikes);spikeCount=0;
frameSteps=unique([firstStep,ceil((firstStep+1)/weightStride)*weightStride:weightStride:lastStep,lastStep]);
weights=zeros(numel(frameSteps),n,n);updates=zeros(numel(frameSteps),4,6);
weights(1,:,:)=st.W;updates(1,:,:)=st.updateCounts;nextFrame=2;
traceSteps=firstStep+1: min(lastStep,floor(p.traceUntilMs/dt));traceSteps=traceSteps(mod(traceSteps,traceStride)==0);
voltage=zeros(numel(traceSteps),n);current=zeros(numel(traceSteps),n);gate=zeros(numel(traceSteps),n);traceIndex=1;
isLearning=~strcmp(p.arm,'N');off=~eye(n);clock=tic;
for stepIndex=firstStep:lastStep-1
 t=stepIndex*dt;tNext=(stepIndex+1)*dt;
 arrival=lif_read_arrivals(st.buffer,stepIndex,x.delaySteps);
 st.s=st.s*(1-dt/p.synTauMs)+p.synJump*arrival;
 if strcmp(p.stdpClock,'arrival')
  st.lastArrival(arrival)=t;
  if isLearning && any(arrival,'all')
   lag=st.lastSpike(:)-t;mask=arrival & lag<0 & lag>=-p.stdpCutoffMs;
   dw=-p.learningScale*p.Aminus*exp(lag/p.tauMinusMs).*x.kernel.*mask;
   [st.W,c]=lif_apply_update(st.W,dw,p.gmax,x.distanceBins);st.updateCounts=st.updateCounts+c;
  end
 end
 G=sum(st.W.*st.s,2)';I=G.*(p.Esyn-st.V);
 vNew=st.V+dt*(-p.gL*(st.V-p.EL)+p.Iext+I)/p.C;
 refractory=t<st.refractoryUntil-1e-10;vNew(refractory)=p.Vreset;
 spiking=vNew>p.Vthreshold & ~refractory;
 vNew(spiking)=p.Vreset;st.refractoryUntil(spiking)=tNext+p.refractoryMs;
 st.lastSpike(spiking)=tNext;st.V=vNew;
 st.buffer(mod(stepIndex+1,size(st.buffer,1))+1,:)=spiking;
 if any(spiking)
  nodeIds=find(spiking);m=numel(nodeIds);ii=spikeCount+(1:m);times(ii)=tNext;ids(ii)=nodeIds;spikeCount=spikeCount+m;
  if isLearning
   if strcmp(p.stdpClock,'somatic')
    dw=lif_somatic_update(p,st.lastSpike,spiking,x.kernel);
   else
    lag=tNext-st.lastArrival;mask=spiking(:) & off & lag>0 & lag<=p.stdpCutoffMs;
    dw=p.learningScale*p.Aplus*exp(-lag/p.tauPlusMs).*x.kernel.*mask;
   end
   [st.W,c]=lif_apply_update(st.W,dw,p.gmax,x.distanceBins);st.updateCounts=st.updateCounts+c;
  end
 end
 st.step=stepIndex+1;
 if nextFrame<=numel(frameSteps) && st.step==frameSteps(nextFrame)
  weights(nextFrame,:,:)=st.W;updates(nextFrame,:,:)=st.updateCounts;nextFrame=nextFrame+1;
 end
 if traceIndex<=numel(traceSteps) && st.step==traceSteps(traceIndex)
  voltage(traceIndex,:)=st.V;current(traceIndex,:)=I;gate(traceIndex,:)=G;traceIndex=traceIndex+1;
 end
 assert(all(isfinite(st.V)) && all(st.s(:)>=0),'Invalid LIF state');
end
elapsed=toc(clock);times=times(1:spikeCount);ids=ids(1:spikeCount);
windowStart=max(p.burnInMs,firstStep*dt);duration=p.endMs-windowStart;
rate=nan(n,1);if duration>0,rate=accumarray(ids(times>windowStart)',1,[n,1])/(duration/1000);end
result=struct('parameters',p,'inputs',x,'timesMs',times,'nodeIds',ids,'rateHz',rate, ...
 'frameTimesMs',frameSteps*dt,'weights',weights,'cumulativeUpdates',updates, ...
 'traceTimesMs',traceSteps*dt,'currentTimesMs',traceSteps*dt-dt, ...
 'voltage',voltage,'synapticCurrent',current,'inputConductance',gate, ...
 'finalState',st,'runtimeSeconds',elapsed,'steps',nSteps);
end
