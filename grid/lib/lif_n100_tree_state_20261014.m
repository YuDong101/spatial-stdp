function r=lif_n100_tree_state_20261014(q,t,old,cfg)
% Five physical leaves; detection sensitivity is a separate diagnostic.
% Codes: 0 unresolved,1 dispersed,2 phase clusters,3 partial,4 lock+slip,5 dispersed slip.
if nargin<4,cfg=struct;end
defaults=struct('coreSpan',.35,'coreWidth',.35,'persistence',.8,'minCore',5,...
    'halfWidth',.175,'escapeNet',.20,'escapeTail',.03,'captureTail',.02,...
    'dispersionWidth',.5,'sampleStride',50);
for f=string(fieldnames(defaults))',if ~isfield(cfg,f),cfg.(f)=defaults.(f);end,end
q=double(q);t=double(t(:)')/1000;[n,T]=size(q);ix=unique([1:cfg.sampleStride:T T]);
rel=q-q(old.reference,:);span=max(rel,[],2)-min(rel,[],2);D=span>1;
r=struct('code',0,'reason',"unresolved_locking_evidence",'reference',old.reference,...
 'D',D,'phase_range',span,'exact_one',sum(span==1),'near_one',sum(abs(span-1)<=1e-10),...
 'old_D_changes',sum(D~=old.D),'support_mask',false(n,1),'core_count_by_band',zeros(5,1),...
 'core_support_by_band',zeros(5,1),'candidate_count',0,'verified_core_count',0,...
 'core_diagnostics',struct([]),'macroK',nan,'macro_id',zeros(n,1),...
 'escape_mask',false(n,1),'capture_mask',false(n,1),'unresolved_mask',false(n,1),...
 'partition_codes',zeros(0,1),'partition_diagnostics',{{}},'density_sensitive',false,...
 'q',q(:,ix),'time',t(ix)*1000,'inputHash',old.inputHash,'config',cfg);
% Pairwise full-window excursions are independent of the density labels.
pairspan=nan(n);live=find(~D);
if any(cellfun(@(g)g.K>0,old.groups))
 for a=live',z=q(live,:)-q(a,:);pairspan(a,live)=max(z,[],2)-min(z,[],2);end
end
cores=cell(5,1);supportByBand=false(n,5);diag=struct([]);
for b=1:5
 ids=old.groups{b}.cluster_id;cohorts={};
 for k=1:max(ids)
  original=find(ids==k & ~D);
  if numel(original)<cfg.minCore,continue,end
  r.candidate_count=r.candidate_count+1;m=original;
  % Remove the member with most pairwise violations; deterministic node tie.
  while numel(m)>=cfg.minCore
   bad=pairspan(m,m)>cfg.coreSpan+1e-12;
   if ~any(bad,'all'),break,end
   degree=sum(bad,2);[~,remove]=max(degree);m(remove)=[];
  end
  accepted=false;w90=nan;w100=nan;frac=[0 0 0];ps=nan;
  if numel(m)>=cfg.minCore
   w90=arcwidth(q(m,:),.9);w100=arcwidth(q(m,:),1);
   pass=w90<=cfg.coreWidth+1e-12;mid=floor(T/2);
   frac=[mean(pass),mean(pass(1:mid)),mean(pass(mid+1:end))];
   ps=max(pairspan(m,m),[],'all');accepted=all(frac>=cfg.persistence-1e-12);
  end
  d=struct('band',b,'subgroup',k,'original_n',numel(original),'retained_n',numel(m),...
   'members',m,'accepted',accepted,'pair_span_max',ps,'width90_max',max(w90),...
   'width100_max',max(w100),'width_pass_whole_early_late',frac);
  diag=[diag;d]; %#ok<AGROW>
  if accepted,cohorts{end+1}=m;supportByBand(m,b)=true;end %#ok<AGROW>
 end
 cores{b}=cohorts;r.core_count_by_band(b)=numel(cohorts);r.core_support_by_band(b)=sum(supportByBand(:,b));
end
r.core_diagnostics=diag;r.support_mask=any(supportByBand,2);
r.verified_core_count=sum(r.core_count_by_band);r.density_sensitive=any(r.core_count_by_band==0)&&any(r.core_count_by_band>0);
% Positive dispersion evidence, independently of mean frequency or slip fraction.
w=arcwidth(q(:,ix),.9);half=floor(numel(ix)/2);
r.dispersion_width90_median=median(w);
r.dispersion_fraction=[mean(w>cfg.dispersionWidth),mean(w(1:half)>cfg.dispersionWidth),mean(w(half+1:end)>cfg.dispersionWidth)];
if ~any(r.support_mask)
 noCandidates=all(cellfun(@(g)g.K==0,old.groups));
 if noCandidates && all(r.dispersion_fraction>=cfg.persistence)
  if any(D),r.code=5;r.reason="dispersed_with_observed_slip";
  else,r.code=1;r.reason="dispersed_without_observed_slip";end
 elseif noCandidates,r.reason="no_core_but_dispersion_unconfirmed";
 else,r.reason="candidate_core_not_verified";end
 return
end
if any(D),r.code=4;r.reason="verified_core_plus_observed_slip";return,end
% Different cluster counts do not veto a common physical leaf. Only a direct
% phase-cluster versus persistent-escape conflict leaves the state unresolved.
parts={};seen={};
for b=[4 3 5 2 1]
 if isempty(cores{b}),continue,end
 % The macro geometry uses only independently verified fixed cores. Pruned
 % moving members cannot shift its centers or veto the core partition.
 ids=zeros(n,1);for k=1:numel(cores{b}),ids(cores{b}{k})=k;end
 a=lif_phase_macro_candidate(q(:,ix),ids,.5,.35,false(max(ids)));
 if ~isfinite(a.candidate_K),continue,end
 macro=a.member_id;M=a.candidate_K;coreid=zeros(n,1);
 for k=1:M,coreid(macro==k & supportByBand(:,b))=k;end
 if any(arrayfun(@(k)sum(coreid==k)<cfg.minCore,1:M)),continue,end
 if any(cellfun(@(v)isequal(v,coreid),seen)),continue,end
 seen{end+1}=coreid; %#ok<AGROW>
 p=member_evidence(q,t,coreid,cfg);p.band=b;p.macroK=M;p.core_id=coreid;
 parts{end+1}=p; %#ok<AGROW>
end
r.partition_diagnostics=parts;
if isempty(parts),r.reason="verified_core_but_member_partition_unresolved";return,end
r.partition_codes=cellfun(@(p)p.code,parts)';positive=unique(r.partition_codes(r.partition_codes>0));
if numel(positive)>1,r.reason="conflicting_direct_member_evidence";return,end
if isempty(positive),[~,pick]=min(cellfun(@(p)sum(p.unresolved),parts));r.reason="member_dynamics_unresolved";
else,r.code=positive;pick=find(r.partition_codes==r.code,1);
 if r.code==3,r.reason="verified_core_plus_persistent_escape";else,r.reason="verified_phase_cluster_membership";end
end
p=parts{pick};r.macroK=p.macroK;r.macro_id=p.membership;r.escape_mask=p.escape;
r.capture_mask=p.capture;r.unresolved_mask=p.unresolved;r.member_diagnostics=p;
end

function p=member_evidence(q,t,coreid,cfg)
n=size(q,1);M=max(coreid);early=t<=t(1)+.5;late=t>=t(end)-.5;
ta=t>=19&t<=19.2;tb=t>=19.8;last=t>=19;v=exp(2i*pi*q);
centers=zeros(M,size(q,2));resultant=zeros(M,size(q,2));cost=zeros(n,M);
for k=1:M
 z=mean(v(coreid==k,:),1);centers(k,:)=unwrap(angle(z))/(2*pi);resultant(k,:)=abs(z);
 circular=abs(mod(q(:,early)-centers(k,early)+.5,1)-.5);cost(:,k)=median(circular,2);
end
[~,membership]=min(cost,[],2);membership(coreid>0)=coreid(coreid>0);
res=zeros(size(q));
for node=1:n
 m=find(coreid==membership(node));m(m==node)=[];
 c=unwrap(angle(mean(v(m,:),1)))/(2*pi);z=q(node,:)-c;
 res(node,:)=z-round(median(z(early)));
end
net=median(res(:,late),2)-median(res(:,early),2);tail=median(res(:,tb),2)-median(res(:,ta),2);
residence=mean(abs(res)<=cfg.halfWidth,2);endResidence=mean(abs(res(:,last))<=cfg.halfWidth,2);
% Count a sustained outside episode, not just repeated isolated crossings.
outside=abs(res)>cfg.halfWidth;longest=zeros(n,1);run=zeros(n,1);
for j=1:size(q,2),run=(run+1).*outside(:,j);longest=max(longest,run);end
duration=longest*median(diff(t));
escape=residence<cfg.persistence & abs(net)>=cfg.escapeNet & abs(tail)>=cfg.escapeTail & net.*tail>0 & duration>=.5;
gridCounts=zeros(27,1);j=0;
for width=[.15 .175 .2]
 run=zeros(n,1);dur=zeros(n,1);outsideGrid=abs(res)>width;
 for k=1:size(q,2),run=(run+1).*outsideGrid(:,k);dur=max(dur,run);end
 dur=dur*median(diff(t));
 for shift=[.15 .2 .25],for speed=[.02 .03 .05]
  j=j+1;gridCounts(j)=sum(mean(abs(res)<=width,2)<.8 & abs(net)>=shift & abs(tail)>=speed & net.*tail>0 & dur>=.5);
 end,end
end
% Captured members must end inside and approach the fixed core. They are
% reported separately; slow unexplained exterior motion is not waived.
capture=residence<cfg.persistence & endResidence>=cfg.persistence & abs(tail)<cfg.captureTail & ...
 median(abs(res(:,late)),2)<median(abs(res(:,early)),2);
inside=residence>=cfg.persistence & endResidence>=cfg.persistence;unresolved=~(inside|capture|escape);
initialMembership=membership;captureDestination=zeros(n,1);
% An initially exterior node may approach a different fixed destination.
% Verify that destination over the whole trace; never reassign per frame.
for node=find(unresolved & coreid==0)'
 best=inf;destination=0;
 for k=1:M
  z=q(node,:)-centers(k,:);z=z-round(median(z(late)));
  terminal=mean(abs(z(last))<=cfg.halfWidth);speed=median(z(tb))-median(z(ta));
  approach=median(abs(z(late)))<median(abs(z(early)));
  if terminal>=cfg.persistence && abs(speed)<cfg.captureTail && approach
   cost=median(abs(z(late)));if cost<best,best=cost;destination=k;end
  end
 end
 if destination>0,capture(node)=true;unresolved(node)=false;captureDestination(node)=destination;membership(node)=destination;end
end
code=0;if any(escape),code=3;elseif ~any(unresolved),code=2;end
p=struct('code',code,'membership',membership,'escape',escape,'capture',capture,'unresolved',unresolved,...
 'initial_membership',initialMembership,'capture_destination',captureDestination,...
 'residence',residence,'end_residence',endResidence,'net',net,'tail',tail,'outside_duration',duration,...
 'escape_counts27',gridCounts,'core_resultant_min',min(resultant,[],'all'),...
 'centers_sample',centers(:,1:50:end),'residual_sample',res(:,1:50:end));
end

function w=arcwidth(q,fraction)
p=sort(mod(q,1),1);n=size(p,1);count=ceil(fraction*n);pp=[p;p+1];
w=min(pp(count:count+n-1,:)-p,[],1);
end
