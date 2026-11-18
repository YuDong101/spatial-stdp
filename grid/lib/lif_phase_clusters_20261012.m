function out = lif_phase_clusters_20261012(q, supplied)
% Persistent phase clusters, then full-cycle slips. q is node x time in cycles.
% Cluster detection never uses weights, D, or the historical <=.1 lock mask.
% Uniform node rotations preserve each temporal trace. Their maximal circular
% density over the entire sampled window gives a conservative null cutoff.
cfg=struct('bandwidth',.03,'bins',256,'sampleStride',50,'nullDraws',199,...
    'nullQuantile',.99,'nullSeed',61012,'persistence',.8,'minimumMembers',5,...
    'maximumRadius',.15,'minimumValleyDrop',.20,'reference',[],...
    'methodVersion','persistent_density_coassignment_20261013_v3_safe_merge');
if nargin>1,keys=fieldnames(supplied);for j=1:numel(keys),cfg.(keys{j})=supplied.(keys{j});end,end
cfg.methodVersion='persistent_density_coassignment_20261013_v3_safe_merge';
q=double(q);[n,T]=size(q);assert(n>=2 && T>=2 && all(isfinite(q),'all'));
assert(cfg.bins>=64 && cfg.persistence>.5 && cfg.persistence<=1);
ix=unique([1:cfg.sampleStride:T T]);phase=mod(q(:,ix)-q(1,ix),1);
kernel=exp(-.5*(min((0:cfg.bins-1)',(cfg.bins:-1:1)')/cfg.bins/cfg.bandwidth).^2);
kernel=kernel/sum(kernel)*cfg.bins;kernelFFT=fft(kernel);
density=circular_density(phase,kernelFFT,cfg.bins);
rs=RandStream('mt19937ar','Seed',cfg.nullSeed);nullMaximum=zeros(cfg.nullDraws,1);
for b=1:cfg.nullDraws
    rotated=mod(phase+rand(rs,n,1),1);
    surrogateDensity=circular_density(rotated,kernelFFT,cfg.bins);
    nullMaximum(b)=max(surrogateDensity,[],'all');
end
sortedNull=sort(nullMaximum);rank=min(cfg.nullDraws,ceil((cfg.nullDraws+1)*cfg.nullQuantile));
cutoff=sortedNull(rank);
instant=zeros(n,numel(ix));instantK=zeros(1,numel(ix));
for t=1:numel(ix)
    f=density(:,t);peaks=find(f>circshift(f,1)&f>=circshift(f,-1)&f>cutoff);
    % Merge shallow adjacent modes; retain the higher peak deterministically.
    while numel(peaks)>1
        heights=f(peaks);drop=ones(size(peaks));
        for j=1:numel(peaks)
            k=mod(j,numel(peaks))+1;a=peaks(j);b=peaks(k);
            if b<=a,b=b+cfg.bins;end
            valley=min(f(mod((a:b)-1,cfg.bins)+1));
            drop(j)=1-valley/min(heights(j),heights(k));
        end
        [shallow,j]=min(drop);if shallow>=cfg.minimumValleyDrop,break,end
        k=mod(j,numel(peaks))+1;
        if f(peaks(j))>=f(peaks(k)),peaks(k)=[];else,peaks(j)=[];end
    end
    if isempty(peaks),continue,end
    centers=(peaks-1)/cfg.bins;
    distance=abs(mod(phase(:,t)-centers'+.5,1)-.5);
    [nearest,label]=min(distance,[],2);label(nearest>cfg.maximumRadius)=0;
    for j=1:numel(peaks),if sum(label==j)<cfg.minimumMembers,label(label==j)=0;end,end
    live=unique(label(label>0));instantK(t)=numel(live);
    for j=1:numel(live),instant(label==live(j),t)=j;end
end
% Coassignment is invariant to per-frame mode numbering and global rotation.
P=zeros(n);for t=1:numel(ix),a=instant(:,t);P=P+double(a==a' & a>0 & a'>0);end
P=P/numel(ix);presence=diag(P);
% Average-linkage consensus permits breathing clouds while avoiding a
% nearest-center partition of nodes with no supported instantaneous cluster.
distance=max(0,min(1,1-P));distance(1:n+1:end)=0;
tree=linkage(squareform(distance),'average');group=cluster(tree,'cutoff',1-cfg.persistence,'criterion','distance');
clusterId=zeros(n,1);K=0;
for g=unique(group)'
    members=find(group==g & presence>=cfg.persistence);
    if numel(members)<cfg.minimumMembers,continue,end
    % Mean agreement with the other members must meet the same persistence.
    while numel(members)>=cfg.minimumMembers
        agreement=(sum(P(members,members),2)-diag(P(members,members)))/(numel(members)-1);
        bad=agreement<cfg.persistence-1e-12;if ~any(bad),break,end
        members(bad)=[];
    end
    if numel(members)>=cfg.minimumMembers,K=K+1;clusterId(members)=K;end
end
% Overlapping groups may be joined only if the union preserves reliability
% for every member. A slipping group can otherwise bridge unrelated stable
% groups and cause all of them to be removed after a transitive union.
[clusterId,unresolvedSeparation]=lif_phase_safe_merge(P,clusterId,cfg.persistence);
K=max(clusterId);
% Stable numeric labels: descending population, then smallest original node.
if K>0
    sizes=arrayfun(@(k)sum(clusterId==k),(1:K)');first=arrayfun(@(k)find(clusterId==k,1),(1:K)');
    [~,order]=sortrows([-sizes first],[1 2]);old=clusterId;
    for k=1:K,clusterId(old==order(k))=k;end
    unresolvedSeparation=unresolvedSeparation(order,order);
end
membershipConfidence=zeros(n,1);
for k=1:K
    m=find(clusterId==k);membershipConfidence(m)=(sum(P(m,m),2)-diag(P(m,m)))/(numel(m)-1);
end
% Only after K and membership exist do we obtain the inherited D reference.
if isempty(cfg.reference)
    span=zeros(n,n);for i=1:n,d=q-q(i,:);span(i,:)=max(d,[],2)-min(d,[],2);end
    [~,reference]=max(sum(span<=.1,2));referencePolicy='legacy_max_excursion_0p1_cohort_reference_only';
else,reference=cfg.reference;referencePolicy='explicit_common_parent_reference';end
rel=q-q(reference,:);excursion=max(rel,[],2)-min(rel,[],2);slip=excursion>=1;
displayId=clusterId;displayId(slip)=-1;
out=struct('K',K,'cluster_id',clusterId,'D',slip,'display_id',displayId,...
    'unclustered',clusterId==0 & ~slip,'reference',reference,'reference_policy',referencePolicy,...
    'excursion_cycles',excursion,'membership_confidence',membershipConfidence,...
    'cluster_presence',presence,'instantaneous_K',instantK,'instantaneous_labels',instant,...
    'sample_indices',ix,'density_null_cutoff',cutoff,'null_max_density',nullMaximum,...
    'coassignment',P,'subgroup_separation_unresolved',unresolvedSeparation,'config',cfg);
end

function f=circular_density(p,kernelFFT,bins)
[n,T]=size(p);bin=mod(round(p*bins),bins)+1;
indices=bin+(0:T-1)*bins;counts=reshape(accumarray(indices(:),1,[bins*T 1]),bins,T)/n;
f=real(ifft(fft(counts,[],1).*kernelFFT,[],1));
end
