function out = lif_phase_macro_candidate(q, ids, gapRatio, width, unresolved)
% EXPLORATORY ONLY: circular-gap scale diagnostic, not a production classifier.
% q contains sampled node phases in cycles. ids are persistent MICRO groups.
% The gap scale has not been approved as a scientific classification rule.
% Every original member is retained. Neither weights nor physical D are used.
q=double(q);ids=double(ids(:));K=max(ids);T=size(q,2);
assert(size(q,1)==numel(ids) && all(isfinite(q),'all'));
assert(gapRatio>0 && gapRatio<=1 && width>0 && width<.5);
assert(all(ids>=0 & ids==fix(ids)) && isequal(unique(ids(ids>0))',(1:K)));
if nargin<5,unresolved=false(K);end
assert(isequal(size(unresolved),[K K]));
out=struct('gap_ratio',gapRatio,'width_limit',width,'subgroup_K',K,...
    'proposed_K',0,'candidate_K',nan,'partition_fraction',nan,...
    'joint_fraction',nan,'geometry_pass',false,'conflict',any(unresolved,'all'),...
    'micro_conflict',any(unresolved,'all'),'cross_macro_conflict',false,...
    'within_unresolved_pairs',0,'cross_unresolved_pairs',0,...
    'modal_tie',false,'minimum_center_resultant',nan,...
    'subgroup_partition',zeros(K,1),'member_id',zeros(size(ids)),...
    'status',"no_supported_subgroups",'centers',zeros(K,T),...
    'instant_partition',zeros(K,T),'normalized_gaps',zeros(K,T));
if K==0,return,end
v=exp(2i*pi*q);resultants=zeros(K,T);
for k=1:K
    z=mean(v(ids==k,:),1);out.centers(k,:)=mod(angle(z)/(2*pi),1);resultants(k,:)=abs(z);
end
out.minimum_center_resultant=min(resultants,[],'all');
for t=1:T
    [c,order]=sort(out.centers(:,t));g=diff([c;c(1)+1]);
    out.normalized_gaps(:,t)=sort(g/max(g));
    cut=g>=gapRatio*max(g)-1e-12;
    % Walk from an actual cut, then canonicalize by subgroup identity.
    first=mod(find(cut,1),K)+1;labels=zeros(K,1);label=1;
    for j=0:K-1
        p=mod(first+j-1,K)+1;labels(order(p))=label;
        if cut(p),label=label+1;end
    end
    out.instant_partition(:,t)=canonical(labels);
end
[patterns,~,index]=unique(out.instant_partition','rows');
counts=accumarray(index,1);[count,win]=max(counts);out.modal_tie=sum(counts==count)>1;
partition=patterns(win,:)';M=max(partition);out.proposed_K=M;
out.subgroup_partition=partition;out.partition_fraction=count/T;
out.cross_macro_conflict=any(unresolved & partition~=partition','all');
out.cross_unresolved_pairs=nnz(triu(unresolved & partition~=partition',1));
out.within_unresolved_pairs=nnz(triu(unresolved & partition==partition',1));
out.conflict=out.cross_macro_conflict;
out.member_id(ids>0)=partition(ids(ids>0));
out.center_width=zeros(M,T);out.member90_width=zeros(M,T);out.member100_width=zeros(M,T);
starts=zeros(M,T);out.arc_start_tied=false(M,T);
for m=1:M
    out.center_width(m,:)=arc(out.centers(partition==m,:),1);
    [out.member90_width(m,:),starts(m,:),out.arc_start_tied(m,:)]=arc(q(out.member_id==m,:),.9);
    out.member100_width(m,:)=arc(q(out.member_id==m,:),1);
end
% Extra diagnostic only, not a user-approved classification gate. If an arc
% has several equally short starts, overlap depends on the chosen start.
out.bulk_arcs_disjoint=true(1,T);
for a=1:M,for b=a+1:M
    for shift=-1:1
        lo=max(starts(a,:),starts(b,:)+shift);
        hi=min(starts(a,:)+out.member90_width(a,:),starts(b,:)+shift+out.member90_width(b,:));
        out.bulk_arcs_disjoint=out.bulk_arcs_disjoint & hi-lo<=1e-12;
    end
end,end
same=all(out.instant_partition==partition,1);
out.joint_samples=same & all(out.center_width<=width+1e-12,1) & ...
    all(out.member90_width<=width+1e-12,1);
out.joint_fraction=mean(out.joint_samples);out.geometry_pass=out.joint_fraction>=.8-1e-12 && ...
    ~out.modal_tie && out.minimum_center_resultant>1e-12;
out.status="unresolved_geometry_or_partition";
if out.geometry_pass
    out.status="candidate_geometry_pass";
    if out.conflict,out.status="unresolved_coassignment";else,out.candidate_K=M;end
end
end

function labels=canonical(labels)
old=labels;live=unique(old,'stable');
for k=1:numel(live),labels(old==live(k))=k;end
end

function [w,start,tied]=arc(p,fraction)
p=sort(mod(p,1),1);n=size(p,1);count=ceil(fraction*n);pp=[p;p+1];
widths=pp(count:count+n-1,:)-p;[w,index]=min(widths,[],1);
start=p(sub2ind(size(p),index,1:size(p,2)));
tied=false(1,size(p,2));
for t=1:size(p,2)
    starts=p(abs(widths(:,t)-w(t))<=1e-12,t);
    tied(t)=numel(starts)>1 && max(starts)-min(starts)>1e-12;
end
end
