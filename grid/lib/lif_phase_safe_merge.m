function [labels, unresolved] = lif_phase_safe_merge(P, labels, persistence)
% Merge reliable subgroups only when every existing member remains reliable.
% Incompatible overlap is retained as unresolved, never a transitive union
% followed by deletion of previously supported members.
labels=labels(:);threshold=1-persistence;
while true
    live=unique(labels(labels>0));candidates=zeros(0,3);
    for a=1:numel(live),for b=a+1:numel(live)
        shared=mean(P(labels==live(a),labels==live(b)),'all');
        if shared>threshold,candidates(end+1,:)=[shared live(a) live(b)];end %#ok<AGROW>
    end,end
    if isempty(candidates),break,end
    candidates=sortrows(candidates,[-1 2 3]);accepted=false;
    for pair=1:size(candidates,1)
        a=candidates(pair,2);b=candidates(pair,3);members=find(labels==a|labels==b);
        agreement=(sum(P(members,members),2)-diag(P(members,members)))/(numel(members)-1);
        if all(agreement>=persistence-1e-12)
            labels(labels==b)=a;accepted=true;break
        end
    end
    if ~accepted,break,end
end
live=unique(labels(labels>0));old=labels;
for j=1:numel(live),labels(old==live(j))=j;end
unresolved=false(numel(live));
for a=1:numel(live),for b=a+1:numel(live)
    shared=mean(P(labels==a,labels==b),'all');
    unresolved(a,b)=shared>threshold;unresolved(b,a)=unresolved(a,b);
end,end
end
