function [w,counts] = lif_apply_update(w,dw,gmax,bins)
%LIF_APPLY_UPDATE Requested/applied LTP,LTD and discarded clipping by distance.
counts=zeros(4,6);idx=find(dw~=0);
if isempty(idx),return;end
requested=dw(idx);old=w(idx);new=min(gmax,max(0,old+requested));applied=new-old;w(idx)=new;
a=[max(requested,0),max(-requested,0),max(applied,0),max(-applied,0),max(old+requested-gmax,0),max(-(old+requested),0)];
for j=1:6,counts(:,j)=accumarray(bins(idx),a(:,j),[4,1]);end
end
