function a = lif_read_arrivals(buffer,stepIndex,delaySteps)
%LIF_READ_ARRIVALS Read spikes at t-delay before writing future spikes.
L=size(buffer,1);n=size(buffer,2);assert(L>max(delaySteps,[],'all'));
row=mod(stepIndex-delaySteps,L)+1;
idx=row+L*repmat(0:n-1,n,1);a=buffer(idx);a(1:n+1:end)=false;
end
