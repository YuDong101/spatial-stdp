function dw = lif_somatic_update(p,lastSpike,spiking,kernel)
%LIF_SOMATIC_UPDATE Nearest recorded somatic spikes; simultaneous pairs are zero.
n=numel(lastSpike);delta=lastSpike(:)-lastSpike(:)';
eligible=isfinite(delta) & abs(delta)<=p.stdpCutoffMs & (spiking(:)|spiking(:)');
dw=zeros(n);positive=eligible & delta>0;negative=eligible & delta<0;
dw(positive)=p.Aplus*exp(-delta(positive)/p.tauPlusMs);
dw(negative)=-p.Aminus*exp(delta(negative)/p.tauMinusMs);
dw=dw.*kernel*p.learningScale;dw(1:n+1:end)=0;
end
