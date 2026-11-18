function plot_figure1(alpha, gmax, winner, fraction, phases, outputDir)
% Render the modal map and the six finite-window relative phase trajectories.
% Inputs for b--g are 100 x 1001, ordered as in the manuscript.
palette = load(fullfile(fileparts(mfilename('fullpath')), 'figure1_palette.mat'));
[alphaLevels, ~, column] = unique(alpha(:));
[conductances, ~, row] = unique(gmax(:));
gridSize = [numel(conductances), numel(alphaLevels)];
index = sub2ind(gridSize, row, column);
assert(isequal(gridSize,[50 60]) && numel(unique(index)) == 3000 && numel(index) == 3000);
assert(all(ismember(winner,0:5)) && all(fraction > 0 & fraction <= 1));

% Flat shading takes the lower-left vertex color of each cell.
rgb = ones([gridSize + 1, 3]);
colors = palette.stateCmap(winner(:) + 1,:) .* fraction(:) + (1 - fraction(:));
for channel = 1:3
    layer = zeros(gridSize);
    layer(index) = colors(:,channel);
    rgb(1:end-1,1:end-1,channel) = layer;
end

f = figure('Visible','off','Units','points','Position',[30 30 585 565], ...
    'Color','white','Name','Figure 1');
closeOnExit = onCleanup(@() close(f)); %#ok<NASGU>
map = newAxes(f,[55 210 299.38 299.38]);
xEdges = logEdges(alphaLevels);
yEdges = logEdges(conductances);
surface(map,xEdges,yEdges,zeros(gridSize + 1),rgb,'EdgeColor','none','FaceColor','flat');
set(map,'Tag','StateMap','XScale','log','YScale','log', ...
    'XLim',xEdges([1 end]),'YLim',yEdges([1 end]), ...
    'XTick',[0.1 1 10],'YTick',[0.0001 0.001 0.01 0.1]);
colormap(map,palette.stateCmap);
xlabel(map,'$\alpha$','Interpreter','latex','FontSize',9.6);
ylabel(map,'$g_{\max}$ ($\mu$S)','Interpreter','latex','FontSize',9.6);
panelLabel(map,'a');
labels = 'bcdefg';
examples = [1 .0001; .3 .01; 1 .001; 1.34 .002; 2 .02; 5 .01];
locations = [425 441.30 130 68.076; 425 325.65 130 68.076; ...
    425 210 130 68.076; 55 66 130 71; 240 66 130 71; 425 66 130 71];
for k = 1:6
    trace = phases{k};
    assert(isequal(size(trace),[100 1001]) && all(isfinite(trace),'all'));
    x = examples(k,1);
    y = examples(k,2);
    line(map,x,y,'LineStyle','none','Marker','o','MarkerSize',12, ...
        'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',0.8,'Clipping','off');
    text(map,x,y,labels(k),'FontName','Arial','FontSize',10.5,'FontWeight','bold', ...
        'HorizontalAlignment','center','VerticalAlignment','middle','Clipping','off');
    a = newAxes(f,locations(k,:));
    image(a,'XData',[15.001 20],'YData',[1 100],'CData',trace, ...
        'CDataMapping','scaled','Tag',['phase_' labels(k)]);
    set(a,'CLim',[0 1],'XLim',[14.9985 20.0025],'YLim',[.5 100.5], ...
        'XTick',[15 17.5 20],'YTick',[1 50 100],'Layer','top');
    colormap(a,palette.phaseCmap);
    xlabel(a,'Time (s)','FontSize',9.6);
    ylabel(a,'Node rank','FontSize',9.6);
    text(a,65,locations(k,4)+7,sprintf('{\\it\\alpha} = %g; {\\itg}_{max} = %g',x,y), ...
        'Units','points','FontName','Arial','FontSize',9,'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom');
    panelLabel(a,labels(k));
end
phaseScale = newAxes(f,[175 17 260 5]);
image(phaseScale,'XData',[0 1],'YData',[0 1],'CData',linspace(0,1,257), ...
    'CDataMapping','scaled');
set(phaseScale,'XLim',[-.5/256 1+.5/256],'YLim',[-1 2], ...
    'XTick',[0 .5 1],'YTick',[],'CLim',[0 1],'Layer','top');
colormap(phaseScale,palette.phaseCmap);
text(phaseScale,130,12,'Relative phase (cycle)','Units','points', ...
    'FontName','Arial','FontSize',9,'HorizontalAlignment','center','VerticalAlignment','bottom');

key = axes(f,'Units','points','Position',[60 218 134 53], ...
    'XLim',[0 134],'YLim',[0 53],'Visible','off','Tag','StateKey');
rectangle(key,'Position',[0 0 134 53],'FaceColor','w','EdgeColor','none');
states = {'Dispersion without full-cycle slips','Phase clusters','Partial locking', ...
    'Locking-slip coexistence','Dispersed slipping'};
for k = 1:5
    y = 57 - 10*k;
    rectangle(key,'Position',[4 y-2.3 7.5 4.6],'FaceColor',palette.stateCmap(k+1,:), ...
        'EdgeColor','none');
    text(key,15,y,states{k},'FontName','Arial','FontSize',7.7,'VerticalAlignment','middle');
end
drawnow;
savefig(f,fullfile(outputDir,'Figure1.fig'));
exportgraphics(f,fullfile(outputDir,'Figure1.pdf'),'ContentType','vector','BackgroundColor','white');
exportgraphics(f,fullfile(outputDir,'Figure1.png'),'Resolution',200,'BackgroundColor','white');
end

function a = newAxes(f,position)
a = axes(f,'Units','points','Position',position,'FontName','Arial','FontSize',9, ...
    'LineWidth',.65,'Box','on','TickDir','out','NextPlot','add');
end

function edges = logEdges(centers)
edges = exp(interp1(1:numel(centers),log(centers),.5:numel(centers)+.5,'linear','extrap'));
end

function panelLabel(a,letter)
text(a,-29,a.Position(4)+10,['(' letter ')'],'Units','points', ...
    'FontName','Arial','FontSize',10.5,'FontWeight','bold','VerticalAlignment','bottom');
end
