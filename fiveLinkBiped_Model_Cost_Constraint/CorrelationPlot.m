% Produce the input lower triangular matrix data

function CorrelationPlot(feat_matrix, feat_names)
% C = -1 + 2.*rand(12,12);
% C = tril(C,-1);
% C(logical(eye(size(C)))) = 1;
CorrColormap = customcolormap(linspace(0,1,11), {'#68011d','#b5172f','#d75f4e','#f7a580','#fedbc9','#f5f9f3','#d5e2f0','#93c5dc','#4295c1','#2265ad','#062e61'});


len_cost = size(feat_matrix,2);
for i = 1:len_cost
    for j = i:len_cost
        
        C(i,j) = corr2(feat_matrix(:,i),feat_matrix(:,j));
        
    end
end

% C = reshape(C',1,len_cost*len_cost);
C = transpose(C);
% Set [min,max] value of C to scale colors
clrLim = [-1,1];
% load('CorrColormap.mat') % Uncomment for custom CorrColormap
% Set the  [min,max] of diameter where 1 consumes entire grid square
diamLim = [0.1, 1];
% myLabel = {'ICA','Elev','Pr','Rmax','Rmin','Srad','Wspd','Tmin','Tmax','VPD','ET_o','AW'};
myLabel = feat_names;
% Compute center of each circle
% This assumes the x and y values were not entered in imagesc()
x = 1 : 1 : size(C,2); % x edges
y = 1 : 1 : size(C,1); % y edges
[xAll, yAll] = meshgrid(x,y);
xAll(C==0)=nan; % eliminate cordinates for zero correlations
% Set color of each rectangle
% Set color scale
% cmap = winter(256);
cmap = CorrColormap; % Uncomment for CorrColormap
Cscaled = (C - clrLim(1))/range(clrLim); % always [0:1]
colIdx = discretize(Cscaled,linspace(0,1,size(cmap,1)));
% Set size of each circle
% Scale the size between [0 1]
Cscaled = (abs(C)*.9 - 0)/1;
diamSize = Cscaled * range(diamLim) + diamLim(1);
% Create figure
fh = figure(321);
% subplot(1,2,1)
% set(gcf,'Color','white','renderer','Painters','Units','Centimeters','Position',[17 1 17 12]);
ax = axes(fh);
hold(ax,'on')
% colormap(ax,'winter');
colormap(CorrColormap) %Uncomment for CorrColormap
tickvalues = 1:length(C);
x = zeros(size(tickvalues));
text(x, tickvalues, myLabel, 'HorizontalAlignment', 'right', 'FontSize', 15);
x(:) = length(C)+1;
text(tickvalues, x, myLabel, 'HorizontalAlignment', 'right','Rotation',90,  'FontSize', 15);
% text(tickvalues+.3, x, myLabel, 'HorizontalAlignment', 'right',  'FontSize', 15);
% Create circles
theta = linspace(0,2*pi,50); % the smaller, the less memory req'd.
h = arrayfun(@(i)fill(diamSize(i)/2 * cos(theta) + xAll(i), ...
    diamSize(i)/2 * sin(theta) + yAll(i), cmap(colIdx(i),:),'LineStyle','none'),1:numel(xAll));
axis(ax,'equal')
axis(ax,'tight')
set(ax,'YDir','Reverse')
colorbar('Ticks',[-1 -.5 0 .5, 1],...
         'TickLabels',{-1 -.5 0 .5, 1})
caxis(clrLim);
axis off
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.002 .002] , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.7       ,...
  'FontSize'    , 16);

