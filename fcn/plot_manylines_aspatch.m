function [ h , pp_spans , cmap ] = plot_manylines_aspatch(inlines,varargin)
% inlines should be matrix, where each column is data you want to plot

pp = 0:25 ; 
percentiles = flipud([ (pp)' 100-(pp)' ]) ;
pp_spans = percentiles(:,2) - percentiles(:,1) ; 

x = 1:size(inlines,1) ;

% make the low density more visibile by making it ~10%
lowbump = floor(size(percentiles,1)*0.1) ; 

cmap = [1 1 1 ; 0 0.447 0.741 ]; % white to blue
cmap = interp1(linspace(0,1,size(cmap,1)),cmap, ...
    linspace(0,1,size(percentiles,1)+1+lowbump)); % +1 so that the end isn't white
cmap = flipud(cmap) ; 
cmap = cmap(1:size(percentiles,1),:) ; 

for idx = fliplr(1:length(percentiles))
    % disp(idx)
    y = prctile(inlines',percentiles(idx,:)) ; 

    hold on
%     fill([x fliplr(x)], [y(1,:) fliplr(y(2,:))], [0 0.447 0.741] , 'facealpha', 0.1, 'edgealpha', 0 )
    fill([x fliplr(x)], [y(1,:) fliplr(y(2,:))], cmap(idx,:), 'edgealpha', 0 )

end

hold off

h = gca ;

