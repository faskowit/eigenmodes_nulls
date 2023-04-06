function [ h ] = plot_manylines(inlines,varargin)
% inlines should be matrix, where each column is data you want to plot

x = repmat((1:size(inlines,1))',1,size(inlines,2)) ; 
x =[ x ; nan(1,size(inlines,2)) ] ; 

y = [ inlines ; nan(1,size(inlines,2)) ] ; 

% make inlines into one line, then plot eet
h = plot(x(:),y(:),varargin{:}) ;  
