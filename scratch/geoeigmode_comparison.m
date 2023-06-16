clc 
clearvars

%% Load relevant repository MATLAB functions

addpath(genpath('./NSBLab_repo/functions_matlab'));
addpath('./fcn/')

%% Load surface files for visualization

surface_interest = 'fsLR_32k';
hemisphere = 'lh';
mesh_interest = { 'sphere' 'midthickness' } ;

% Load cortex mask
cortex = logical(dlmread(sprintf('./NSBLab_repo/data/template_surfaces_volumes/%s_cortex-%s_mask.txt', surface_interest, hemisphere)));

disp('loaded surfaces')

%% First load paper eigenmodes

hemisphere = 'lh';
num_modes = 200;

modes_str = struct() ; 

for idx = 1:length(mesh_interest)

    disp([ 'loading: ' mesh_interest{idx} ])

    switch mesh_interest{idx}
        case { 'sphere' 'veryinflated' 'pial' 'white' }
            eigenmodes = dlmread(sprintf('./gen_data/fsLR_32k_%s-%s_emode_%i.txt', mesh_interest{idx}, hemisphere, num_modes));
        otherwise
            eigenmodes = dlmread(sprintf('./osf_dl/template_eigenmodes/fsLR_32k_%s-%s_emode_%i.txt', mesh_interest{idx}, hemisphere, num_modes));
    end

    modes_str.(mesh_interest{idx}).eigenmodes = eigenmodes ; 

end

disp('loaded eigenmodes')

%% how to conpare the sets of eigenmodes

pbins = 10 ;

disc_mthick = cell2mat(arrayfun(@(x_) ...
    discretize(modes_str.midthickness.eigenmodes(cortex,x_),pbins),1:num_modes,UniformOutput=false)) ;

disc_sphere = cell2mat(arrayfun(@(x_) ...
    discretize(modes_str.sphere.eigenmodes(cortex,x_),pbins),1:num_modes,UniformOutput=false)) ;

%%

costMat = 1-abs(corr(modes_str.midthickness.eigenmodes, ...
                modes_str.sphere.eigenmodes,'type','p')) ; 


% force the trivial first mode to match
costMat(1,:) = inf(size(costMat,1),1) ; 
costMat(:,1) = inf(size(costMat,1),1) ; 
costMat(1,1) = 0 ; 

[a,b] = munkres(costMat) ; 

% subplot(1,2,1)
% imagesc(costMat)
% subplot(1,2,2)
% imagesc(costMat(a,a))

modes_str.sphere.match_em = modes_str.sphere.eigenmodes(:,a) ; 

modes_corr_emp = abs(arrayfun(@(i_) corr( ...
    modes_str.sphere.match_em(cortex,i_),...
    modes_str.midthickness.eigenmodes(cortex,i_)), 1:num_modes )) ; 

%%

filename = sprintf('./gen_data/spininds_%s-%s.mat',surface_interest,hemisphere) ; 
ll = load(filename) ; 

nperms = 5000 ; 

rng(4242)
spin_inds_smaller = ll.spin_inds(:,randperm(size(ll.spin_inds,2),nperms)) ; 

%% do some spin testing

filename = sprintf('./gen_data/sphere-midthick_comp_%s-%s.mat',surface_interest,hemisphere) ; 

if ~isfile(filename)

    modes_corr_rand = nan(num_modes,nperms) ; 
    for idx = 1:nperms
    
        disp(idx)
    
        spind = spin_inds_smaller(:,idx) ; 
    
        % only compare in corex & not spun medial wall area
        spmask = (cortex) & (cortex(spind)) ; 
            
        spundat = modes_str.sphere.match_em(spind,:) ; 
    
        modes_corr_rand(:,idx) = arrayfun(@(i_) corr( ...
            spundat(spmask,i_),...
            modes_str.midthickness.eigenmodes(spmask,i_)), 1:num_modes ) ; 
    
    end

    save(filename,'modes_corr_rand')
else
    load(filename)
end

%%

tiledlayout(2,1)

nexttile
plot_manylines_aspatch(abs(modes_corr_rand))
hold on
plot(abs(modes_corr_emp),'Color','r')
hold off

xlim([2 200])
ylim([0 1])

title('emp. spatial corr (red), spun spatial corr (blue)')

nexttile

pvals = ( sum(bsxfun(@gt,abs(modes_corr_rand),modes_corr_emp'),2) + 1) ./ (nperms+1) ; 

plot(pvals,'Color','r')
xlim([2 200])
ylim([0 1])
title('p-val')

