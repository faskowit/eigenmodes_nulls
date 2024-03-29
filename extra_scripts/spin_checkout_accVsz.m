clc 
clearvars

%% Load relevant repository MATLAB functions

addpath(genpath('./NSBLab_repo/functions_matlab'));
addpath('./fcn/')

%% Load surface files for visualization

surface_interest = 'fsLR_32k';
hemisphere = 'lh';
mesh_interest = 'midthickness';

% Load midthickness
[vertices, faces] = read_vtk(sprintf('./NSBLab_repo/data/template_surfaces_volumes/%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';

% Load sphere
mesh_interest = 'sphere';
[vertices, faces] = read_vtk(sprintf('./NSBLab_repo/data/template_surfaces_volumes/%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere));
surface_sphere.vertices = vertices';
surface_sphere.faces = faces';

% Load cortex mask
cortex = logical(dlmread(sprintf('./NSBLab_repo/data/template_surfaces_volumes/%s_cortex-%s_mask.txt', surface_interest, hemisphere)));

disp('loaded surfaces')

%% Load connectome modes

hemisphere = 'lh';
num_modes = 200;

if num_modes == 200
    eigenmodes = dlmread(sprintf('./osf_dl/template_eigenmodes/fsLR_32k_midthickness-%s_emode_%i.txt', hemisphere, num_modes));
elseif num_modes == 50
    eigenmodes = dlmread(sprintf('./NSBLab_repo/data/examples/fsLR_32k_midthickness-%s_emode_%i.txt', hemisphere, num_modes));
else
    error('dont have that number of modes, please')
end

disp('loaded eigenmodes')

%% Setup for BrainSpace

addpath(genpath('./BrainSpace/matlab'))

% make a matlab surface that BrainSpace is happy with
brainspace_sphere = struct() ; 
brainspace_sphere.tri = surface_sphere.faces ;  
brainspace_sphere.coord = surface_sphere.vertices' ; 

%% Generate spin inds, which we can use for all of the tests

ncoords = size(brainspace_sphere.coord,2) ; 

filename = sprintf('./gen_data/spininds_%s-%s.mat',surface_interest,hemisphere) ; 

if ~isfile(filename)
    disp('generating spin inds')
    spin_inds = spin_permutations((1:ncoords)',...
        brainspace_sphere,2.5e4,...
        'random_state',42) ;
    assert(ncoords<intmax('int16'))
    spin_inds = int16(squeeze(spin_inds{1})) ; 
    save(filename,'spin_inds')
    disp('finished making spin nulls')
else
    load(filename)
    disp('loaded spin inds')
end

%% Loop over zstat maps and do spin tests

nperms = 5e3; 
% nperms = 1000 ;  

trpvals = ceil(logspace(log10(2),log10(200),10)) ; 
trpvals(end) = 200 ; 

loaded_data = load('./NSBLab_repo/data/figures_Nature/Figure1.mat') ;
map_names = fieldnames(loaded_data.task_map_emp) ;

perm_acc_v_sz = struct() ; 

for map_idx = 1:length(map_names)

    %% Get the empirical accuracy  
    data_to_reconstruct = loaded_data.task_map_emp.(map_names{map_idx}) ; 
    
    % make a smaller version to send to parfor
    spin_inds_smaller = spin_inds(:,randperm(size(spin_inds,2),nperms)) ; 
    
    perm_acc_n_mwsize = nan(nperms,length(trpvals)+1) ; 

    % measure reconstruct accuracy on spin data
    parfor idx = 1:nperms
    
        disp([num2str(map_idx) ' - ' num2str(idx)])
    
        tmp_surr_data = data_to_reconstruct(spin_inds_smaller(:,idx)) ;
    
        % dont predict where the medial wall has 'spun' into
        perm_mask = (cortex & ~isnan(tmp_surr_data)) ; 

        tmp = calc_eigdecomp_recon_acc(tmp_surr_data(perm_mask),...
                                           eigenmodes(perm_mask,:),num_modes) ; 
        
        perm_acc_n_mwsize(idx,:) = [ arrayfun(@(x_) trapz(tmp(1:x_)),trpvals) ...
            sum(perm_mask,'all') ] ;

    end

    perm_acc_v_sz(map_idx).acc = perm_acc_n_mwsize(:,1:length(trpvals)) ;
    perm_acc_v_sz(map_idx).sz = perm_acc_n_mwsize(:,length(trpvals)+1) ;

end

%%

filename = sprintf('./gen_data/spin_accVsz_%s-%s.mat',surface_interest,hemisphere) ; 
% and then save it
save(filename,"perm_acc_v_sz")

%% viz it

tiledlayout(3,length(map_names))
set(gcf,'Position', [200 200 1300 600]);

map_names_better = { 'Social' 'Motor' 'Gambling' 'WM' ...
    'Language' 'Emotion' 'Relational' } ;

floorsz = length(cortex) - (sum(~cortex)*2) ; 
maxsz = sum(cortex) ; 

cmap = flipud(winter(10)) ; 
piktrap = [ 4 9 10 ] ; 

for idx = 1:length(map_names)

    for jdx = 1:length(piktrap) 

        xx = idx + ((jdx-1)*length(map_names))  
        nexttile(xx)


        pp = piktrap(jdx) ; 

        x = perm_acc_v_sz(idx).acc(:,pp) ; 
        y = perm_acc_v_sz(idx).sz ; 
    
        [s,tt] = scatter_w_rho(x,y,'filled','o','MarkerFaceColor',cmap(jdx,:)) ; 
        tt.FontSize = 7 ; 
        tt.Position = [0.015,0.05] ;
        s.MarkerFaceAlpha = 0.2 ;

%         h = histogram2(x(:),y(:),[10 10],'DisplayStyle','tile');
%         h.EdgeColor = 'none' ; 
%         grid off
        
        yline(floorsz,'Color','r','LineWidth',2)
        yline(maxsz,'Color','c','LineWidth',2)

        ylim([26.25e3 30e3])

        if jdx == 1
            title( {map_names_better{idx} ' ' } )
        end

        if idx == 4 
            xlabel([ 'Recon acc. area after ' num2str(trpvals(pp)) ' modes'])
        end

        if idx == 1
            ylabel('Num. vertices in recon.')
        end

    end
    
end

%%

outfile = './figures/spin_checkout_accVsz.pdf' ; 
orient(gcf,'landscape')
print(gcf,'-dpdf',outfile,'-vector','-bestfit')

outfile = './figures/spin_checkout_accVsz.png' ; 
orient(gcf,'landscape')
print(gcf,'-dpng',outfile)