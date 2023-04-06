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

%% Load eigenmodes

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

loaded_data = load('./NSBLab_repo/data/figures_Nature/Figure1.mat') ;
map_names = fieldnames(loaded_data.task_map_emp) ;

for map_idx = 1:length(map_names)

    spin_results = struct() ; 
    spin_results.name = map_names{map_idx} ; 

    %% Get the empirical accuracy  
    data_to_reconstruct = loaded_data.task_map_emp.(map_names{map_idx}) ; 
    
    % viz it?
    % draw_surface_bluewhitered_dull(surface_midthickness, data_to_reconstruct, hemisphere, find(cortex==0), 1);
    
    %% do spins and recort results
    
    filename = sprintf('./gen_data/spinres_%s_%s-%s.mat',map_names{map_idx},surface_interest,hemisphere) ; 
    
    if ~isfile(filename)

        % compute original accuracy
        spin_results.recon_acc = calc_eigdecomp_recon_acc(data_to_reconstruct(cortex),eigenmodes(cortex,:),num_modes) ; 
        
        % make a smaller version to send to parfor
        spin_inds_smaller = spin_inds(:,randperm(5e4,nperms)) ; 
        
        perm_acc_spin = nan(num_modes,nperms) ; 

        % measure reconstruct accuracy on spin data
        parfor idx = 1:nperms
        
        %     if mod(idx,10) == 1
        %         disp([ num2str(idx) ' of ' num2str(nperms) ])
        %     end
        
            disp(idx)
        
            tmp_surr_data = data_to_reconstruct(spin_inds_smaller(:,idx)) ;
        
            % dont predict where the medial wall has 'spun' into
            perm_mask = (cortex & ~isnan(tmp_surr_data)) ; 
        
            perm_acc_spin(:,idx) = calc_eigdecomp_recon_acc(tmp_surr_data(perm_mask),...
                                               eigenmodes(perm_mask,:),num_modes) ; 
        
        end
    
        % record results in the struct 
        spin_results.perm_acc = perm_acc_spin ; 

        % save it
        save(filename,'spin_results')

    end

end

%% viz it

tiledlayout(2,length(map_names))

for idx = 1:length(map_names)

    filename = sprintf('./gen_data/spinres_%s_%s-%s.mat',map_names{idx},surface_interest,hemisphere) ; 

    ll = load(filename) ; 

    perm_acc = ll.spin_results.perm_acc ;
    recon_acc = ll.spin_results.recon_acc ;

    nexttile(idx)
    plot_manylines(perm_acc,'Color',[0 0.4470 0.7410 0.05],'LineWidth',2) 
    hold on 
    plot(recon_acc,'r','LineWidth',2)
    hold off
    
    xlim([1 50])
    ylim([-0.25 1])

    if idx == 1
        ylabel('recon accuracy')
    end

    title(map_names{idx},'Interpreter','none')
    
    nexttile(idx+length(map_names))
    pvals = ( sum(bsxfun(@gt,perm_acc,recon_acc'),2) + 1) ./ (nperms+1) ; 
    plot(pvals,'mo-','LineWidth',2)
    
    xlim([1 50])
    ylim([0 1])

    % title([ 'p-value' map_names{idx}],'Interpreter','none')
    
    if idx == 1
        ylabel('p-value')
    end

    if idx == 4
        xlabel('eigenmodes used for recon')
    end

end




