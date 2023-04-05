clc 
clearvars

%% Load relevant repository MATLAB functions

addpath(genpath('./BrainEigenModes/functions_matlab'));
addpath('./fcn/')

%% Load surface files for visualization

surface_interest = 'fsLR_32k';
hemisphere = 'lh';
mesh_interest = 'midthickness';

% Load midthickness
[vertices, faces] = read_vtk(sprintf('./BrainEigenModes/data/template_surfaces_volumes/%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';

% Load cortex mask
cortex = logical(dlmread(sprintf('./BrainEigenModes/data/template_surfaces_volumes/%s_cortex-%s_mask.txt', surface_interest, hemisphere)));

disp('loaded surfaces')

%% Load eigenmodes

hemisphere = 'lh';
num_modes = 50;

if num_modes == 200
    eigenmodes = dlmread(sprintf('./osf_dl/template_eigenmodes/fsLR_32k_midthickness-%s_emode_%i.txt', hemisphere, num_modes));
elseif num_modes == 50
    eigenmodes = dlmread(sprintf('./BrainEigenModes/data/examples/fsLR_32k_midthickness-%s_emode_%i.txt', hemisphere, num_modes));
else
    error('dont have that number of modes, please')
end

disp('loaded eigenmodes')

%% Setup for BrainSpace

addpath(genpath('./BrainSpace/matlab'))

% make a matlab surface that BrainSpace is happy with
brainspace_midthick = struct() ; 
brainspace_midthick.tri = surface_midthickness.faces ;  
brainspace_midthick.coord = surface_midthickness.vertices' ; 

%% Make geodesic

G = surface_to_graph(brainspace_midthick,'mesh',logical(~cortex),true) ; 
shortest_paths = distances(G,'Method','unweighted') ; % computes shortest paths
closeness = 1./shortest_paths ; 
closeness(isinf(closeness)) = 0 ;

dist_thr = prctile(shortest_paths,20,'all') ; % arbitrary param, how much 
                                             % local sp neighbourhood to
                                             % capture in graph...
W = closeness .* (shortest_paths < dist_thr) ; 

clear closeness shortest_paths

% tic
% MEM = compute_mem(geodesic_sim_thr) ; 
% toc

% center it
W = full(W - mean(W) - mean(W,2) + mean(W,"all",'omitnan')); 

assert(issymmetric(W),'not symmetric,problem')

tic
[MEM,lambda] = eigs(W,num_modes*10); % not the full decomp, to save time
                                      % still waaaay less
                                      % dimensionality that the mesh
toc

lambda = diag(lambda) ; 

% Sort eigenvectors and values.
[~, idx] = sort(lambda,'descend');
MEM = MEM(:,idx);

%% Loop over zstat maps and do moran tests

nperms = 5e3; 

loaded_data = load('./BrainEigenmodes/data/figures/Figure1.mat') ;
map_names = fieldnames(loaded_data.task_map_emp) ;

for map_idx = 1:length(map_names)

    moran_results = struct() ; 
    moran_results.name = map_names{map_idx} ; 

    %% Get the empirical accuracy  
    data_to_reconstruct = loaded_data.task_map_emp.(map_names{map_idx}) ; 
    
    % viz it?
    % draw_surface_bluewhitered_dull(surface_midthickness, data_to_reconstruct, hemisphere, find(cortex==0), 1);
    
    %% do spins and recort results
    
    filename = sprintf('./gen_data/moranres_%s_%s-%s.mat',map_names{map_idx},surface_interest,hemisphere) ; 
    
    if ~isfile(filename)

        % compute original accuracy
        moran_results.recon_acc = calc_eigdecomp_recon_acc(data_to_reconstruct(cortex),eigenmodes(cortex,:),num_modes) ; 
        
        % make the randomized data
        surr_data_tmp = moran_randomization(data_to_reconstruct(~~cortex),...
                                MEM,nperms,'procedure','singleton') ;
        surr_data = nan(length(cortex),nperms) ; 
        for idx = 1:nperms
            surr_data(~~cortex,idx) = surr_data_tmp(:,1,idx) ;  
        end
        clear surr_data_tmp

        perm_acc_moran = nan(num_modes,nperms) ; 

        % measure reconstruct accuracy on spin data
        parfor idx = 1:nperms
        
        %     if mod(idx,10) == 1
        %         disp([ num2str(idx) ' of ' num2str(nperms) ])
        %     end
        
            disp(idx)
        
            perm_acc_moran(:,idx) = calc_eigdecomp_recon_acc(surr_data(cortex,idx),...
                                               eigenmodes(cortex,:),num_modes) ; 
        
        end
    
        % record results in the struct 
        moran_results.perm_acc = perm_acc_moran ; 

        % save it
        save(filename,'moran_results')

    end

end

%% viz it

tiledlayout(2,length(map_names))

for idx = 1:length(map_names)

    filename = sprintf('./gen_data/moranres_%s_%s-%s.mat',map_names{idx},surface_interest,hemisphere) ; 

    ll = load(filename) ; 

    perm_acc = ll.moran_results.perm_acc ;
    recon_acc = ll.moran_results.recon_acc ;

    nexttile(idx)
    plot(perm_acc,'Color',[0 0.4470 0.7410 0.05],'LineWidth',2) 
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