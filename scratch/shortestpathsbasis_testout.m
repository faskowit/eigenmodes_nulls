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
brainspace_midthick = struct() ; 
brainspace_midthick.tri = surface_midthickness.faces ;  
brainspace_midthick.coord = surface_midthickness.vertices' ; 

%% Make shortest paths basis

sp_levels = [ 1 5:5:50 ] ; 

filename = sprintf('./gen_data/spbasis_%s_%s-%s.mat',num2str(num_modes),surface_interest,hemisphere) ; 

if ~isfile(filename)

    sp_basis = struct() ; 

    G = surface_to_graph(brainspace_midthick,'mesh',logical(~cortex),true) ; 

    shortest_paths = distances(G,'Method','unweighted') ; % computes shortest paths
    
    closeness = 1./shortest_paths ; 
    closeness(isinf(closeness)) = 0 ;
    
    sp_thrs = prctile(shortest_paths,sp_levels,'all') ; % arbitrary param, how much 
                                                     % local sp neighbourhood to
                                                     % capture in graph...

    for jdx = 1:length(sp_levels)
        
        G_mask = (shortest_paths < sp_thrs(jdx)) ;  
    
        W = closeness .* G_mask ; 
                
        % center it
        W = full(W - mean(W) - mean(W,2) + mean(W,"all",'omitnan')); 
        
        % enforce symmetry
        W = triu(W,1) + triu(W)' ; 
        assert(issymmetric(W),'not symmetric,problem')
        
        disp([ 'decomp ' num2str(jdx) '- out of  ' num2str(length(sp_levels))])
        tic
        [VV,DD] = eigs(W,num_modes); 
                                                                         
        toc
        disp('finished decomp')
    
        % Sort eigenvectors and values.
        [DDsort, ss] = sort(diag(DD),'descend');
        VV = VV(:,ss);
    
        sp_basis(jdx).VV = VV ;
        sp_basis(jdx).DD = DDsort ; 
        sp_basis(jdx).thr = sp_thrs(jdx) ;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    save(filename,'sp_basis')
else
    load(filename)
end

clear closeness shortest_paths

%% Loop over zstat maps and do moran tests

filename = sprintf('./gen_data/spbasis_pred_%s_%s-%s.mat',num2str(num_modes),surface_interest,hemisphere) ; 

if ~isfile(filename)

    loaded_data = load('./NSBLab_repo/data/figures_Nature/Figure1.mat') ;
    map_names = fieldnames(loaded_data.task_map_emp) ;

    spbasis_results = struct() ; 

    for map_idx = 1:length(map_names)
    
        spbasis_results(map_idx).name = map_names{map_idx} ; 
    
        %% Get the empirical accuracy  
        data_to_reconstruct = loaded_data.task_map_emp.(map_names{map_idx}) ; 
        
        % viz it?
        % draw_surface_bluewhitered_dull(surface_midthickness, data_to_reconstruct, hemisphere, find(cortex==0), 1);
        
        %% do spins and recort results
        
        % compute original accuracy
        spbasis_results(map_idx).emp = calc_eigdecomp_recon_acc(data_to_reconstruct(cortex),eigenmodes(cortex,:),num_modes) ; 
    
        sp_res = nan(num_modes,length(sp_levels)) ; 
        dtr = data_to_reconstruct(cortex) ;

        parfor idx = 1:length(sp_levels)
    
            disp([num2str(map_idx) ' ' num2str(idx)])
    
            sp_res(:,idx) = calc_eigdecomp_recon_acc(dtr,...
                sp_basis(idx).VV,num_modes) ; 
        
        end
    
        spbasis_results(map_idx).spres = sp_res ; 
    
    end
    
    save(filename,"spbasis_results")

else
    load(filename)
end

%% viz it

tiledlayout(2,length(map_names))

for idx = 1:length(map_names)

    %filename = sprintf('./gen_data/moranres_%s_%s-%s.mat',map_names{idx},surface_interest,hemisphere) ; 

    %ll = load(filename) ; 

    perm_acc = spbasis_results(idx).spres ;
    recon_acc = spbasis_results(idx).emp ; 

    nexttile(idx)
    plot_manylines(perm_acc,'Color',[0 0.4470 0.7410 0.9],'LineWidth',2) ; 
    hold on 
    plot(recon_acc,'r','LineWidth',2)
    hold off
    
    xlim([1 num_modes])
    ylim([-0.25 1])

    if idx == 1
        ylabel('recon accuracy')
    end

    title(map_names{idx},'Interpreter','none')
    
    nexttile(idx+length(map_names))
    pvals = ( sum(bsxfun(@gt,perm_acc,recon_acc'),2) + 1) ./ (size(perm_acc,2)+1) ; 
    plot(pvals,'mo-','LineWidth',2)
    
    xlim([1 num_modes])
    ylim([0 1])

    % title([ 'p-value' map_names{idx}],'Interpreter','none')
    
    if idx == 1
        ylabel('p-value')
    end

    if idx == 4
        xlabel('eigenmodes used for recon')
    end

end