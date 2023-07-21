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

%% Setup for BrainSpace

addpath(genpath('./BrainSpace/matlab'))

% make a matlab surface that BrainSpace is happy with
brainspace_midthick = struct() ; 
brainspace_midthick.tri = surface_midthickness.faces ;  
brainspace_midthick.coord = surface_midthickness.vertices' ; 

%% Make MEM the BrainSpace way

hemisphere = 'lh';
num_modes = 200;

tic
filename = sprintf('./gen_data/moranMEM_BS_modes-%s_%s-%s.mat',num2str(num_modes),surface_interest,hemisphere) ; 
if ~isfile(filename)
    disp('computing MEM')
    datestr(now,'HH:MM:SS')
    MEM = compute_mem(brainspace_midthick,...
        'mask',~cortex,'n_ring',1e4,...
        'distance_metric','geodesic') ; 
    save(filename,"MEM")
else
    disp('loading MEM')
    load(filename)
end
toc
disp('finished')

%% Load eigenmodes

if num_modes == 200
    eigenmodes = dlmread(sprintf('./osf_dl/template_eigenmodes/fsLR_32k_midthickness-%s_emode_%i.txt', hemisphere, num_modes));
elseif num_modes == 50
    eigenmodes = dlmread(sprintf('./NSBLab_repo/data/examples/fsLR_32k_midthickness-%s_emode_%i.txt', hemisphere, num_modes));
else
    error('dont have that number of modes, please')
end

disp('loaded eigenmodes')

%% Make geodesic
% 
% hemisphere = 'lh';
% num_modes = 200;
% 
% filename = sprintf('./gen_data/moranMEM_modes-%s_%s-%s.mat',num2str(num_modes),surface_interest,hemisphere) ; 
% 
% if ~isfile(filename)
% 
%     G = surface_to_graph(brainspace_midthick,'mesh',logical(~cortex),true) ; 
%     shortest_paths = distances(G,'Method','unweighted') ; % computes shortest paths
%     closeness = 1./shortest_paths ; 
%     closeness(isinf(closeness)) = 0 ;
%     
%     dist_thr = prctile(shortest_paths,20,'all') ; % arbitrary param, how much 
%                                                  % local sp neighbourhood to
%                                                  % capture in graph...
%     W = closeness .* (shortest_paths < dist_thr) ; 
%     
%     clear closeness shortest_paths
%     
%     % center it
%     W = full(W - mean(W) - mean(W,2) + mean(W,"all",'omitnan')); 
%     
%     % enforce symmetry
%     W = triu(W,1) + triu(W)' ; 
%     assert(issymmetric(W),'not symmetric,problem')
%     
%     eigsize = 'full' ; 
% 
%     disp('doing decomp')
% 
%     if strcmp(eigsize,'partial')
% 
%         tic
%         [MEM,lambda] = eigs(W,num_modes*10); % not the full decomp, to save time
%                                               % still waaaay less
%                                               % dimensionality that the mesh
%         toc
%         lambda = diag(lambda) ; 
%     
%     elseif strcmp(eigsize,'full')
%         
%         tic
%         % Full eigenvalue decomposition of W. 
%         % [MEM,lambda] = eig(full(W),'vector');
%         W = sparse(W) ; 
%         [MEM,lambda] = eigs(W,size(W,1)) ; 
%         lambda = diag(lambda) ;
%         toc
% 
%     else 
%         error('neeed valid eigsize var')
%     end
% 
%     % Remove zero eigenvector
%     zdx = find(abs(lambda) < 1e-10); 
% 
%     % See supplemental info 3 of Ref 1, function scores.listw().
%     w = [ones(size(MEM,1),1),MEM(:,zdx)];
%     Q = qr(w);
%     MEM(:,zdx) = Q(:,1:end-1);
%     MEM(:,zdx(1)) = [];
%     lambda(zdx(1)) = []; 
% 
%     disp('finished decomp')
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     % Sort eigenvectors and values.
%     [~, ss] = sort(lambda,'descend');
%     MEM = MEM(:,ss);
% 
%     save(filename,'MEM')
%     clear W
% else
%     load(filename)
% end

%%

% Loop over zstat maps and do moran tests

% nperms = 1e3 ; 
nperms = 5e3; 
% nperms  = 100 ; 

loaded_data = load('./NSBLab_repo/data/figures_Nature/Figure1.mat') ;
map_names = fieldnames(loaded_data.task_map_emp) ;

eigenmodes_nocort = eigenmodes(cortex,:) ;

clear eigenmodes

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
        moran_results.recon_acc = calc_eigdecomp_recon_acc(data_to_reconstruct(cortex),eigenmodes_nocort,num_modes) ; 
        
        % make the randomized data
        rng(map_idx)
        disp('generating moran dat')
        surr_data_tmp = moran_randomization(data_to_reconstruct(~~cortex),...
                                MEM,nperms,'procedure','singleton') ;
        surr_data = nan(length(cortex),nperms) ; 
        for idx = 1:nperms
            surr_data(~~cortex,idx) = surr_data_tmp(:,1,idx) ;  
        end
        clear surr_data_tmp

        perm_acc_moran = nan(num_modes,nperms) ; 

        surr_data_nocort = surr_data(cortex,:) ;

        clear surr_data

        % measure reconstruct accuracy on spin data
        parfor idx = 1:nperms
        
        %     if mod(idx,10) == 1
        %         disp([ num2str(idx) ' of ' num2str(nperms) ])
        %     end
        
            disp([ num2str(map_idx) ' -     ' num2str(idx) ] )
        
            perm_acc_moran(:,idx) = calc_eigdecomp_recon_acc(surr_data_nocort(:,idx),...
                                               eigenmodes_nocort,num_modes) ; 
        
        end
    
        % record results in the struct 
        moran_results.perm_acc = perm_acc_moran ; 

        % moran_results.surr_dat = surr_data_nocort ;

        % save it
        save(filename,'moran_results','surr_data_nocort')

    end

end

% %% viz it
% 
% tiledlayout(2,length(map_names))
% 
% for idx = 1:length(map_names)
% 
%     filename = sprintf('./gen_data/moranres_%s_%s-%s.mat',map_names{idx},surface_interest,hemisphere) ; 
% 
%     ll = load(filename) ; 
% 
%     perm_acc = ll.moran_results.perm_acc ;
%     recon_acc = ll.moran_results.recon_acc ;
% 
%     nexttile(idx)
%     plot_manylines(perm_acc,'Color',[0 0.4470 0.7410 0.05],'LineWidth',2) ; 
%     hold on 
%     plot(recon_acc,'r','LineWidth',2)
%     hold off
%     
%     xlim([1 num_modes])
%     ylim([-0.25 1])
% 
%     if idx == 1
%         ylabel('recon accuracy')
%     end
% 
%     title(map_names{idx},'Interpreter','none')
%     
%     nexttile(idx+length(map_names))
%     pvals = ( sum(bsxfun(@gt,perm_acc,recon_acc'),2) + 1) ./ (nperms+1) ; 
%     plot(pvals,'mo-','LineWidth',2)
%     
%     xlim([1 num_modes])
%     ylim([0 1])
% 
%     % title([ 'p-value' map_names{idx}],'Interpreter','none')
%     
%     if idx == 1
%         ylabel('p-value')
%     end
% 
%     if idx == 4
%         xlabel('eigenmodes used for recon')
%     end
% 
% end

%%

surface_interest = 'fsLR_32k';
hemisphere = 'lh';
num_modes = 200 ; 

loaded_data = load('./NSBLab_repo/data/figures_Nature/Figure1.mat') ;
map_names = fieldnames(loaded_data.task_map_emp) ;

tiledlayout(2,length(map_names))
set(gcf,'Position', [200 200 1200 600]);

map_names_better = { 'Social' 'Motor' 'Gambling' 'WM' ...
    'Language' 'Emotion' 'Relational' } ;

for idx = 1:length(map_names)

    filename = sprintf('./gen_data/moranres_%s_%s-%s.mat',map_names{idx},surface_interest,hemisphere) ; 

    ll = load(filename,'moran_results') ; 

    perm_acc = ll.moran_results.perm_acc ;
    recon_acc = ll.moran_results.recon_acc ;

    nexttile(idx)
    [h, ppspan, cmap] = plot_manylines_aspatch(perm_acc) ; 
    hold on 
    plot(recon_acc,'r','LineWidth',1)
    hold off
       
    xlim([1 num_modes])
    ylim([-0.25 1])

    h.XTick = [ 1 h.XTick] ; 

    if idx == 1
        ylabel('Reconstruction accuracy')
    end

    if idx == 4
        xlabel('Geometric eigenmodes used for reconstruction')
    end

    title(map_names_better{idx},'Interpreter','none')
    
    nexttile(idx+length(map_names))
    pvals = ( sum(bsxfun(@gt,perm_acc,recon_acc'),2) + 1) ./ (size(perm_acc,2)+1) ; 
    plot(pvals,'.-','Color',[0.6350 0.0780 0.1840]	,'LineWidth',1,'MarkerSize',2) ; 
    
    xlim([1 num_modes])
    ylim([0 1])

    h = gca ;
    h.XTick = [ 1 h.XTick] ; 

    % title([ 'p-value' map_names{idx}],'Interpreter','none')

    % plot sig?
    sigpoints = (pvals < 0.05) ;
    % disp([ num2str(sum(sigpoints)) ' low of ' num2str(min(pvals))])
    hold on 
    h = scatter(find(sigpoints),ones(sum(sigpoints),1).*.98, ...
        10,[0.8500 0.3250 0.0980],'filled','diamond') ; 
    hold off
    
    % put a lil text of how many sig
    text(204,0.98,num2str(sum(sigpoints)),'FontSize',8,'Color',[0.8500 0.3250 0.0980])

    if idx == 1
        ylabel('p-value')
    end

    if idx == 4
        xlabel('Geometric eigenmodes used for reconstruction')
    end

end

%% save figure 

outfile = './figures/moranres_eigenmodes_patchversion.pdf' ; 
orient(gcf,'landscape')
print(gcf,'-dpdf',outfile,'-bestfit','-vector')
close(gcf)


%% example moran data
% 
% G = surface_to_graph(brainspace_midthick,'geodesic',~cortex,1) ; 
% ddd = data_to_reconstruct(cortex) ; 
% targ_std = mean(arrayfun(@(x_) std(ddd(G.neighbors(x_))),1:size(G.Nodes,1))) ;

% h = quick_trisurf(surface_midthickness,data_to_reconstruct) ; 
% h.EdgeColor = 'none'
% view([-45 12 -5])

%% look at some stuff

% tiledlayout(3,3)
% set(gcf,'Position', [200 200 800 800]);
% 
% %figure
% for idx = 42:50
%     nexttile
% 
%     disp(idx)
% 
%     ss = nan(length(cortex),1) ;
%     s = eigenmodes_nocort(:,idx) ; 
%     ss(cortex) = s ; 
%     h = quick_trisurf(surface_midthickness,ss)
%     h.EdgeColor = 'none' ; 
%     view([-75 11 -5])
%     material shiny
%     % camlight headlight
%     lighting gouraud
%     xticks('') ; yticks('') ; zticks('')
%     axis square
% 
%     tmp = mean(arrayfun(@(x_) std(s(G.neighbors(x_))),1:size(G.Nodes,1)),'omitnan') ; 
% 
%     text(0.05,0.1,0,num2str(tmp),'Units','normalized')
% end

%%

% mmm = moran_randomization(data_to_reconstruct(~~cortex),MEM,1,'procedure','pair') ;
% 
% ss = nan(length(cortex),1) ;
% ss(cortex) = mmm ; 
% h = quick_trisurf(surface_midthickness,ss)
% h.EdgeColor = 'none'
% view([-45 12 -5])
% waitforbuttonpress

%% get laplacian of surface

filename = sprintf('./gen_data/surflap_modes-%s_%s-%s.mat',num2str(num_modes),surface_interest,hemisphere) ; 

if ~isfile(filename)

    surf_conn = calc_surface_connectivity(surface_midthickness) ; 
    surf_conn = surf_conn(cortex,cortex) ; 
    [~,surf_lap] = calc_LaplacianMatrix(surf_conn) ; 
    [surflap_em,~] = eigs(surf_lap,num_modes,'smallestabs') ; 
    save(filename,'surflap_em')
else
    load(filename)
end

%% look at the null similarities?


loaded_data = load('./NSBLab_repo/data/figures_Nature/Figure1.mat') ;
map_names = fieldnames(loaded_data.task_map_emp) ;


nreps = 500 ; 
surpower = cell(length(map_names),1) ; 
emppower = cell(length(map_names),1) ; 

for idx = 1:length(map_names)

    filename = sprintf('./gen_data/moranres_%s_%s-%s.mat',map_names{idx},surface_interest,hemisphere) ; 

    ll = load(filename,'surr_data_nocort') ; 

    dtr = loaded_data.task_map_emp.(map_names{idx}) ; 
    dtr = dtr(cortex) ;

    % the emp
    [~,tmp_coefs] = calc_eigdecomp_recon_acc(dtr,surflap_em,num_modes) ; 
    [~,emppower{idx}] = calc_power_spectrum(tmp_coefs(:,num_modes)) ;  

    % and the surr

    sur_dat_loop = ll.surr_data_nocort(:,randi(5000,[nreps 1])) ; 
    tmp_betas = nan(200,nreps) ; 

    parfor jdx = 1:nreps

        disp([num2str(idx) ' - ' num2str(jdx)])
  
        [~,sur_coefs] = calc_eigdecomp_recon_acc(sur_dat_loop(:,jdx),surflap_em,num_modes) ; 
    
        [~,tmp_betas(:,jdx)] = calc_power_spectrum(sur_coefs(:,num_modes)) ;  

    end

    surpower{idx} = tmp_betas ; 

end

%% map sim

map_sim = cell(length(map_names),1) ; 

for idx = 1:length(map_names)

    disp(idx)

    filename = sprintf('./gen_data/moranres_%s_%s-%s.mat',map_names{idx},surface_interest,hemisphere) ; 
    ll = load(filename,'surr_data_nocort') ; 

    dtr = loaded_data.task_map_emp.(map_names{idx}) ; 
    dtr = dtr(cortex) ;

    map_sim{idx} = arrayfun(@(x_) corr(ll.surr_data_nocort(:,x_),dtr), 1:size(ll.surr_data_nocort,2) ) ; 

end

%% take a looksie

tiledlayout(2,length(map_names))
set(gcf,'Position', [200 200 1200 600]);
cmap = turbo(length(map_names)+8) ; 
cmap = cmap((end-7):end,:); 

for idx = 1:length(map_names)
    
    nexttile(idx)

    ccc = cmap(idx,:) ; 
    plot(log(surpower{idx}),'Color',[ ccc 0.1 ])
    hold on
    plot(log(emppower{idx}),'Color',[0.5 0.5 0.5 0.5],'LineWidth',0.5)
    hold off

    ylim([-35 0])

    title(map_names_better{idx},'Interpreter','none')

    if idx == 1
        ylabel('Normalized power (log scale) ')
    end

    if idx == 4
        xlabel('Mode')
    end

    nexttile(idx+length(map_names))

    h = histogram(map_sim{idx},25,"Normalization","count") ; 
    h.FaceColor = ccc ; 
    h.EdgeAlpha = 0.1 ; 
    xlim([-0.7 0.7])
    ylim([0 500])

    if idx == 1
        ylabel('Count')
    end

    if idx == 4
        xlabel('Correlation of surrograte to empirical map')
    end
end

%%

outfile = './figures/moranres_spatial_freq.pdf' ; 
orient(gcf,'landscape')
print(gcf,'-dpdf',outfile,'-bestfit','-vector')
close(gcf)

