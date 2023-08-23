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

clear vertices faces
disp('loaded surfaces')

%% Setup for BrainSpace

addpath(genpath('./BrainSpace/matlab'))

% make a matlab surface that BrainSpace is happy with
brainspace_midthick = struct() ; 
brainspace_midthick.tri = surface_midthickness.faces ;  
brainspace_midthick.coord = surface_midthickness.vertices' ; 

%% load moran mem

hemisphere = 'lh';
num_modes = 200;

tic
filename = sprintf('./gen_data/moranMEM5000_BS_modes-%s_%s-%s.mat',num2str(num_modes),surface_interest,hemisphere) ; 
if ~isfile(filename)
    
    filenameMEM = sprintf('./gen_data/moranMEM_BS_modes-%s_%s-%s.mat',num2str(num_modes),surface_interest,hemisphere) ; 
    if ~isfile(filenameMEM)
        disp('computing MEM')
        datestr(now,'HH:MM:SS')
        MEM = compute_mem(brainspace_midthick,...
            'mask',~cortex,'n_ring',1e4,...
            'distance_metric','geodesic') ; 
        save(filename,"MEM")
    else
        load(filenameMEM)
    end
    disp('saving MEM 5000')
    MEMredu = MEM(:,1:5000) ; 
    save(filenameMEM,'MEMredu')
    clear MEM
else
    disp('loading MEM 5000')
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

%%

% Loop over zstat maps and do variogram tests

% nperms = 1e2 ;
nperms = 5e3; 
% nperms  = 100 ; 

loaded_data = load('./NSBLab_repo/data/figures_Nature/Figure1.mat') ;
map_names = fieldnames(loaded_data.task_map_emp) ;

eigenmodes_nocort = eigenmodes(cortex,:) ;

for map_idx = 1:length(map_names)

    signflip_results = struct() ; 
    signflip_results.name = map_names{map_idx} ; 

    %% Get the empirical accuracy  
    data_to_reconstruct = loaded_data.task_map_emp.(map_names{map_idx}) ; 
    
    % viz it?
    % draw_surface_bluewhitered_dull(surface_midthickness, data_to_reconstruct, hemisphere, find(cortex==0), 1);
    
    %% do spins and recort results
    
    filename = sprintf('./gen_data/signflipNamp_%s_%s-%s.mat',map_names{map_idx},surface_interest,hemisphere) ; 
    
    if ~isfile(filename)

        % compute original accuracy
        signflip_results.recon_acc = calc_eigdecomp_recon_acc(data_to_reconstruct(cortex),eigenmodes(cortex,:),num_modes) ; 
        
        % make the randomized data (takes a lil bit of time)
        disp('making signflip dat')
        surr_data_nocort = basis_signflipNamp_surr(...
            data_to_reconstruct(cortex),MEMredu(:,1:5000),nperms) ; % block 100

        perm_acc = nan(num_modes,nperms) ; 

        % measure reconstruct accuracy on spin data
        parfor idx = 1:nperms
        
        %     if mod(idx,10) == 1
        %         disp([ num2str(idx) ' of ' num2str(nperms) ])
        %     end
        
            disp([ num2str(map_idx) ' -     ' num2str(idx) ] )
        
            perm_acc(:,idx) = calc_eigdecomp_recon_acc(surr_data_nocort(:,idx),...
                                               eigenmodes_nocort,num_modes) ; 
        
        end
    
        % record results in the struct 
        signflip_results.perm_acc = perm_acc ; 

        % save it
        save(filename,'signflip_results','surr_data_nocort')

    end

end

%% vizz it

surface_interest = 'fsLR_32k';
hemisphere = 'lh';
num_modes = 200 ; 

loaded_data = load('./NSBLab_repo/data/figures_Nature/Figure1.mat') ;
map_names = fieldnames(loaded_data.task_map_emp) ;

tiledlayout(2,length(map_names))
set(gcf,'Position', [200 200 1200 600]);

% nperms = 5e3; 

map_names_better = { 'Social' 'Motor' 'Gambling' 'WM' ...
    'Language' 'Emotion' 'Relational' } ;

for idx = 1:length(map_names)

    filename = sprintf('./gen_data/signflipNamp_%s_%s-%s.mat',map_names{idx},surface_interest,hemisphere) ; 

    ll = load(filename,'signflip_results') ; 

    perm_acc = ll.signflip_results.perm_acc ;
    recon_acc = ll.signflip_results.recon_acc ;

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

outfile = './figures/signflipNamp_eigenmodes_patchversion.pdf' ; 
orient(gcf,'landscape')
print(gcf,'-dpdf',outfile,'-bestfit','-vector')
close(gcf)

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

nreps = 500 ; 
surpower = cell(length(map_names),1) ; 
emppower = cell(length(map_names),1) ; 

for idx = 1:length(map_names)

    filename = sprintf('./gen_data/signflipNamp_%s_%s-%s.mat',map_names{idx},surface_interest,hemisphere) ; 

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

    filename = sprintf('./gen_data/signflipNamp_%s_%s-%s.mat',map_names{idx},surface_interest,hemisphere) ; 
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
    % plot(log(surpower{idx}),'Color',[ ccc 0.1 ])
    % hold on
    % plot(log(emppower{idx}),'Color',[0.5 0.5 0.5 0.5],'LineWidth',0.5)
    % hold off
    % 
    % ylim([-35 0])
    % 
    % title(map_names_better{idx},'Interpreter','none')
    % 
    % if idx == 1
    %     ylabel('Normalized power (log scale) ')
    % end

    if idx == 4
        xlabel('Mode')
    end

    nexttile(idx+length(map_names))

    h = histogram(map_sim{idx},25,"Normalization","count") ; 
    h.FaceColor = ccc ; 
    h.EdgeAlpha = 0.1 ; 
    xlim([-0.7 0.7])
    ylim([0 600])

    if idx == 1
        ylabel('Count')
    end

    if idx == 4
        xlabel('Correlation of surrograte to empirical map')
    end
end

%%

outfile = './figures/signflipNamp_spatial_freq.pdf' ; 
orient(gcf,'landscape')
print(gcf,'-dpdf',outfile,'-bestfit','-vector')
close(gcf)

%% look at some example surrogate data


idx=1;
filename = sprintf('./gen_data/signflipNamp_%s_%s-%s.mat',map_names{idx},surface_interest,hemisphere) ; 
ll = load(filename,'surr_data_nocort') ; 
dtr = loaded_data.task_map_emp.(map_names{idx}) ; 

%%

tiledlayout(4,4)
set(gcf,'Position', [200 200 800 800]);

nexttile
quick_trisurf(surface_midthickness,dtr)
view(-92.7,7.2)
material shiny
% camlight right
lighting gouraud
xticks('') ; yticks('') ; zticks('')
grid off
axis equal
axis off
title('empirical')

rng(42)
rr = randi(5000,[1 15]);
for  idx = rr
    
    ss = ll.surr_data_nocort(:,idx) ;
    dd = cortex .*1 ; 
    dd(cortex) = ss ;

    nexttile
    quick_trisurf(surface_midthickness,dd)
    view(-92.7,7.2)
    material shiny
    % camlight right
    lighting gouraud
    xticks('') ; yticks('') ; zticks('')
    grid off
    axis equal
    axis off

end

%%

outfile = './figures/signflipNamp_example_surrs.png' ; 
orient(gcf,'landscape')
print(gcf,'-dpng',outfile)
outfile = './figures/signflipNamp_example_surrs.pdf' ; 
print(gcf,'-dpdf',outfile,'-bestfit','-vector')
close(gcf)

