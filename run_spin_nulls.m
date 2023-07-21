clc 
clearvars

%% Load relevant repository MATLAB functions

addpath(genpath('./NSBLab_repo/functions_matlab'));
addpath('./fcn/')

%% setup some variables

nperms = 5e3; 
surface_interest = 'fsLR_32k';
hemisphere = 'lh';
mesh_interest = 'midthickness';
num_modes = 200;

map_names_better = { 'Social' 'Motor' 'Gambling' 'WM' ...
    'Language' 'Emotion' 'Relational' } ;

%% Load surface files for visualization

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
        spin_inds_smaller = spin_inds(:,randperm(size(spin_inds,2),nperms)) ; 
        
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VIZ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% draw some spun data

filename = sprintf('./gen_data/spininds_%s-%s.mat',surface_interest,hemisphere) ; 
spindat = load(filename) ; 

loaded_data = load('./NSBLab_repo/data/figures_Nature/Figure1.mat') ;
map_names = fieldnames(loaded_data.task_map_emp) ;
dtr = loaded_data.task_map_emp.(map_names{1}) ; 
dtr(~cortex) = min(dtr) * 1.1 ; 

%% view

% Load sphere
mesh_interest = 'sphere';
[vertices, faces] = read_vtk(sprintf('./NSBLab_repo/data/template_surfaces_volumes/%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere));
surface_sphere.vertices = vertices';
surface_sphere.faces = faces';

rng(42)
spind = spindat.spin_inds(:,randperm(size(spindat.spin_inds,2),16)) ; 

set(gcf,'Position', [200 200 1200 1200]);
tiledlayout(4,4)
for idx = 1:16

    nexttile
    h = quick_trisurf(surface_sphere,dtr(spind(:,idx))) ; 
    h.EdgeColor = "none";
    material shiny
    % camlight right
    lighting gouraud
    xticks('') ; yticks('') ; zticks('')

end

outfile = './figures/spin_data_example.pdf' ; 
print(gcf,'-dpdf',outfile,'-vector','-bestfit')

%% other supporting pictures

% normal data on midthickness
fig = figure('Position', [200 200 600 300]);
h = quick_trisurf(surface_midthickness,dtr)
h.EdgeColor = "none";
view([-100 0 10]);
material shiny
% camlight right
lighting gouraud
xticks('') ; yticks('') ; zticks('')
axis equal
axis fill
axis off

outfile = './figures/unspun_midthickness.pdf' ; 
print(gcf,'-dpdf',outfile,'-vector','-bestfit')


% spun data on midthickness
fig = figure('Position', [200 200 600 300]);
h = quick_trisurf(surface_midthickness,dtr(spindat.spin_inds(:,4242)))
h.EdgeColor = "none";
view([-100 0 10]);
material shiny
% camlight right
lighting gouraud
xticks('') ; yticks('') ; zticks('')
axis equal
axis fill
axis off

outfile = './figures/spun_midthickness.pdf' ; 
print(gcf,'-dpdf',outfile,'-vector','-bestfit')

% un-spun data on sphere
h = quick_trisurf(surface_sphere,dtr) ; 
h.EdgeColor = "none";
material shiny
% camlight right
axis equal
axis fill
lighting gouraud
xticks('') ; yticks('') ; zticks('')

outfile = './figures/unspun_sphere.pdf' ; 
print(gcf,'-dpdf',outfile,'-vector','-bestfit')

% spun data on sphere
h = quick_trisurf(surface_sphere,dtr(spindat.spin_inds(:,4242))) ; 
h.EdgeColor = "none";
material shiny
% camlight right
lighting gouraud
axis equal
axis fill
xticks('') ; yticks('') ; zticks('')

outfile = './figures/spun_sphere.pdf' ; 
print(gcf,'-dpdf',outfile,'-vector','-bestfit')

%% viz it as PATCH! PREFERRED METHOD

loaded_data = load('./NSBLab_repo/data/figures_Nature/Figure1.mat') ;
map_names = fieldnames(loaded_data.task_map_emp) ;

tiledlayout(2,length(map_names))
set(gcf,'Position', [200 200 1200 600]);

for idx = 1:length(map_names)

    filename = sprintf('./gen_data/spinres_%s_%s-%s.mat',map_names{idx},surface_interest,hemisphere) ; 

    ll = load(filename) ; 

    perm_acc = ll.spin_results.perm_acc ;
    recon_acc = ll.spin_results.recon_acc ;

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
    pvals = ( sum(bsxfun(@gt,perm_acc,recon_acc'),2) + 1) ./ (nperms+1) ; 
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

outfile = './figures/spinnulls_eigenmodes_patchversion.pdf' ; 
orient(gcf,'landscape')
print(gcf,'-dpdf',outfile,'-bestfit','-vector')

%% visualize the colormap
figure
imagesc(ppspan) ; cb = colorbar ; colormap(cmap)
cb.Label.String = 'Data range' ; 

outfile = './figures/spinnulls_eigenmodes_patchversion_cb.pdf' ; 
orient(gcf,'landscape')
print(gcf,'-dpdf',outfile,'-bestfit','-vector')

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

%% get some power based on laplacian basis
% a lil extra analysis followed by more viz

nreps = 500 ; 
surpower = cell(length(map_names_better),1) ; 
emppower = cell(length(map_names_better),1) ; 

loaded_data = load('./NSBLab_repo/data/figures_Nature/Figure1.mat') ;
map_names = fieldnames(loaded_data.task_map_emp) ;

for idx = 1:length(map_names_better)


    dtr = loaded_data.task_map_emp.(map_names{idx}) ; 
    %dtr = dtr(cortex) ;

    % the emp
    [~,tmp_coefs] = calc_eigdecomp_recon_acc(dtr(cortex),surflap_em,num_modes) ; 
    [~,emppower{idx}] = calc_power_spectrum(tmp_coefs(:,num_modes)) ;  

    % and the surr

    tmp_betas = nan(200,nreps) ; 
    % get nreps rotated data
    spin_inds_smaller = spin_inds(:,randperm(size(spin_inds,2),nreps)) ; 

    sur_dat_loop = cell2mat(...
        arrayfun(@(x_) dtr(spin_inds_smaller(:,x_)), 1:nreps, ...
        'UniformOutput', false)) ;

    parfor jdx = 1:nreps

        disp([num2str(idx) ' - ' num2str(jdx)])
  
        tmp_sur = sur_dat_loop(:,jdx) ;

        perm_mask = (cortex & ~isnan(tmp_sur)) ; 

        [~,sur_coefs] = calc_eigdecomp_recon_acc(...
            tmp_sur(perm_mask),surflap_em(perm_mask(cortex),:),num_modes) ; 
    
        [~,tmp_betas(:,jdx)] = calc_power_spectrum(sur_coefs(:,num_modes)) ;  

    end

    surpower{idx} = tmp_betas ; 

end

%% map sim

map_sim = cell(length(map_names),1) ; 

for idx = 1:length(map_names)

    disp(idx)

    dtr = loaded_data.task_map_emp.(map_names{idx}) ; 
    % dtr = dtr(cortex) ;

    spin_inds_smaller = spin_inds(:,randperm(size(spin_inds,2),nperms)) ; 
    sur_dat_loop = cell2mat(...
        arrayfun(@(x_) dtr(spin_inds_smaller(:,x_)), 1:nperms, ...
        'UniformOutput', false)) ;


    map_sim{idx} = arrayfun(@(x_) ...
        corr(sur_dat_loop(:,x_),dtr,'rows','complete'), 1:nperms ) ; 

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

    ylim([-30 0])

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
    ylim([0 1000])

    if idx == 1
        ylabel('Count')
    end

    if idx == 4
        xlabel('Correlation of surrograte to empirical map')
    end
end

%%

outfile = './figures/spinnulls_spatial_freq.pdf' ; 
orient(gcf,'landscape')
print(gcf,'-dpdf',outfile,'-bestfit','-vector')
close(gcf)
