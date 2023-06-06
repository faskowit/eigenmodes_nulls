clc 
clearvars

%% Load relevant repository MATLAB functions

addpath(genpath('./NSBLab_repo/functions_matlab'));
addpath('./fcn/')

%% Load surface files for visualization

surface_interest = 'fsLR_32k';
hemisphere = 'lh';

% Load midthickness
[vertices, faces] = read_vtk(sprintf('./NSBLab_repo/data/template_surfaces_volumes/%s_%s-%s.vtk',surface_interest , 'midthickness', hemisphere));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';

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

%% load some data to reconstruct

loaded_data = load('./NSBLab_repo/data/figures_Nature/Figure1.mat') ;
map_names = fieldnames(loaded_data.task_map_emp) ;

for map_idx = 1:length(map_names)

    blob_results = struct() ; 
    blob_results.name = map_names{map_idx} ; 

    %% Get the empirical accuracy  
    data_to_reconstruct = loaded_data.task_map_emp.(map_names{map_idx}) ; 
   
    %% do the regular recon
    
    blob_results.surfrecon = calc_eigdecomp_recon_acc(data_to_reconstruct(cortex),...
        eigenmodes(cortex,:), ...
        num_modes) ; 

    %% and do the iterations with the blobs
    
    filename = sprintf('./gen_data/randsurfs_%s_%s-%s.mat',map_names{map_idx},surface_interest,hemisphere) ; 
    
    if ~isfile(filename)
    
        blobrecons = nan(200,100) ; 
        
        parfor idx = 1:100
            disp(idx)
        
            loadsurf = load(sprintf('gen_data/randsurfs_eigen/rand_surf_emode_200_%03g.txt',idx)) ; 
            
            blobrecons(:,idx) = calc_eigdecomp_recon_acc(data_to_reconstruct(cortex),...
                loadsurf(cortex,:), ...
                num_modes) ; 
        end
    
        blob_results.blobrecons = blobrecons ;

        save(filename,'blob_results')
    else
        load(filename)
    end

end

%%

tiledlayout(2,length(map_names))
set(gcf,'Position', [200 200 1200 600]);

map_names_better = { 'Social' 'Motor' 'Gambling' 'Working memory' ...
    'Language' 'Emotion' 'Relational' } ;

surface_interest = 'fsLR_32k';
hemisphere = 'lh';
num_modes = 200 ; 

for idx = 1:length(map_names)

    filename = sprintf('./gen_data/randsurfs_%s_%s-%s.mat',map_names{idx},surface_interest,hemisphere) ;

    ll = load(filename) ; 

    perm_acc = ll.blob_results.blobrecons ;
    recon_acc = ll.blob_results.surfrecon ;

    nexttile(idx)
    plot_manylines(perm_acc,'Color',[0 0.4470 0.7410 0.05],'LineWidth',2) 
    hold on 
    plot(recon_acc,'r','LineWidth',2)
    hold off
    
    xlim([1 num_modes])
    ylim([-0.25 1])

    if idx == 1
        ylabel('recon accuracy')
    end

    if idx == 4
        xlabel('eigenmodes used for recon')
    end

    title(map_names_better{idx},'Interpreter','none')
    
    nexttile(idx+length(map_names))
    pvals = ( sum(bsxfun(@gt,perm_acc,recon_acc'),2) + 1) ./ (size(ll.blob_results.blobrecons,2)+1) ; 
    plot(pvals,'mo-','LineWidth',2,'MarkerSize',3)
    
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

% figure
% plot_manylines(blobrecons)
% hold on 
% plot(surfrecon)
% xlim([1 num_modes])

%% save it

outfile = './figures/randsurfs_eigenmodes.pdf' ; 
orient(gcf,'landscape')
print(gcf,'-dpdf',outfile,'-bestfit','-vector')

%% viz random shapes?

set(gcf,'Position', [200 200 1200 1200]);

tiledlayout(4,4)

for idx = 1:16

    nexttile()

    [randsurf.vertices, randsurf.faces] = read_vtk(sprintf('./gen_data/rand_surfs/rand_surf_%03g.vtk',idx));
    
    h = quick_trisurf(randsurf) ;
    h.EdgeColor = "none";
    material shiny
    % camlight headlight
    lighting gouraud
    xticks('') ; yticks('') ; zticks('')

end

outfile = './figures/randsurfs_view_shapes.pdf' ; 
orient(gcf,'landscape')
print(gcf,'-dpdf',outfile,'-vector','-bestfit')
