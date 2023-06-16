clc 
clearvars

%% Load relevant repository MATLAB functions

addpath(genpath('./NSBLab_repo/functions_matlab'));
addpath('./fcn/')

%% Load surface files for visualization

surface_interest = 'fsLR_32k';
hemisphere = 'lh';
mesh_interest = { 'midthickness' 'sphere' } ;

% Load midthickness
[vertices, faces] = read_vtk(sprintf('./NSBLab_repo/data/template_surfaces_volumes/%s_%s-%s.vtk',surface_interest , 'midthickness', hemisphere));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';

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

%% load 10k maps

ll = load('./osf_dl/empirical/fsLR_32k_neurovault_maps_N=10000_lh.mat') ; 
neurovault_maps = ll.neurovault_maps ; 

e_midthick = modes_str.midthickness.eigenmodes(cortex,:) ; 
e_sphere = modes_str.sphere.eigenmodes(cortex,:) ; 

%%

num_modes = 150 ; 
filename = sprintf('./gen_data/predict10k_midthick-sphere_modes%s_%s-%s.mat',...
    num2str(num_modes),surface_interest,hemisphere) ; 

if ~isfile(filename)

    nmaps = size(ll.neurovault_maps,2) ; 
    vault_results = struct() ; 
    recon_mid = nan(nmaps,num_modes) ; 
    recon_sphere = nan(nmaps,num_modes) ; 
    
    parfor idx = 1:nmaps
    
        disp(idx)
    
        data_to_reconstruct = neurovault_maps(:,idx) ; 
    
        recon_mid(idx,:) = calc_eigdecomp_recon_acc(...
            data_to_reconstruct(cortex),...
            e_midthick,num_modes) ; 
    
        recon_sphere(idx,:) = calc_eigdecomp_recon_acc(...
            data_to_reconstruct(cortex),...
            e_sphere,num_modes) ; 
    
    end

    vault_results.recon_mid = recon_mid ;
    vault_results.recon_sphere = recon_sphere ; 

    save(filename,'vault_results') ;

else
    load(filename)
end

%%

tiledlayout(1,3)
 
nexttile
plot_manylines_aspatch(vault_results.recon_mid',[0 0.447 0.741],[0:25])
ylim([-.2 1])
xlim([1 150])
title('Midthickness reconstruction')
ylabel('Reconstruction accuracy')

nexttile
plot_manylines_aspatch(vault_results.recon_sphere',[0.92 0.694 0.125],[0:25])
ylim([-.2 1])
xlim([1 150])
title('Sphere reconstruction')
xlabel('Number of modes')

nexttile
plot_manylines_aspatch(vault_results.recon_mid',[0 0.447 0.741],[25 25],...
    'facealpha', 0.1,'edgecolor',[0 0.447 0.741],'edgealpha',0.5,'LineWidth',1)
hold on
plot_manylines_aspatch(vault_results.recon_sphere',[0.92 0.694 0.125],[25 25],...
    'facealpha', 0.1,'edgecolor',[0.92 0.694 0.125],'edgealpha',0.5,'LineWidth',1)
hold on
plot(mean(vault_results.recon_mid,'omitnan'),...
    'Color',[0 0.447 0.741 0.7],'LineWidth',3)
plot(mean(vault_results.recon_sphere,'omitnan'),...
    'Color',[0.92 0.694 0.125 0.7],'LineWidth',3)
hold off
ylim([-.2 1])
xlim([1 150])
title('Middle 50% and mean reconstruction')

outfile = './figures/predict10k_midthick_sphere.pdf' ; 
orient(gcf,'landscape')
print(gcf,'-dpdf',outfile,'-bestfit','-vector')
