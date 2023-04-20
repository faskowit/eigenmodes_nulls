clc 
clearvars

%% Load relevant repository MATLAB functions

addpath(genpath('./NSBLab_repo/functions_matlab'));
addpath('./fcn/')

%% Load surface files for visualization

surface_interest = 'fsLR_32k';
hemisphere = 'lh';
mesh_interest = { 'midthickness' 'sphere' 'veryinflated' 'pial' 'white' } ;

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

%% load some data to reconstruct

loaded_data = load('./NSBLab_repo/data/figures_Nature/Figure1.mat') ;
map_names = fieldnames(loaded_data.task_map_emp) ;
data_to_reconstruct = loaded_data.task_map_emp.(map_names{1}) ; 

%% reconstruct

recon_acc = nan(num_modes,length(mesh_interest)) ;

for idx = 1:length(mesh_interest)

    disp(idx)
    recon_acc(:,idx) = calc_eigdecomp_recon_acc(data_to_reconstruct(cortex),...
        modes_str.(mesh_interest{idx}).eigenmodes(cortex,:), ...
        num_modes) ; 

end

%% plot it

plot(recon_acc,'LineWidth',2)
legend(mesh_interest,'Location','southeast')
xlabel('number of modes')
ylabel('reconstruction accuracy')
