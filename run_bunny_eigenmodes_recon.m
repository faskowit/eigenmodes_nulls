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

% load bunny
[vertices, faces] = read_vtk('gen_data/bunny.vtk');
surface_bunny.vertices = vertices';
surface_bunny.faces = faces';

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

%% load bunnymodes

bunnymodes = load('gen_data/bunny_emode_200.txt') ; 

% viz for fun
trisurf(surface_bunny.faces,surface_bunny.vertices(:,1),...
    surface_bunny.vertices(:,2),surface_bunny.vertices(:,3),...
    bunnymodes(:,5))

%% load some data to reconstruct

loaded_data = load('./NSBLab_repo/data/figures_Nature/Figure1.mat') ;
map_names = fieldnames(loaded_data.task_map_emp) ;
data_to_reconstruct = loaded_data.task_map_emp.(map_names{1}) ; 

%%

surfrecon = calc_eigdecomp_recon_acc(data_to_reconstruct(cortex),...
        eigenmodes(cortex,:), ...
        num_modes) ; 

bunnyrecon = calc_eigdecomp_recon_acc(data_to_reconstruct(cortex),...
        bunnymodes(cortex,:), ...
        num_modes) ; 

%% and do the iterations with the 'fake medial wall'

filename = 'gen_data/bunny_modes.mat' ; 

if ~isfile(filename)

    bunnyrecons = nan(200,100) ; 
    
    for idx = 1:100
        disp(idx)
    
        loadbunny = load(sprintf('gen_data/bunny_medial_eigen/bunny_emode_200_%03g.txt',idx)) ; 
        loadwall = load(sprintf('gen_data/bunny_medial/bunny_med_wall_%03g.txt',idx)) ; 
    
        intersect_wall = loadwall & cortex ; 
    
        bunnyrecons(:,idx) = calc_eigdecomp_recon_acc(data_to_reconstruct(intersect_wall),...
            loadbunny(intersect_wall,:), ...
            num_modes) ; 
    end

    save(filename,'bunnyrecons')
else
    load(filename)
end

%%

plot(bunnyrecons)

