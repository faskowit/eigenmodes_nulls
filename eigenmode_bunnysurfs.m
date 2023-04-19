
clc
clearvars

%% add paths

addpath(genpath('./NSBLab_repo/functions_matlab'));
addpath(genpath('./fcn/'))

%% read fsLR

surface_interest = 'fsLR_32k';
hemisphere = 'lh';
mesh_interest = 'midthickness';

% Load midthickness
[vertices, faces] = read_vtk(sprintf('./NSBLab_repo/data/template_surfaces_volumes/%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';

% and cortex
cortex = logical(dlmread(sprintf('./NSBLab_repo/data/template_surfaces_volumes/%s_cortex-%s_mask.txt', surface_interest, hemisphere)));

% what fraction of the midthick surface is the "hole", so we can approx.
medialwall_frac = sum(cortex==0)/length(cortex) ;

%% read bunny

[F,V] = stanford_bunny() ;

%% Display bunny

clf;
trisurf(F,V(:,1),V(:,2),V(:,3),flipud(1:length(F)))

%% write the bunny surface as a gifti

% make a brainspace compatible surface
brainspace_bunny = struct() ;
brainspace_bunny.tri = F ;
brainspace_bunny.coord = V' ;

write_surface(brainspace_bunny,'./gen_data/bunny.gii')

%% calculate places a 'fake medial wall' might go

addpath(genpath('./BrainSpace/matlab'))

targ_medialwall_sz = floor(length(V)*medialwall_frac) ; 

bunny_g = surface_to_graph(brainspace_bunny,'mesh') ;
bunny_nei = get_nei_mat(bunny_g) ; 

rng(42)
rand_seeds = randi(length(brainspace_bunny.coord),100,1) ; 

odir = './gen_data/bunny_medial/' ; 
mkdir(odir)

% make a bunch of fake medial walls
for idx = 1:length(rand_seeds)

    disp(idx)

    rs = rand_seeds(idx) ; 

    fake_med_wall = ones(length(V),1) ; 
    % seed it
    fake_med_wall(rs) = 2 ;
    
    % dilate once
    fake_med_wall = dil_surf_parc(fake_med_wall,bunny_nei,1) ; 
    % now dilate until exceed target size
    while sum(fake_med_wall==2) < targ_medialwall_sz
        fake_med_wall = dil_surf_parc(fake_med_wall,bunny_nei,1,1) ;
    end
    
    % and now a lil trick to make the borders a lil more smooth
    fake_med_wall = dil_surf_parc(fake_med_wall,bunny_nei,2,1) ; % erode 
    fake_med_wall = dil_surf_parc(fake_med_wall,bunny_nei,1,1) ; % dilate
    fake_med_wall = dil_surf_parc(fake_med_wall,bunny_nei,2,1) ; % erode

    writematrix(fake_med_wall-1,sprintf('%s/bunny_med_wall_%03g.txt',odir,idx))

end

