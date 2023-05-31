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

% Load sphere
mesh_interest = 'sphere';
[vertices, faces] = read_vtk(sprintf('./NSBLab_repo/data/template_surfaces_volumes/%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere));
surface_sphere.vertices = vertices';
surface_sphere.faces = faces';

% Load cortex mask
cortex = logical(dlmread(sprintf('./NSBLab_repo/data/template_surfaces_volumes/%s_cortex-%s_mask.txt', surface_interest, hemisphere)));

disp('loaded surfaces')

%% perturb sphere

addpath(genpath('./BrainSpace/matlab'))

G = surface_to_graph(surface_sphere,'mesh') ; 
shortest_paths = distances(G,'Method','unweighted') ; % computes shortest paths
closeness = 1./shortest_paths ; 
closeness(isinf(closeness)) = 0 ;

dist_thr = prctile(shortest_paths,20,'all') ; % arbitrary param, how much 
                                             % local sp neighbourhood to
                                             % capture in graph...
W = closeness .* (shortest_paths < dist_thr) ; 

clear closeness shortest_paths

% center it
W = full(W - mean(W) - mean(W,2) + mean(W,"all",'omitnan')); 

% enforce symmetry
W = triu(W,1) + triu(W)' ; 
assert(issymmetric(W),'not symmetric,problem')

disp('doing decomp')

tic
W = sparse(W) ; 
[MEM,lambda] = eigs(W,10) ; 
lambda = diag(lambda) ;
toc

disp('finished decomp')

%% generate random shapes! 

rng(42)

odir = 'gen_data/rand_surfs' ;
mkdir(odir)

newsurf = struct() ;

for idx = 1:100

    disp(idx)
    
    xyz_ranges = range(surface_sphere.vertices) ; 
    
    % pick three at random
    rand_mem = MEM(:,randi(size(MEM,2),1,3)) ; 
    
    % move the vertices
    newverts = surface_sphere.vertices + (normalize(rand_mem,'range',[-1 1]).*(xyz_ranges.*.15)) ; 
    
    newsurf.coord = newverts' ;
    newsurf.tri = surface_sphere.faces ; 
    
    write_surface(newsurf,sprintf('%s/rand_surf_%03g.gii',odir,idx))

end

%% just a quick viz

quick_trisurf(read_surface(sprintf('%s/rand_surf_%03g.gii',odir,4)))
