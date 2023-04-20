
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

%% now try to upsample the bunny

[VV,FF] = upsample(V,F,'Iterations',3) ; 

% code copied from:
% https://www.numerical-tours.com/matlab/meshwav_3_simplification/

faces1 = FF';
vertex1 = VV';
vertex0 = VV' ;
while sum(~all(isinf(vertex1))) > length(cortex)
    edges = compute_edges(faces1);
    D = vertex0(:,edges(1,:)) - vertex0(:,edges(2,:));
    D = sum(D.^2,1);
    [~,k] = min(D);
    e = edges(:,k);
    vertex1(:,e(1)) = mean( vertex1(:,e),2 );
    vertex1(:,e(2)) = Inf;
    faces1(faces1==e(2)) = e(1);
    a = sum( diff(sort(faces1))==0 );
    faces1(:,a>0) = [];
    sum(~all(isinf(vertex1)))
end

F = faces1' ; 
V = vertex1' ; 

%% rectify the inf rows and faces

%uniquify the vertices
[newverts,i,j] = unique(V, 'rows');
newverts = newverts(1:end-1,:) ; % get rid of the infs
newfaces = [ j(F(:,1)) j(F(:,2)) j(F(:,3)) ] ; 

trisurf(newfaces,newverts(:,1),newverts(:,2),newverts(:,3))

%% write the bunny surface as a gifti

addpath(genpath('./BrainSpace/matlab'))

% make a brainspace compatible surface
brainspace_bunny = struct() ;
brainspace_bunny.tri = newfaces ;
brainspace_bunny.coord = newverts' ;

write_surface(brainspace_bunny,'./gen_data/bunny.gii')
% after you write the surface, you need to use freesurfer to conver to vtk
% using mris_convert

%% calculate places a 'fake medial wall' might go

targ_medialwall_sz = floor(length(cortex)*medialwall_frac) ; 

bunny_g = surface_to_graph(brainspace_bunny,'mesh') ;
bunny_nei = get_nei_mat(bunny_g) ; 

rng(42)
rand_seeds = randi(length(cortex),100,1) ; 

odir = './gen_data/bunny_medial/' ; 
mkdir(odir)

% make a bunch of fake medial walls
for idx = 1:length(rand_seeds)

    disp(idx)

    rs = rand_seeds(idx) ; 

    fake_med_wall = ones(length(cortex),1) ; 
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

    writematrix(~(fake_med_wall-1),sprintf('%s/bunny_med_wall_%03g.txt',odir,idx))

end

%% viz some 

trisurf(newfaces,newverts(:,1),newverts(:,2),newverts(:,3),...
    load(sprintf('%s/bunny_med_wall_%03g.txt',odir,51)))


