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

for map_idx = 1:length(map_names)

    filename = sprintf('./gen_data/altsurfs_%s_%s-%s.mat',map_names{map_idx},surface_interest,hemisphere) ;

    if ~isfile(filename)

        data_to_reconstruct = loaded_data.task_map_emp.(map_names{map_idx}) ;
    
        %% reconstruct
        
        alt_recon_acc = nan(num_modes,length(mesh_interest)) ;
        
        parfor idx = 1:length(mesh_interest)
        
            disp(idx)
            alt_recon_acc(:,idx) = calc_eigdecomp_recon_acc(data_to_reconstruct(cortex),...
                modes_str.(mesh_interest{idx}).eigenmodes(cortex,:), ...
                num_modes) ; 
        
        end
    
        save(filename,"alt_recon_acc") ;
    
    else
        disp('already generated')
    end

end

%% plot it

loaded_data = load('./NSBLab_repo/data/figures_Nature/Figure1.mat') ;
map_names = fieldnames(loaded_data.task_map_emp) ;

tiledlayout(1,length(map_names))
set(gcf,'Position', [200 200 2000 400]);

map_names_better = { 'Social' 'Motor' 'Gambling' 'WM' ...
    'Language' 'Emotion' 'Relational' } ;

for idx = 1:length(map_names)

    filename = sprintf('./gen_data/altsurfs_%s_%s-%s.mat',map_names{idx},surface_interest,hemisphere) ;

    ll = load(filename) ; 

    recon_acc = ll.alt_recon_acc ;

    nexttile(idx)
    hold on 
    pp = plot(recon_acc,'LineWidth',1) ;
    
    for jdx = 1:length(pp)
        pp(jdx).Color = [ pp(jdx).Color 0.5 ] ;
    end

    hold off
    
    for pdx = 1:length(pp)
        pp(pdx).Color(4) = 0.8 ;
    end

    xlim([1 num_modes])
    ylim([-0.25 1])

    h = gca ;
    h.XTick = [ 1 h.XTick] ; 

    if idx == 1
        ylabel('Reconstruction accuracy')
    end

    if idx == 4
        xlabel('Geometric eigenmodes used for reconstruction')
    end

    if idx == 7
        legend(mesh_interest,'Location','southeastoutside')
    end

    title(map_names_better{idx},'Interpreter','none')
    

end

%% save it

outfile = './figures/alt_eigenmodes.pdf' ; 
orient(gcf,'landscape')
print(gcf,'-dpdf',outfile,'-bestfit','-vector')

%% and some viz of the other figures

% some texture
% [u,v] = pca(vertices,'NumComponents',2) ; 

addpath('fcn/external/')

loopover = { 'sphere' 'veryinflated' 'pial' 'midthickness' 'white' } ; 

% compute curvature stuff
curvStruct = struct() ; 

for idx = 1:length(loopover)
        
    disp(idx)

    sss = loopover{idx} ;

    % Load midthickness
    [vertices, faces] = read_vtk(sprintf('./NSBLab_repo/data/template_surfaces_volumes/%s_%s-%s.vtk',surface_interest , sss, hemisphere));
    surface_view.vertices = vertices';
    surface_view.faces = faces';

    % u = compute_normal(surface_midthickness.vertices,surface_midthickness.faces) ; 
    opts.curvature_smoothing = 10 ; 
        opts.verb = 0 ; 

    [~,~,curvStruct(idx).cmin,curvStruct(idx).cmax,...
        ~,curvStruct(idx).cgaus,curvStruct(idx).normal] = ...
        compute_curvature(surface_view.vertices,...
        surface_view.faces, opts) ; 

end

%%

big_mins = [ curvStruct.cmin ] ;
big_maxs = [ curvStruct.cmax ] ;

tmp = abs(big_mins(:,2:5)) + abs(big_maxs(:,2:5)) ; 
clim_min = mean(prctile(tmp,1));
clim_max = mean(prctile(tmp,99));

tiledlayout(1,5)

for idx = 1:length(loopover)

    nexttile
        
    sss = loopover{idx} ;

    % Load midthickness
    [vertices, faces] = read_vtk(sprintf('./NSBLab_repo/data/template_surfaces_volumes/%s_%s-%s.vtk',surface_interest , sss, hemisphere));
    surface_view.vertices = vertices';
    surface_view.faces = faces';

%     h = quick_trisurf(surface_midthickness,cgaus) ; 
    c1 = curvStruct(idx).cmin ; 
    c2 = curvStruct(idx).cmax ;
    h = quick_trisurf(surface_view,abs(c1)+abs(c2)) ; 
    h.EdgeColor = "none";
%     view([-100 0 10]);
    view([-90 0 0]);
    clim([clim_min clim_max])
    material shiny
    % camlight right
    lighting gouraud
    xticks('') ; yticks('') ; zticks('')
    axis equal
    axis off

    title(sss)

end

outfile = './figures/alt_shapes.pdf' ; 
orient(gcf,'landscape')
print(gcf,'-dpdf',outfile,'-bestfit','-vector')

%% visualize the curvatures... but on the sphere

tiledlayout(4,5)

% Load sphere
sss = 'sphere' ;
[vertices, faces] = read_vtk(sprintf('./NSBLab_repo/data/template_surfaces_volumes/%s_%s-%s.vtk',surface_interest , sss, hemisphere));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';

for idx = 1:length(loopover)

    nexttile
        
    sss = loopover{idx} ;

%     h = quick_trisurf(surface_midthickness,cgaus) ; 
    c1 = curvStruct(idx).cmin ; 
    c2 = curvStruct(idx).cmax ;
    h = quick_trisurf(surface_midthickness,abs(c1)+abs(c2)) ; 
    h.EdgeColor = "none";
%     view([-100 0 10]);
    view([-90 0 0]);
    clim([clim_min clim_max])
    material shiny
    % camlight right
    lighting gouraud
    xticks('') ; yticks('') ; zticks('')
%     axis equal
%     axis off

    title(sss)

end

aaa = { 'x' 'y' 'z' } ; 

for jdx = 1:3
% and now to the normals
for idx = 1:length(loopover)

    nexttile
        
    sss = loopover{idx} ;

    h = quick_trisurf(surface_midthickness,curvStruct(idx).normal(jdx,:)) ; 
    h.EdgeColor = "none";
%     view([-100 0 10]);
    view([-90 0 0]);
    clim([-1 1])
    material shiny
    % camlight right
    lighting gouraud
    xticks('') ; yticks('') ; zticks('')
%     axis equal
%     axis off

    %title(sss)
    if idx == 1
        zlabel([ 'normals ' aaa{jdx} ])
    end

end
end 

%%

outfile = './figures/shape_view_on_spheres.png' ; 
orient(gcf,'landscape')
% print(gcf,'-dpdf',outfile,'-bestfit','-vector')
print(gcf,'-dpng',outfile)

%% plot inds 

% Load midthickness
[vertices, faces] = read_vtk(sprintf('./NSBLab_repo/data/template_surfaces_volumes/%s_%s-%s.vtk',surface_interest , 'sphere', hemisphere));
surface_view_sphere.vertices = vertices';
surface_view_sphere.faces = faces';

for jdx = 1:3
        
    tiledlayout(2,5)

    disp(jdx)

for idx = 1:length(loopover)

    disp(idx)

    nexttile(idx)
        
    sss = loopover{idx} ;

    % Load midthickness
    [vertices, faces] = read_vtk(sprintf('./NSBLab_repo/data/template_surfaces_volumes/%s_%s-%s.vtk',surface_interest , sss, hemisphere));
    surface_view.vertices = vertices';
    surface_view.faces = faces';

    h = quick_trisurf(surface_view,surface_view.vertices(:,jdx)) ; 
    h.EdgeColor = "none";
%     view([-100 0 10]);
    view([-90 0 0]);
    material shiny
    % camlight right
    lighting gouraud
    xticks('') ; yticks('') ; zticks('')
    
    nexttile(idx+5)

    h = quick_trisurf(surface_view_sphere,surface_view.vertices(:,jdx)) ; 
    h.EdgeColor = "none";
%     view([-100 0 10]);
    view([-90 0 0]);
    material shiny
    % camlight right
    lighting gouraud
    xticks('') ; yticks('') ; zticks('')
end

set(gcf,'Position', [200 200 1200 800]);

outfile = ['./figures/shape_inds_spheres_' aaa{jdx} '.png'] ; 
orient(gcf,'landscape')
% print(gcf,'-dpdf',outfile,'-bestfit','-vector')
print(gcf,'-dpng',outfile)

close(gcf)

end
