function [ recon_acc, recon_beta] = calc_eigdecomp_recon_acc(data_to_recon,eigenmodes,num_modes,corrtype)
% recon_acc: vector of correlation coeffs btwn data and recon'ed data
%
% Josh Faskowitz

if nargin < 4
    corrtype = 'p' ; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit the betas

% recon_beta = zeros(num_modes, num_modes);
% for mode = 1:num_modes
%     basis_tmp = eigenmodes(:, 1:mode);
%     recon_beta(1:mode,mode) = calc_eigendecomposition(data_to_recon, basis_tmp, 'matrix');
% end

% arrayfun it
recon_beta_cell = arrayfun(@(ind_) ...
    calc_eigendecomposition(data_to_recon, ...
        eigenmodes(:, 1:ind_), 'matrix' ...
    ), 1:num_modes, 'UniformOutput',false) ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assess correlation btwn data and recon

% recon_acc = zeros(1, num_modes);               
% for mode = 1:num_modes
%     recon_temp = eigenmodes(:, 1:mode)*recon_beta(1:mode,mode);
%     recon_acc(mode) = corr(data_to_recon, recon_temp);
% end

% arrayfun it
recon_acc = arrayfun(@(ind_) ...
    corr(data_to_recon, ...
         eigenmodes(:, 1:ind_)*recon_beta_cell{ind_}, ...
         'type',corrtype ...
    ), 1:num_modes) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% provide recon betas if wanted

if nargout > 1
    % make the recon beta matrix if needed
    recon_beta = zeros(num_modes, num_modes) ; 
    for idx = 1:num_modes
        recon_beta(1:idx,idx) = recon_beta_cell{idx} ;
    end
end
