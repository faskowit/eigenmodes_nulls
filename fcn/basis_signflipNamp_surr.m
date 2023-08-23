function [ surr_data , yrecon, surr_coefs, coefs ] = basis_signflipNamp_surr(...
    y,basisset,nperms,randwei,num2flip) 

dimtouse = size(basisset,2) ; 

if nargin < 4
    % weights when picking bases to flip
    randwei = 'uniform' ; 
end

if nargin < 5
    % how many bases to flip
    num2flip = floor(dimtouse/2) ;
end

if length(y) ~= size(basisset,1)
    error('dimensions of y and basis are off')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fit the coefficients once
coefs = calc_eigendecomposition(y,basisset) ; 
yrecon = basisset*coefs ; 

switch randwei
    case 'uniform'
        weivec = ones(dimtouse,1) ; 
    case 'power'
        weivec = abs(coefs).^2 ./ sum(abs(coefs).^2) ; 
    otherwise
        error('rand_wei not implemented')
end

% preallocate
surr_data = nan(size(basisset,1),nperms) ; 
surr_coefs = nan(dimtouse,nperms) ; 

% loop
for idx = 1:nperms

    disp(idx)

    % signflip the
    flipvec = ones(dimtouse,1) ; 
    % flipvec(randi(dimtouse,floor(dimtouse/2),1)) = -1 ; 
  
    % pick inds without replacement
    inds_shuff = datasample(...
        1:dimtouse,num2flip,'Replace',false,'Weights',weivec) ;

    % flip some signs
    flipvec(inds_shuff) = -1 ; 

    ncoefs = surr_logpolyfit(abs(coefs),'slide50') ; % samo inds frm norm dist
                                                     % w/ mu 50, to move inds 
                                                     % up or down

    % make the rand data
    surr_data(:,idx) = basisset*(ncoefs.*flipvec) ;
    % save the coeffs
    surr_coefs(:,idx) = ncoefs.*flipvec ; 
end

