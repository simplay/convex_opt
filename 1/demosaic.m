function [demosaicedImg] = demosaic(mosaiced, Omega, lambda, iterations)
% DEMOSAIC demosaices a given mosaiced image by solving a convex opt.
%          problem relying on gradient descend using a Symmetric BC 
%          and a fixed numer of iteration.
%
% @author: MICHAEL SINGLE
%          08-917-445
%
% @param mosaiced M x N x 3 Color Image encoded as a mosaiced image
%        according to bayer filter theory.
% @param Omega bayer filter tensor M x N x 3
%        Omega(:,:,1) is RED color mask
%        Omega(:,:,2) is GREEN color mask
%        Omega(:,:,3) is BLUE color mask
% @param lambda [Float] regularization factor to control exactness of result and the smoothness 
%        when solving the demosaicing convex minimization problem.
% @param iterations [Integer] number of iterations
% @return demosaicedImg M x N x 3 Color Image.

    % learning rate for gradient descent
    alpha = 0.001; 
    
    % in order to avoid zero-division
    eps = 1e-6;
    
    % initial guesses 
    uR = mosaiced(:,:,1); 
    uG = mosaiced(:,:,2);
    uB = mosaiced(:,:,3);
    
    h = waitbar(0,'Progress Gradient Descend');
    for i = 1:iterations
       waitbar(i/iterations) 
       [uR] = gradDescendStep(uR, mosaiced(:,:,1), Omega(:,:,1), lambda, alpha, eps);
       [uG] = gradDescendStep(uG, mosaiced(:,:,2), Omega(:,:,2), lambda, alpha, eps);
       [uB] = gradDescendStep(uB, mosaiced(:,:,3), Omega(:,:,3), lambda, alpha, eps);
    end
    close(h); 
    demosaicedImg = mat2Img(uR,uG,uB);
end

function [u_next] = gradDescendStep(u, g, omega, lambda, alpha, eps)
% gradDescendStep perform one gradient descent iteration
% that uses a symmetric BC. The derivations of all expressions are
% documented in the report.
%
% @param u current iterative solution of one particular color channel
% @param g initial mosaiced image of color channel corresponding to given u.
% @param Omega bayer filter mask M x N for a certain color channel
% @param lambda [Float] regularization factor to control exactness of result and the smoothness 
%        when solving the demosaicing convex minimization problem.
% @param alpha learning rate for gradient descend.
% @param eps small perturbation to avoid division by zero issues.
% @return next iteration of color channel u.

       % U = u + BC: mirrored bondaries pad with 
       U = padarray(u,[2,2], 'symmetric');

       % construct tau based on image U
       % first non-boundary index is N=3 since U has a [2,2] boundary.
       
       % since we will need the value of tau_i-1,j and tau_i_j-1 we have to compute the value 
       % the value of tau at U's [1,1], we have to start to compute at
       % index M = N-1
       
       % u_ip1j_M_u_ij := u_i+1,j - u_i,j
       u_ip1j_M_u_ij = U(3:end,2:end-1) - U(2:end-1,2:end-1);
       
       % u_ijp1_M_u_ij := u_i,j+1 - u_i,j
       u_ijp1_M_u_ij = U(2:end-1,3:end) - U(2:end-1,2:end-1);
       
       % tau2 := tau squared
       tau2 = u_ip1j_M_u_ij.^2 + u_ijp1_M_u_ij.^2;
       tau = (tau2+eps).^(0.5);

       % Compute d||grad(u)||_2/du_i,j = d tau_i,j /d u_i,j + d tau_i-1,j/d u_i,j + d tau_i,j-1/ d u_i,j
       dtau_ij_du_ij = (2*U(3:end-2,3:end-2) - U(4:end-1,3:end-2) - U(3:end-2,4:end-1))./tau(2:end-1,2:end-1);
       dtau_im1j_du_ij = (U(3:end-2,3:end-2) - U(2:end-3,3:end-2))./tau(1:end-2,2:end-1);
       dtau_ijm1_du_ij = (U(3:end-2,3:end-2) - U(3:end-2,2:end-3))./tau(2:end-1,1:end-2);
       d_grad_u_du_ij = dtau_ij_du_ij + dtau_im1j_du_ij + dtau_ijm1_du_ij;

       % gradient descend update rule
       u_next = u - alpha*(lambda*(omega.*(u-g)) + d_grad_u_du_ij);
end