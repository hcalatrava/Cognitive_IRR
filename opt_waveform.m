function trSNRmatrix = opt_waveform(A, L, scat_coeff, EST_Phi, Sigma_power_half, sgm_sqr)

Amatrix = diag(A);
% SNRmatrix = inv(bigR*( eye(N) - EST_Phi' * pinv(EST_Phi*EST_Phi') * EST_Phi )*bigR') * (Amatrix * scat_coeff * EST_Phi) * (Amatrix * scat_coeff * EST_Phi)';
SNRmatrix = inv(Amatrix*Sigma_power_half*Sigma_power_half*Amatrix' + sgm_sqr*eye(L,L)) * (Amatrix * scat_coeff * EST_Phi) * (Amatrix * scat_coeff * EST_Phi)';
trSNRmatrix = -real(trace(SNRmatrix)) + (real(trace(SNRmatrix)/L)) * (real(trace(Amatrix * Amatrix')) - L);