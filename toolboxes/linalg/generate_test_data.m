%Generate Test Data for linalg toolbox

%Test data for matrix-matrix mult
M = 100;
N = 234;
K = 40;
avgs = 100;

A = single(complex(randn(M,K),randn(M,K)));
B = single(complex(randn(K,N),randn(K,N)));
C1 = single(zeros(M,N));
C2 = single(A*B+C1);

write_mr_raw(A.','A.cplx');
write_mr_raw(B.','B.cplx');
write_mr_raw(C1.','C1.cplx');
write_mr_raw(C2.','C2.cplx');

S = zeros(K);
for a=1:avgs,
    tmp_noise = complex(randn(K,1),randn(K,1));
    S = S + tmp_noise * tmp_noise';
end
clear tmp_noise;
S = S/avgs;

write_mr_raw(S.', 'S.cplx');
S_chol = chol(S,'lower');
write_mr_raw(S_chol.', 'S_chol.cplx');
S_chol_inv = inv(S_chol);
write_mr_raw(S_chol_inv.', 'S_chol_inv.cplx');
