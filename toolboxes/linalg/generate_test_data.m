%Generate Test Data for linalg toolbox

%Test data for matrix-matrix mult
M = 100;
N = 234;
K = 40;

A = single(complex(randn(M,K),randn(M,K)));
B = single(complex(randn(K,N),randn(K,N)));
C1 = single(zeros(M,N));
C2 = single(A*B+C1);

write_mr_raw(A.','A.cplx');
write_mr_raw(B.','B.cplx');
write_mr_raw(C1.','C1.cplx');
write_mr_raw(C2.','C2.cplx');

