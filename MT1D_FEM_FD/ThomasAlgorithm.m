function x=ThomasAlgorithm(A,d)
%https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
[N,M]=size(A);
if N~=M
    error('The number of rows is not equal to that of columns');
end
% c=zeros(1,N-1);
%c=A(N+1:N:N*N)
c=A(N+1:N+1:N*N);
u=zeros(N,1);
L=zeros(N,1);
y=zeros(N,1);
x=zeros(N,1);
u(1)=A(1,1);
L(2)=A(2,1)/u(1);
y(1)=d(1);
for j=2:N
    L(j)=A(j,j-1)/u(j-1);
    u(j)=A(j,j)-L(j)*c(j-1);
    y(j)=d(j)-L(j)*y(j-1);
end
x(N)=y(N)/u(N);
for j=N-1:-1:1
    x(j)=(y(j)-c(j)*x(j+1))/u(j);
end