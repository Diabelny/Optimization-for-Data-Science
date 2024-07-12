data = load("Data_Energy_community\inputs.txt");
training=data';
predictions=load("Data_Energy_community\groundtruth.txt");
predictions=predictions';
N=size(predictions);
N=2*N(1);
eps=10^-2;% iperparametro formula SVR
epsilon=10^-2;
kernel='Radiale';
max_iter=1000;
tol=10^-12;
%array_C=logspace(-2,3,6);
%array_sigma=logspace(-3,3,7);
param_kernel=100;
C=10^-2;
lambda_0=C*ones(1,N);
%lambda_0=C/2*ones(1,N);
%lambda_0=zeros(1,N);

main(training,predictions,C,N,eps,epsilon,kernel,param_kernel,lambda_0,max_iter,tol);
mainL1(training,predictions,C,N,eps,epsilon,kernel,param_kernel,lambda_0,max_iter,tol);
      

