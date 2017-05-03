function[Solution,iteration_table,Error_val]=Gauss_Seidel_V1(U,Error_bound,F)
A=U;
%n1=length(U(:,1));
%C=[solution vector];
C=F;
n=length(C);
X=zeros(n,1);
Error_calc=1;

iteration=0;
while (Error_calc)>Error_bound
    iteration=iteration+1;
    Z=X;
    for i=1:n
        j=1:n;
        j(i)=[]; %this eliminates the unknown coefficient
        X_temp=X;
        X_temp(i)=[];% this eliminates the solution for the value being calculated for
        X(i)=(C(i)-sum(A(i,j)*X_temp))/(A(i,j));
    end 
    X_solution(:,iteration)=X;
    %Calculation of the error for the first unknown, this is due because
    %the error of this is usually higher than that of the other values
    Error_calc=abs((X(1,1)-Z(1,1))/(X(1,1)));
end 
iteration_table=[1:iteration;X_solution];
Solution=X;
Error_val=Error_calc;