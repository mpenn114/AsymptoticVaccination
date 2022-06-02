function Out = SIRODE(Input,U,beta1,beta2,beta3,beta4,mu1,mu2,N)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Note that the expected input into this function is of the form:
    %(V_1,V_2,...,V_n)
    %where
    %V_i = (S_i,I_i,R_i,S^V_i,I^V_i,R^V_i,W_i)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Groups = length(Input)/7;
    InfIndices = 2:7:(Groups*7 - 5);
    InfVIndices = 5:7:(Groups*7 - 2);
    if Groups ~= ceil(Groups)
        disp("Error - Input array does not have dimension divisible by 7")
    end
    if length(U) ~= Groups
        disp('Error - Vaccination policy does not have same length as number of groups')
    end
    Out = zeros(length(Input),1);
    for i = 1:Groups
        
        Out(7*(i-1) + 1) = -Input(7*(i-1) + 1)*(dot(beta1(i,:),Input(InfIndices))+dot(beta2(i,:),Input(InfVIndices))) - (Input(7*(i-1) + 1)*U(i))/(N(i) - Input(7*i));
        Out(7*(i-1) + 2) = Input(7*(i-1) + 1)*(dot(beta1(i,:),Input(InfIndices))+dot(beta2(i,:),Input(InfVIndices))) - mu1(i)*Input(7*(i-1) + 2);
        Out(7*(i-1) + 3) = mu1(i)*Input(7*(i-1) + 2);
        Out(7*(i-1) + 4) = -Input(7*(i-1) + 4)*(dot(beta3(i,:),Input(InfIndices))+dot(beta4(i,:),Input(InfVIndices))) + (Input(7*(i-1) + 1)*U(i))/(N(i) - Input(7*i));
        Out(7*(i-1) + 5) = Input(7*(i-1) + 4)*(dot(beta3(i,:),Input(InfIndices))+dot(beta4(i,:),Input(InfVIndices))) -mu2(i)*Input(7*(i-1) + 5);
        Out(7*(i-1) + 6) = mu2(i)*Input(7*(i-1) + 5);
        Out(7*(i-1) + 7) = U(i);
    end
end