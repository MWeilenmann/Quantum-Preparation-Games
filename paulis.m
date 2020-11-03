%-------------------------------------------------------------------------%
%      Return all combinations of local Pauli measurements for n parties. 
%      Each measurement is an element of a cell array (size 3^n, 2^n) where
%      the first index labels the inputs and the second one the outputs
%-------------------------------------------------------------------------%

function A = paulis(n)

% pauli measurements
 pX=[0,1;1,0];
 pY=[0,-1i;1i,0];
 pZ=[1,0;0,-1];
 
% povm from paulis
 pXYZ = cell(3,2);
 pXYZ(1,1) = mat2cell((eye(2)+pX)/2, [2]);
 pXYZ(1,2) = mat2cell((eye(2)-pX)/2, [2]);
 pXYZ(2,1) = mat2cell((eye(2)+pY)/2, [2]);
 pXYZ(2,2) = mat2cell((eye(2)-pY)/2, [2]);
 pXYZ(3,1) = mat2cell((eye(2)+pZ)/2, [2]);
 pXYZ(3,2) = mat2cell((eye(2)-pZ)/2, [2]);
 
 A = cell(3^n, 2^n);
 for i=1:3^n
     for j=1:2^n
         A(i,j) = pXYZ(mod(i-1,3)+1, mod(j-1,2)+1);   
     end
 end

 for k=2:n
     for i=1:3^n
         for j=1:2^n
             A(i,j) = mat2cell(kron(cell2mat(pXYZ(mod(floor((i-1)/(3^(k-1))),3)+1, mod(floor((j-1)/(2^(k-1))),2)+1)),cell2mat(A(i,j))), [2^k], [2^k]);
         end
     end
 end 
end