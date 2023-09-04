function M=calmem(X,c,L)
SUM=0;
for k=1:L-1
   SUM=SUM+X(L-k)*c(k);
end
M=SUM;