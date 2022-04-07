function Y = FFTm(X)    %recursive function
%The size of x needs to be a power of 2

  N=max(size(X)); 
  p=log2(N);  
  
  if p == 1
  Y =  DFT(X);
  return 
   
  else 
  Even  = FFTm(X(1:2:end-1));    
  Odd   = FFTm(X(2:2:end)); 
  
  Wm   = exp(-1i*2*pi/N);
  Y    = nan(1,N);

for i=0:N-1   
    if i<=N/2-1
       Y(i+1)= Even(i+1)+Wm^(i)*Odd(i+1);
    else
       Y(i+1)= Even(i-N/2+1)+Wm^(i)*Odd(i-N/2+1);
    end        
end
  end
 
  end
     
function series=DFT(X)    % DFT function

M = max(size(X));
Wm=exp(-1i*2*pi/M);
series=nan(1,M);
int=nan(1,M);
for k=0:M-1           
    for x=0:M-1    
        int(x+1)=Wm^(x*k)*X(x+1);
    end  
    series(k+1)=sum(int);      
end

end