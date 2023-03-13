function y_shift = circshift_lininterp(y,l)

  N     = length(y);
  
% take shift modulo N to account for negative shifts and shifts exceeding vector size  
  l_eff = mod(l,N)
  
% enforce column vector
  SWAP = 0;
  if size(y,2) > size(y,1)
    SWAP = 1; % need's to be swapped again in the end
    y = y';
  end

% support points and weights  
  i1 = floor(l_eff);
  i2 = i1+1;

  w1 = i1+1-l_eff;
  w2 = 1 - w1;

% realize circshift via matrix multiplication  
  h_mat =  w1*circshift(speye(N),i1) + w2*circshift(speye(N),i2);
  
  full(h_mat)

  y_shift = h_mat * y;
  
  

  if SWAP == 1
      y_shift = y_shift';
  end
  
%{  
figure(1);clf;hold all;
plot(y,'ro-')
plot(y_shift,'bo-')
plot(circshift(y,floor(l)),'b--')
plot(circshift(y,ceil(l)),'b--')
  %}
  
  
  
  
  
return
%% START
N = length(a)

l_eff = mod(l,N);

%circshift_realval(X, lag)
i1 = floor(l_eff)
i2 = i1+1

w1 = i1+1-l_eff
w2 = 1 - w1

h = zeros(size(a))
h(i1+1) = w1
h(i2+1) = w2

a
crcshifta = circshift(a,floor(l))

if size(a,1) > size(a,2)
h_conv = [h(2:N);h(1:N)];
else
h_conv = [h(2:N),h(1:N)];
end
a_shift = conv(h_conv, a, 'valid')

figure(1);clf;hold all;
plot(a,'ro-')
plot(crcshifta,'go-')
plot(a_shift,'Bo-')

sum(a_shift)
sum(a)


%% matrix  
   h_mat =  w1*circshift(speye(N),i1) + w2*circshift(speye(N),i2)

   h_mat * a

return








a_shift2 = interp1(0:length(a)-1,a, mod([0:length(a)-1]-l, length(a)),'linear')


%%

  y = 0:19;
  l = 2;

% enforce row vector
  SWAP = 0;
  if size(y,1) > size(y,2)
    SWAP = 1; % need's to be swapped again in the end
    y = y';
  end
  
% periodic append
  N = length(y);
  a_mod  = [y(end),y];
  x_mod  = [-1,0:N-1];
  x_shft = mod(x_mod - l, N)
  y_shft = interp1(x_mod,y, x_shft,'linear')
  
  circshift(y,floor(l))
  

a
SWAP
