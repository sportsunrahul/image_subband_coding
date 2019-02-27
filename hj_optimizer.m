
N=1024;
h0 = firhalfband(12, 0.4);
h1 = h0;
h1(2:2:end) = -1*h1(2:2:end);
H0 = fft(h0,N);
H1 = fft(h1,N);
figure(1);
[H0,w0] = freqz(h0);
[H1,w1] = freqz(h1);
l = length(H0);
w = 0:(pi/2)/(l-1):pi/2;
plot(w,abs(H0).^2 + abs(H1).^2);
ylabel("|H0(w)|^{2} + |H1(w)|^{2}");
title("How well the constrained |H0(w)|^{2} + |H1(w)|^{2} = 1 is before optimization?");

E=0;
for i = 1:N/4
     E = E + 2*((abs(H0(i))^2) + (abs(H1(i))^2) - 1)^2;
end

ws=N/4 + 1;
for i = ws:N/2
     E = E + 0.5*(abs(H0(i))^2);
end
figure(2);
freqz(h0);
hold on;
freqz(h1);
hold off;
title("Filter Before Optimization ");

StepSize = 1*ones(length(h0),1);
MinStepSize = StepSize.*1e-8;
half = 7;

[X,BestF,Iters] = hookejeeves((length(h0)-1)/2, h0(1:half), StepSize, MinStepSize, 1e-8, 1000, 'fx1');
h0(1:half)=X;
h0(half+1:end)=flip(X(1:end-1));
h1 = h0;
h1(2:2:end) = -1*h1(2:2:end);
figure(3);
freqz(h0);
hold on;
freqz(h1);
hold off;
title("Filter After Optimization");



figure(4);
[H0,w0] = freqz(h0);
[H1,w1] = freqz(h1);
l = length(H0);
w = 0:(pi/2)/(l-1):pi/2;
plot(w,abs(H0).^2 + abs(H1).^2);
ylabel("|H0(w)|^{2} + |H1(w)|^{2}");
title("How well the constraint |H0(w)|^{2} + |H1(w)|^{2} = 1 is fulfilled? ");


function y = fx1(h,k)
    half = 7;
    h0 = ones(13,1);
    h0(1:half)=h;
    h0(half+1:end)=flip(h(1:end-1));
    h1 = h0;
    h1(2:2:end) = -1*h1(2:2:end);
    N=1024;
    H0 = fft(h0,N);
    H1 = fft(h1,N);
    E=0;
    a = 0.1;
    for i = 1:N/4
         E = E + 2*((abs(H0(i))^2) + (abs(H1(i))^2) - 1)^2;
    end
    E = (1-a)*E;
    ws=N/4 + 1;
    for i = ws:N/2
         E = E + a*(abs(H0(i))^2);
    end
    y=E;
end

function [X,BestF,Iters] = hookejeeves(N, X, StepSize, MinStepSize, Eps_Fx, MaxIter, myFx)
% Function HOOKEJEEVS performs multivariate optimization using the
% Hooke-Jeeves search method.
%
% Input
%
% N - number of variables
% X - array of initial guesses
% StepSize - array of search step sizes
% MinStepSize - array of minimum step sizes
% Eps_Fx - tolerance for difference in successive function values
% MaxIter - maximum number of iterations
% myFx - name of the optimized function
%
% Output
%
% X - array of optimized variables
% BestF - function value at optimum
% Iters - number of iterations
%
Xnew = X;
BestF = feval(myFx, Xnew, N);
LastBestF = 100 * BestF + 100;
bGoOn = true;
Iters = 0;
while bGoOn
  Iters = Iters + 1;
  if Iters > MaxIter
    break;
  end
  X = Xnew;
  for i=1:N
    bMoved(i) = 0;
    bGoOn2 = true;
    while bGoOn2
      xx = Xnew(i);
      Xnew(i) = xx + StepSize(i);
      F = feval(myFx, Xnew, N);
      if F < BestF
        BestF = F;
        bMoved(i) = 1;
      else
        Xnew(i) = xx - StepSize(i);
        F = feval(myFx, Xnew, N);
        if F < BestF
          BestF = F;
          bMoved(i) = 1;
        else
          Xnew(i) = xx;
          bGoOn2 = false;
        end
      end
    end
  end

  bMadeAnyMove = sum(bMoved);
  if bMadeAnyMove > 0
    DeltaX = Xnew - X;
    lambda = 0.5;
    lambda = linsearch(X, N, lambda, DeltaX, myFx);
    Xnew = X + lambda * DeltaX;
  end
  
  BestF = feval(myFx, Xnew, N);
  % reduce the step size for the dimensions that had no moves
  for i=1:N
    if bMoved(i) == 0
      StepSize(i) = StepSize(i) / 2;
    end
  end
  if abs(BestF - LastBestF) < Eps_Fx
    break
  end
  
  LastBest = BestF;
  bStop = true;
  for i=1:N
    if StepSize(i) >= MinStepSize(i)
      bStop = false;
    end
  end
  
  bGoOn = ~bStop;
end
end

function y = myFxEx(N, X, DeltaX, lambda, myFx)
  X = X + lambda * DeltaX;
  y = feval(myFx, X, N);
end

function lambda = linsearch(X, N, lambda, D, myFx)
  MaxIt = 100;
  Toler = 0.000001;
  iter = 0;
  bGoOn = true;
  while bGoOn
    iter = iter + 1;
    if iter > MaxIt
      lambda = 0;
      break
    end

    h = 0.01 * (1 + abs(lambda));
    f0 = myFxEx(N, X, D, lambda, myFx);
    fp = myFxEx(N, X, D, lambda+h, myFx);
    fm = myFxEx(N, X, D, lambda-h, myFx);
    deriv1 = (fp - fm) / 2 / h;
    deriv2 = (fp - 2 * f0 + fm) / h ^ 2;
    diff = deriv1 / deriv2;
    lambda = lambda - diff;
    if abs(diff) < Toler
      bGoOn = false;
    end
  end
end