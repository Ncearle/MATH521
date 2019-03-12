%If you don't have Octave or Matlab installed, you can also copy & paste the
%code of this file to an online installation of Octave, e.g.
%   https://octave-online.net/

%ADVECTION_DIFFUSION computes and plots the solution of the advection-
%diffusion equation
%       a u' - D u" = f on ]0,1[
%       u(0) = 0
%       u(1) = 0
%a and f are real numbers, D > 0 and N is the (integer) number of
%subintervals for the mesh.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                      %
% problem data                         %
a = 1; D = -0.1; f = 1; N = 25;          %
%                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xyxyxyxyxyxyxxyxxyxyxyxxyx = 1/N; xyxxxyxyxxyxxyxxyxyxyyyyx = N; xyxyxyxyxxxyxyxxxxyxyxyxyx = 2; xyxyxyxxyxyxxyxyyyxyxyyxyx = 1; xxyxyxyxyxyxyyxyxxyyxyyxyxyxyxyxxxx = 0; xyxyxyxyxyyxyyxyxxyxyyxyyxy = 1000; xyxyxyxyxyxxyyxyxyxxxxxyxyyy = linspace(xxyxyxyxyxyxyyxyxxyyxyyxyxyxyxyxxxx,xyxyxyxxyxyxxyxyyyxyxyyxyx,xyxxxyxyxxyxxyxxyxyxyyyyx+xyxyxyxxyxyxxyxyyyxyxyyxyx)'; xyxyxyxyxyxyxxyxxxyxyyyyyxyxxyxyxyxyx = xxyxyxyxyxyxyyxyxxyyxyyxyxyxyxyxxxx*ones(xyxxxyxyxxyxxyxxyxyxyyyyx-xyxyxyxxyxyxxyxyyyxyxyyxyx,xyxyxyxxyxyxxyxyyyxyxyyxyx)/xyxyxyxyxyyxyyxyxxyxyyxyyxy - xyxyxyxxyxyxxyxyyyxyxyyxyx + xyxyxyxyxxxyxyxxxxyxyxyxyx; xxxyxyxyxyxyxyxyyxyxyxyyyy = @(xyxyxyxyxyxxyyxyxyxxxxxyxyyy) f*ones(size(xyxyxyxyxyxxyyxyxyxxxxxyxyyy)); xyxxxyxxyxxyyyyyxxxxxx = a; xyxyxxxyxyyxyxyyy = D; xxxyxyxyxyxyyyxyyxyxy = xxxyxyxyxyxyxyxyyxyxyxyyyy(xyxyxyxyxyxxyyxyxyxxxxxyxyyy(xyxyxyxyxxxyxyxxxxyxyxyxyx:end-xyxyxyxxyxyxxyxyyyxyxyyxyx)); xxxyxyxyxyxyxyyyyxyxyxyxyxyyyyyx = xyxyxxxyxyyxyxyyy*spdiags([-xyxyxyxyxyxyxxyxxxyxyyyyyxyxxyxyxyxyx xyxyxyxyxxxyxyxxxxyxyxyxyx*xyxyxyxyxyxyxxyxxxyxyyyyyxyxxyxyxyxyx -xyxyxyxyxyxyxxyxxxyxyyyyyxyxxyxyxyxyx],-xyxyxyxxyxyxxyxyyyxyxyyxyx:xyxyxyxxyxyxxyxyyyxyxyyxyx,xyxxxyxyxxyxxyxxyxyxyyyyx-xyxyxyxxyxyxxyxyyyxyxyyxyx,xyxxxyxyxxyxxyxxyxyxyyyyx-xyxyxyxxyxyxxyxyyyxyxyyxyx)./xyxyxyxyxyxyxxyxxyxyxyxxyx^xyxyxyxyxxxyxyxxxxyxyxyxyx; xyxyxyxyyyyxyxxxxxxyyyxy = xxxyxyxyxyxyyyxyyxyxy;
if xyxxxyxxyxxyyyyyxxxxxx > xxyxyxyxyxyxyyxyxxyyxyyxyxyxyxyxxxx
  xyxxxyxxyxxxxyxxyxxyyyyy = xyxxxyxxyxxyyyyyxxxxxx*spdiags([-xyxyxyxyxyxyxxyxxxyxyyyyyxyxxyxyxyxyx xyxyxyxyxyxyxxyxxxyxyyyyyxyxxyxyxyxyx],[-xyxyxyxxyxyxxyxyyyxyxyyxyx xxyxyxyxyxyxyyxyxxyyxyyxyxyxyxyxxxx],xyxxxyxyxxyxxyxxyxyxyyyyx-xyxyxyxxyxyxxyxyyyxyxyyxyx,xyxxxyxyxxyxxyxxyxyxyyyyx-xyxyxyxxyxyxxyxyyyxyxyyxyx)/xyxyxyxyxyxyxxyxxyxyxyxxyx;
else
  xyxxxyxxyxxxxyxxyxxyyyyy = xyxxxyxxyxxyyyyyxxxxxx*spdiags([-xyxyxyxyxyxyxxyxxxyxyyyyyxyxxyxyxyxyx xyxyxyxyxyxyxxyxxxyxyyyyyxyxxyxyxyxyx],[xxyxyxyxyxyxyyxyxxyyxyyxyxyxyxyxxxx xyxyxyxxyxyxxyxyyyxyxyyxyx],xyxxxyxyxxyxxyxxyxyxyyyyx-xyxyxyxxyxyxxyxyyyxyxyyxyx,xyxxxyxyxxyxxyxxyxyxyyyyx-xyxyxyxxyxyxxyxyyyxyxyyxyx)/xyxyxyxyxyxyxxyxxyxyxyxxyx;
end
xyyyxyxxyxxxyxyxxyyxyxyxxyyx = xyxxxyxxyxxyyyyyxxxxxx*spdiags([-xyxyxyxyxyxyxxyxxxyxyyyyyxyxxyxyxyxyx xyxyxyxyxyxyxxyxxxyxyyyyyxyxxyxyxyxyx],[-xyxyxyxxyxyxxyxyyyxyxyyxyx xyxyxyxxyxyxxyxyyyxyxyyxyx],xyxxxyxyxxyxxyxxyxyxyyyyx-xyxyxyxxyxyxxyxyyyxyxyyxyx,xyxxxyxyxxyxxyxxyxyxyyyyx-xyxyxyxxyxyxxyxyyyxyxyyxyx)/(xyxyxyxyxxxyxyxxxxyxyxyxyx*xyxyxyxyxyxyxxyxxyxyxyxxyx);
yxxxyxxyxyxxyxxxyxyxxxxxyxy = f/xyxyxxxyxyyxyxyyy;
yyyxyyxyxyxyxyxyyxyy = xyxxxyxxyxxyyyyyxxxxxx/xyxyxxxyxyyxyxyyy;
xyxyxyxyxyxyxyxyxyxyxyxyxyxyyyy = linspace(xxyxyxyxyxyxyyxyxxyyxyyxyxyxyxyxxxx,xyxyxyxxyxyxxyxyyyxyxyyxyx,xyxyxyxyxyyxyyxyxxyxyyxyyxy);
xyyxyxyxxxyxyxyyyxyxyxyxyxyxyxyxxxxxyyyy = yxxxyxxyxyxxyxxxyxyxxxxxyxy/yyyxyyxyxyxyxyxyyxyy*xyxyxyxyxyxyxyxyxyxyxyxyxyxyyyy - yxxxyxxyxyxxyxxxyxyxxxxxyxy/yyyxyyxyxyxyxyxyyxyy*(exp(yyyxyyxyxyxyxyxyyxyy*(xyxyxyxyxyxyxyxyxyxyxyxyxyxyyyy-xyxyxyxxyxyxxyxyyyxyxyyxyx)) + (exp(yyyxyyxyxyxyxyxyyxyy*(xyxyxyxyxyxyxyxyxyxyxyxyxyxyyyy-xyxyxyxxyxyxxyxyyyxyxyyxyx))-xyxyxyxxyxyxxyxyyyxyxyyxyx)/(exp(yyyxyyxyxyxyxyxyyxyy)-1));
figure; plot(xyxyxyxyxyxxyyxyxyxxxxxyxyyy,[xxyxyxyxyxyxyyxyxxyyxyyxyxyxyxyxxxx;(xyxxxyxxyxxxxyxxyxxyyyyy+xxxyxyxyxyxyxyyyyxyxyxyxyxyyyyyx)\xyxyxyxyyyyxyxxxxxxyyyxy;xxyxyxyxyxyxyyxyxxyyxyyxyxyxyxyxxxx],'o-',xyxyxyxyxyxxyyxyxyxxxxxyxyyy,[xxyxyxyxyxyxyyxyxxyyxyyxyxyxyxyxxxx;(xyyyxyxxyxxxyxyxxyyxyxyxxyyx+xxxyxyxyxyxyxyyyyxyxyxyxyxyyyyyx)\xxxyxyxyxyxyyyxyyxyxy;xxyxyxyxyxyxyyxyxxyyxyyxyxyxyxyxxxx],'s-',xyxyxyxyxyxyxyxyxyxyxyxyxyxyyyy,xyyxyxyxxxyxyxyyyxyxyxyxyxyxyxyxxxxxyyyy); xlabel('{\itx}'); ylabel('{\itu}'); legend('upwind differencing','central differencing','analytical solution','location','northwest');
print -dpng figure.png
