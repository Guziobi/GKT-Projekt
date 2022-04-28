function res = ekonomi(NTU, Cmm, Ka, U, Tcin, Thin, tdrift, beta)

bosse = Ka./(U.*(Thin - Tcin).*tdrift.*beta);

%Derivatan av epsilon map NTU
dedNTU=(exp(-NTU+Cmm.*NTU)-2.*exp(-NTU+Cmm.*NTU).*Cmm+Cmm.^2.*exp(-NTU+Cmm.*NTU))./(1-Cmm.*exp(-NTU+Cmm.*NTU)).^2;

%fsolve ska lösa ekvationen de/dNTU = bosse eller de/dNTU - bosse = 0
res = dedNTU - bosse;

end