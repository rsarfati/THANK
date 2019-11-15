
[G1, C, impact, eu, SDX, zmat, NY, NX] = modelTHANK(postmode');

Vstates=dlyapgio(G1,impact*SDX*SDX'*impact');
Vobs=diag(zmat*Vstates*zmat');

SHAREobs=zeros(7,7);

for i=1:7
    STDS=zeros(7); STDS(i,i)=SDX(i,i);
    Vstates=dlyapgio(G1,impact*STDS*STDS'*impact');
    SHAREobs(:,i)=diag(zmat*Vstates*zmat')./Vobs;
end



    