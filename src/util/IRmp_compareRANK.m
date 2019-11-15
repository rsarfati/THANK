close all

[G1, C, impact, eu, SDX, zmat, NY, NX] = modelTHANK(postmode');


resp=-impact*SDX;
RESP(:,:,1)=resp;

for j=2:20
    RESP(:,:,j)=G1*RESP(:,:,j-1);
end

pv=postmode';
pv(36:40)=[0.001 1 0 0 0]';
[G1, C, impact, eu, SDX, zmat, NY, NX] = modelTHANK(pv);

respRANK=-impact*SDX;
RESPRANK(:,:,1)=respRANK;

for j=2:20
    RESPRANK(:,:,j)=G1*RESPRANK(:,:,j-1);
end


figure(1)
sgtitle('IR to a 1-std MP shock')

subplot(2,2,1); 
plot([0:19],squeeze(RESP(1,1,:))); hold on
plot([0:19],squeeze(RESPRANK(1,1,:))); 
legend('HANK','RANK')
title('response of y')

subplot(2,2,2); 
plot([0:19],squeeze(RESP(9,1,:))); hold on
plot([0:19],squeeze(RESPRANK(9,1,:)));
title('response of c')

subplot(2,2,3); 
plot([0:19],squeeze(RESP(13,1,:))); hold on
plot([0:19],squeeze(RESPRANK(13,1,:)));
title('response of i')

subplot(2,2,4); 
plot([0:19],squeeze(RESP(10,1,:))); hold on
plot([0:19],squeeze(RESPRANK(10,1,:)));
title('response of R')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
sgtitle('IR to a 1-std technology shock')

subplot(2,2,1); 
plot([0:19],-squeeze(RESP(1,2,:))-cumsum(squeeze(RESP(17,2,:)))); hold on
plot([0:19],-squeeze(RESPRANK(1,2,:))-cumsum(squeeze(RESPRANK(17,2,:)))); 
legend('HANK','RANK')
title('response of y')

subplot(2,2,2); 
plot([0:19],-squeeze(RESP(9,2,:))-cumsum(squeeze(RESP(17,2,:)))); hold on
plot([0:19],-squeeze(RESPRANK(9,2,:))-cumsum(squeeze(RESPRANK(17,2,:)))); 
title('response of c')

subplot(2,2,3); 
plot([0:19],-squeeze(RESP(13,2,:))-cumsum(squeeze(RESP(17,2,:)))); hold on
plot([0:19],-squeeze(RESPRANK(13,2,:))-cumsum(squeeze(RESPRANK(17,2,:)))); 
title('response of i')

subplot(2,2,4); 
plot([0:19],-squeeze(RESP(10,2,:))); hold on
plot([0:19],-squeeze(RESPRANK(10,2,:)));
title('response of R')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
sgtitle('IR to a 1-std investment shock')

subplot(2,2,1); 
plot([0:19],-squeeze(RESP(1,4,:))); hold on
plot([0:19],-squeeze(RESPRANK(1,4,:))); 
legend('HANK','RANK')
title('response of y')

subplot(2,2,2); 
plot([0:19],-squeeze(RESP(9,4,:))); hold on
plot([0:19],-squeeze(RESPRANK(9,4,:))); 
title('response of c')

subplot(2,2,3); 
plot([0:19],-squeeze(RESP(13,4,:))); hold on
plot([0:19],-squeeze(RESPRANK(13,4,:))); 
title('response of i')

subplot(2,2,4); 
plot([0:19],-squeeze(RESP(10,4,:))); hold on
plot([0:19],-squeeze(RESPRANK(10,4,:))); 
title('response of R')
    