close all

[G1, C, impact, eu, SDX, zmat, NY, NX] = modelTHANK(postmode');

Hetero=1;

resp=-impact*SDX;
RESP(:,:,1)=resp;

for j=2:20
    RESP(:,:,j)=G1*RESP(:,:,j-1);
end


figure(1)
sgtitle('IR to a 1-std MP shock')

subplot(2,2,1); 
plot([0:19],squeeze(RESP(1,1,:))); hold on
if Hetero==1; plot([0:19],squeeze(RESP(39,1,:))); end
legend('y','yH')
title('response of y and yH')

subplot(2,2,2); 
plot([0:19],squeeze(RESP(9,1,:)));
title('response of c')

subplot(2,2,3); 
plot([0:19],squeeze(RESP(13,1,:)));
title('response of i')

subplot(2,2,4); 
plot([0:19],squeeze(RESP(10,1,:)));
title('response of R')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
sgtitle('IR to a 1-std technology shock')

subplot(2,2,1); 
plot([0:19],-squeeze(RESP(1,2,:))-cumsum(squeeze(RESP(17,2,:)))); hold on
if Hetero==1; plot([0:19],-squeeze(RESP(39,2,:))-cumsum(squeeze(RESP(17,2,:)))); end
legend('y','yH')
title('response of y and yH')

subplot(2,2,2); 
plot([0:19],-squeeze(RESP(9,2,:))-cumsum(squeeze(RESP(17,2,:))));
title('response of c')

subplot(2,2,3); 
plot([0:19],-squeeze(RESP(13,2,:))-cumsum(squeeze(RESP(17,2,:))));
title('response of i')

subplot(2,2,4); 
plot([0:19],-squeeze(RESP(10,2,:)));
title('response of R')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
sgtitle('IR to a 1-std investment shock')

subplot(2,2,1); 
plot([0:19],-squeeze(RESP(1,4,:))); hold on
if Hetero==1; plot([0:19],-squeeze(RESP(39,4,:))); end
legend('y','yH')
title('response of y and yH')

subplot(2,2,2); 
plot([0:19],-squeeze(RESP(9,4,:)));
title('response of c')

subplot(2,2,3); 
plot([0:19],-squeeze(RESP(13,4,:)));
title('response of i')

subplot(2,2,4); 
plot([0:19],-squeeze(RESP(10,4,:)));
title('response of R')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
sgtitle('IR to a 1-std investment shock')

subplot(3,2,1); 
plot([0:19],-squeeze(RESP(1,4,:))); hold on
if Hetero==1; plot([0:19],-squeeze(RESP(39,4,:))); end
legend('y','yH')
title('response of y and yH')

subplot(3,2,2); 
plot([0:19],-squeeze(RESP(9,4,:)));
title('response of c')

subplot(3,2,3); 
plot([0:19],-squeeze(RESP(13,4,:)));
title('response of i')

subplot(3,2,4); 
plot([0:19],-squeeze(RESP(10,4,:)));
title('response of R')

subplot(3,2,5); 
plot([0:19],-squeeze(RESP(39,4,:)));
title('response of cH')

subplot(3,2,6); 
plot([0:19],-squeeze(RESP(40,4,:)));
title('response of cS')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)
sgtitle('IR to a 1-std G shock')

subplot(3,2,1); 
plot([0:19],-squeeze(RESP(1,3,:))); hold on
if Hetero==1; plot([0:19],-squeeze(RESP(39,3,:))); end
legend('y','yH')
title('response of y and yH')

subplot(3,2,2); 
plot([0:19],-squeeze(RESP(9,3,:)));
title('response of c')

subplot(3,2,3); 
plot([0:19],-squeeze(RESP(13,3,:)));
title('response of i')

subplot(3,2,4); 
plot([0:19],-squeeze(RESP(10,3,:)));
title('response of R')

subplot(3,2,5); 
plot([0:19],-squeeze(RESP(39,3,:)));
title('response of cH')

subplot(3,2,6); 
plot([0:19],-squeeze(RESP(40,3,:)));
title('response of cS')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(6)
sgtitle('IR to a 1-std b shock')

subplot(3,2,1); 
plot([0:19],-squeeze(RESP(1,7,:))); hold on
if Hetero==1; plot([0:19],-squeeze(RESP(39,7,:))); end
legend('y','yH')
title('response of y and yH')

subplot(3,2,2); 
plot([0:19],-squeeze(RESP(9,7,:)));
title('response of c')

subplot(3,2,3); 
plot([0:19],-squeeze(RESP(13,7,:)));
title('response of i')

subplot(3,2,4); 
plot([0:19],-squeeze(RESP(10,7,:)));
title('response of R')

subplot(3,2,5); 
plot([0:19],-squeeze(RESP(39,7,:)));
title('response of cH')

subplot(3,2,6); 
plot([0:19],-squeeze(RESP(40,7,:)));
title('response of cS')


    