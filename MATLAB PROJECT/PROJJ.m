f1=261.6255653;
fs1= 10*f1;
frmsz1=0.5*fs1;
t1=(0:1/fs1:0.5);
x1=cos(2*pi*f1*t1);

figure1=figure;
plot(t1,x1);

f2=293.6647679;
fs2=10*f2;
t2=(0:1/fs2:0.5);
x2=cos(2*pi*f2*t2);
figure2=figure;
plot(t2,x2);

f3=329.6275569;
fs3=10*f3;
frmsz3=0.5*fs3;
t2=(0:1/fs3:0.5);
x3=cos(2*pi*f3*t3);
figure3=figure;
plot(t3,x3);

f4=349.2282314;
fs4=10*f4;
t4=(0:1/fs4:0.5);
x4=cos(2*pi*f4*t4);
figure4=figure;
plot(t4,x4);


fst= 10*f4;
frmsz= round(0.5*fst) ;
T=(0:1:frmsz-1)*(1/fst);
x1=cos(2*pi*f1*T);
x2=cos(2*pi*f2*T);
x3=cos(2*pi*f3*T);
x4=cos(2*pi*f4*T);
xt=[x1,x2,x3,x4];
frmsz=round(4*frmsz);
T=(0:1:frmsz-1)*(1/fst);
figure5=figure;
plot(T,xt);

sound(xt,fst);
filename='sound.wav';
audiowrite(filename,xt,fst);



xt=[x1,x2,x3,x4];
F=(-frmsz/2:1:(frmsz/2)-1)*fst/frmsz;
XF= fft(xt);
figure6=figure;
plot(F,abs(fftshift(XF)));


N=length(xt);
E= sum(abs(xt).^2)*(1/fst);
E_p=sum((abs(XF).^2)/N)/fst;
display(E_p);
display(E);


[y, x] = butter(20 ,300/(fst/2));
y1_t = filter(y, x, xt);

figure7=figure;
freqz(y,x,F,fst);
sound(y1_t,fst);
filename='do.wav';
audiowrite(filename,y1_t,fst);

#figure8= figure;
#plot(y1_t);

N=length(y1_t);
EY1=sum((abs(y1_t)).^2)/fst;#timedomain
display(EY1);
Y1F=fft(y1_t);
figure9=figure;
plot(F,abs(fftshift(Y1F)));

EYP_1=sum((abs(Y1F).^2)/N)/fst;#parseval
display(EYP_1);
[b, a] = butter(20 ,2*300/fst,'high');
freqz(b,a,F,fst);
Y2T=filter(b,a,xt);
figure9=figure;
plot(T,Y2T);
sound(Y2T,fst);
filename='do2.wav';
audiowrite(filename,Y2T,fst);
N=length(Y2T);
EY_2=sum((abs(Y2T)).^2)*(1/fst);
display(EY_2);
F=(-frmsz/2:1:(frmsz/2)-1)*fst/frmsz;
Y2F= fft(Y2T);
figure10= figure;
plot(F,abs(fftshift(Y2F)));
EYP_2=sum((abs(Y2F).^2)/N)*(1/fst);
display(EYP_2);

