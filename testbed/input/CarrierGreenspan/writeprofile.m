clear all
xdx=[0 130 160]
dx=[1 0.05 0.05];
i=1
x(1)=0;
while x(i)<150
    x(i+1)=x(i)+interp1(xdx,dx,x(i));
    i=i+1
end
h3=0.04*x;
fi=fopen('x3.dep','wt')
fprintf(fi,'%8.5f ',x);fprintf(fi,'\n')
fprintf(fi,'%8.5f ',x);fprintf(fi,'\n')
fprintf(fi,'%8.5f ',x);fprintf(fi,'\n')
fclose(fi)
fi=fopen('y3.dep','wt')
y=ones(size(x))*0;
fprintf(fi,'%8.5f ',y);fprintf(fi,'\n')
y=ones(size(x))*1;
fprintf(fi,'%8.5f ',y);fprintf(fi,'\n')
y=ones(size(x))*2;
fprintf(fi,'%8.5f ',y);fprintf(fi,'\n')
fclose(fi)
fi=fopen('h3.dep','wt')
fprintf(fi,'%8.5f ',h3);fprintf(fi,'\n')
fprintf(fi,'%8.5f ',h3);fprintf(fi,'\n')
fprintf(fi,'%8.5f ',h3);fprintf(fi,'\n')
fclose(fi)
length(x)