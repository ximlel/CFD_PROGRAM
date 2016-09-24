line=20;
column=800;
shock=20;

gamma=1.4;
M=6;

f=1/(2/(gamma+1)/M/M+(gamma-1)/(gamma+1));
g=2*gamma/(gamma+1)*M*M-(gamma-1)/(gamma+1);

rho_R=1.4
u_R=0
p_R=1
rho_L=rho_R*f
u_L=(1-1/f)*M
p_L=p_R*g

rho=zeros(column,1);
for i=1:shock
    rho(i)=rho_L;
end
for i=(shock+1):column
    rho(i)=rho_R;
end
fid = fopen('RHO.txt','wt');
for j=1:line
fprintf(fid,'%12.10f\t',rho);
fprintf(fid,'\n');
end
fclose(fid);

u=zeros(column,1);
for i=1:shock
    u(i)=u_L;
end
for i=(shock+1):column
    u(i)=u_R;
end
fid = fopen('U.txt','wt');
for j=1:line
fprintf(fid,'%12.10f\t',u);
fprintf(fid,'\n');
end
fclose(fid);

v=zeros(column,1);
fid = fopen('V.txt','wt');
for j=1:line
fprintf(fid,'%12.10f\t',v);
fprintf(fid,'\n');
end
fclose(fid);

p=zeros(column,1);
for i=1:shock
    p(i)=p_L;
end
for i=(shock+1):column
    p(i)=p_R;
end
fid = fopen('P.txt','wt');
for j=1:line
fprintf(fid,'%12.10f\t',p);
fprintf(fid,'\n');
end
fclose(fid);
