line=160;
column=40;
shock=1;

rho_L=1.4
u_L=8
p_L=1
rho_R=1.4
u_R=8
p_R=1

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
