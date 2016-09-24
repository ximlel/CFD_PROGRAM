line=10; %*3
column=30;

L_x=0.1;
L_y=0.1;


gamma=1.4;

delta_rho=0.0;

M_0=10;

rho_R=1;
u_R=1;
p_m=1/(gamma*M_0^2);

M_R=u_R/sqrt(gamma*p_m/rho_R);
p_m_d=p_m*(2*gamma*M_R^2-(gamma-1))/(gamma+1);
rho_R_d=rho_R*(gamma+1)*M_R^2/((gamma-1)*M_R^2+2);
u_R_d=u_R*rho_R/rho_R_d;

rho_L=rho_R+delta_rho;
u_L=M_R*sqrt(gamma*p_m/rho_L);
rho_L_d=rho_L*(gamma+1)*M_R^2/((gamma-1)*M_R^2+2);
u_L_d=u_L*rho_L/rho_L_d;



rho=zeros(column,1);
fid = fopen('RHO.txt','wt');
for j=1:line
    rho(1:(column/2))=rho_R;
    rho((column/2+1):column)=rho_R_d;
    fprintf(fid,'%12.10f\t',rho);
    fprintf(fid,'\n');

    rho(1:(column/2))=rho_L;
    rho((column/2+1):column)=rho_L_d;
    fprintf(fid,'%12.10f\t',rho);
    fprintf(fid,'\n');

    rho(1:(column/2))=rho_R;
    rho((column/2+1):column)=rho_R_d;
    fprintf(fid,'%12.10f\t',rho);
    fprintf(fid,'\n');
end
fclose(fid);


u=zeros(column,1);
fid = fopen('U.txt','wt');
for j=1:line
    u(1:(column/2))=u_R;
    u((column/2+1):column)=u_R_d;
    fprintf(fid,'%12.10f\t',u);
    fprintf(fid,'\n');

    u(1:(column/2))=u_L;
    u((column/2+1):column)=u_L_d;
    fprintf(fid,'%12.10f\t',u);
    fprintf(fid,'\n');

    u(1:(column/2))=u_R;
    u((column/2+1):column)=u_R_d;
    fprintf(fid,'%12.10f\t',u);
    fprintf(fid,'\n');
end
fclose(fid);


v=zeros(column,1);
fid = fopen('V.txt','wt');
for j=1:(line*3)
    fprintf(fid,'%12.10f\t',v);
    fprintf(fid,'\n');
end
fclose(fid);


p=p_m*ones(column,1);
fid = fopen('P.txt','wt');
for j=1:(line*3)
    p(1:(column/2))=p_m;
    p((column/2+1):column)=p_m_d;
    fprintf(fid,'%12.10f\t',p);
    fprintf(fid,'\n');
end
fclose(fid);


eps=1e-9;
t_all=3;
step=100000;

fid = fopen('config.txt','wt');
fprintf(fid,'%g\t',gamma);
fprintf(fid,'%g\t',t_all);
fprintf(fid,'%g\t',L_x);
fprintf(fid,'%g\t',L_y);
fprintf(fid,'%g\t',eps);
fprintf(fid,'%i\t',step);
fclose(fid);
