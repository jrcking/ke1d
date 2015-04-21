	figure (3)
SA_ro=load("../SEDOV/SEDOV_solver/den.dat");
SA_u=load("../SEDOV/SEDOV_solver/vel.dat");
SA_E=load("../SEDOV/SEDOV_solver/ener.dat");
SA_p=load("../SEDOV/SEDOV_solver/pres.dat");


subplot(3,1,1),plot(r(:,nt),rho(:,nt),'rx',SA_ro(:,1),SA_ro(:,2),'k',"linewidth",1.5)
trm = t(nt);
vrb = sprintf("time = %e seconds",trm);
%title (vrb)
legend("numerical","exact")
grid
ylabel("density")
axis([r(1,nt) r(nr,nt) 0 1100]); axis "autoy";
subplot(3,1,2),plot(r(:,nt),u(:,nt),'rx',SA_u(:,1),SA_u(:,2),'k',"linewidth",1.5)
grid
ylabel("velocity")
axis([r(1,nt) r(nr,nt) 0 40]); axis "autoy";
subplot(3,1,3),plot(r(:,nt),P(:,nt),'rx',SA_p(:,1),SA_p(:,2),'k',"linewidth",1.5)
grid
ylabel("pressure")
axis([r(1,nt) r(nr,nt) 1e5 9e9]); axis "autoy";
xlabel("radius")


subplot(3,1,1),plot(SA_ro(:,1),SA_ro(:,2),'k',"linewidth",1.5)
trm = t(nt);
vrb = sprintf("time = %e seconds",trm);
%title (vrb)
grid
ylabel("density")
axis([r(1,nt) r(nr,nt) 0 1100]); axis "autoy";
subplot(3,1,2),plot(SA_u(:,1),SA_u(:,2),'k',"linewidth",1.5)
grid
ylabel("velocity")
axis([r(1,nt) r(nr,nt) 0 40]); axis "autoy";
subplot(3,1,3),plot(SA_p(:,1),SA_p(:,2),'k',"linewidth",1.5)
grid
ylabel("pressure")
axis([r(1,nt) r(nr,nt) 1e5 9e9]); axis "autoy";
xlabel("radius")
