	figure (4)

for i=1:nt
%	plot (r(:,i),P(:,i),'r',"linewidth",2,ri(:,i),Pi(:,i),'b',"linewidth",2)
	plot (r(:,i),0.3*PovPinf(:,i),'r',"linewidth",2,r(:,i),3*rhonorm(:,i),'b',"linewidth",2,r(:,i),0.5*unorm(:,i),'g',"linewidth",2,
		r(:,i),gam(:,i),'+k',r(:,i),alpha(:,i),'rx')
	axis ([r(1,i) r(nr,i) -5 12]);
	rnow = Rb(outf*i);	
	text (rnow-0.1,-4,"AIR");
	text (rnow+0.05,-4,"WATER");
	hold on
	plot ([rnow,rnow],[7,-5],'k',"linewidth",1.7)
	hold off
	xlabel ("Radius (m)")
	legend ("normalised pressure","normalised density","normalised velocity","gamma","level set","nominal interface location","northeast")
	trm = t(i);
	vrb = sprintf("time = %e seconds",trm);
	title (vrb)
	pause (0.05)
%	b = num2str(i);
%	print ('-dsvg',"%s.svg",b)
	
endfor

%   max(PovPinf(1,1),PovPinf(nr,1))+1

%	axis ([r(1,1) r(nr,1) 0 max(PovPinf(1,1),PovPinf(nr,1))+1]);
%figure(5)
%for i=1:50
%%%%	plot (r(:,i),P(:,i),'r',"linewidth",2,ri(:,i),Pi(:,i),'b',"linewidth",2)
%
%	subplot(2,2,1),plot (r(:,i),PovPinf(:,i),'r',"linewidth",2)
%	axis ([r(1,1) r(nr,1) -2 12]);
%	xlabel ("Radius (m)")
%	ylabel ("normalised pressure")
%	trm = t(i);
%	vrb = sprintf("time = %e seconds",trm);
%	title (vrb)
%%
%	subplot(2,2,2),plot (r(:,i),rhonorm(:,i),'b',"linewidth",2)
%	axis ([r(1,1) r(nr,1) 0 1.2]);
%	xlabel ("Radius (m)")
%	ylabel ("normalised density")
%%
%	subplot(2,2,3),plot (r(:,i),unorm(:,i),'g',"linewidth",2)
%	axis ([r(1,1) r(nr,1) -5 12]);
%	xlabel ("Radius (m)")
%	ylabel ("normalised velocity")
%%
%	subplot(2,2,4),plot (r(:,i),gam(:,i),'k',"linewidth",2)
%	axis ([r(1,1) r(nr,1) 0 10]);
%	xlabel ("Radius (m)")
%	ylabel ("gamma")
%
%	pause (0.05)
%	b = num2str(i);
%	print ('-dsvg',"gfm%s.svg",b)
%	
%endfor
