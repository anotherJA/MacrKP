integer xmax,ymax,zmax,Nc_init,tmax
	
parameter (xmax=3e3,ymax=3e3,zmax=3e3,Nc_init=9999,tmax=8250)
	
INTEGER :: idum,num_agg(1:Nc_init),cell_macr(1:Nc_init,1:2),&
&ntotal(1:Nc_init),color_number,critical_size_data,ncritical,nc_cells,nc_macr&
&,nc_macr0,variable_name1,variable_name0,variable_name2,num

REAL res,rapport,probacellcr,&
	&radius_well,rad,xinit,yinit,zinit,g,rc,amp_noise,x,y,z,&
	&adh,lproteins,radiusi21,x0,y0,z0,x1,y1,z1,nf0,PI,proba_cell,proba_macr,tauagg,&
	&friction1,friction2,friction0,lf,number_corr,distancefusion,lpenetration,f0,f1
	
REAL, DIMENSION(3)  ::	gf,nf,normf,af,new_v,newf_pos
REAL, DIMENSION(Nc_init)  :: pos_x,pos_y,pos_z,v_x,v_y,v_z,radius,radiusaggregate
REAL, DIMENSION(Nc_init,3)  ::	velocity_previous,sum_adh_force
REAL :: mass(1:Nc_init),mass_c,min_relative_distance,max_relative_distance,&
&max_rad,cell_macr_ttl,ttl_area,friction,mass_m
REAL, DIMENSION(Nc_init,Nc_init)  :: relative_distance	
REAL :: positions(Nc_init,3),nmacrtable(10),&
&variable_name,&
&ncelltable(10),phi0,phi1,&
&radius_cell,&
&deltaG0,zmax_init,layer_macr,a,friction_substrat,steps_persistence&
&,init_rad_c,init_rad_m,vitesse
CHARACTER*16 filename
character*8 :: date
character*10 :: time
character*5 :: zone
integer,dimension(8) :: values,values_fin
character(len=8) :: fmt ! format descriptor
character*4 :: x111,x222
character*4 :: x333

call date_and_time(date,time,zone,values)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Definition of the parameters
critical_size_data=10
PI=4.D0*DATAN(1.D0)
nf0=1e3
radius_well=xmax
g=10*1e6*3600*3600;
lproteins=1
friction0=0.0005;friction=0;friction1=0;friction2=30;rc=10;number_corr=0 
Din=1.0;Dext=5.0;finfty=1.0;ncritical=10
ncelltable(1)=156;ncelltable(2)=312;ncelltable(3)=625;ncelltable(4)=1250;ncelltable(5)=2500;
ncelltable(6)=2500;ncelltable(7)=2500;ncelltable(8)=2500;ncelltable(9)=2500;
nmacrtable(1)=0;nmacrtable(2)=0;nmacrtable(3)=0;nmacrtable(4)=0;nmacrtable(5)=0;
nmacrtable(6)=625;nmacrtable(7)=1250;nmacrtable(8)=2500;nmacrtable(9)=5000;
lf=20.0
distancefusion0=10;distancefusion=distancefusion0
dt=0.02
num=9
tau_growth=0.03;timeadh=10;timeadhm=10;tdelayc=20;tdelaym=0;tauagg=0.1
rep_macr=0.3
init_rad_c=2.5e3;
factorg=1.7
lpenetration=200
fcrit=0.3
proba_cell=25;proba_macr=25
friction_cell=3.6e-6;friction_substrat=3.6e-6;friction_macr=0.1
vitesse=60;
amp_noise=1;amp_noise2=friction_macr*22*vitesse;
deltaG0=3.8;
phi0=0.25;phi1=0.9;phimm=0.75;
radius_cell=7.7/phi1**(1./3.);radius_macr=10.2;rmMM=radius_macr/phimm**(1./3.);rmm=radius_macr/phimm**(1./3.);
rcMM=radius_cell/(phi0)**(1./3.);rcm=radius_cell/(phi1)**(1./3.)
rapportc=0.05;rapportm=0.08
steps_persistence=10
mass_c=(rapportc)*(1e3)*(1e-18)*(4*3.14159/3.)*(radius_cell**3);mass_m=(rapportm)*(1e3)*(1e-18)*(4*3.14159/3.)*(radius_macr**3) 
probacellcr=proba_cell
zmax_init=0;layer_macr=radius_macr
variable_name1=int(friction_macr*100);variable_name2=int(deltaG0*100)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! The simulation is run for different initial quantities of cells

do 51 i51=num,num

nc_cells=ncelltable(i51)
nc_macr=nmacrtable(i51)
init_rad_m=init_rad_c/(2500.)**(1./3.)*(real(nc_macr))**(1./3.)
variable_name0=int(nc_macr)
print*, nc_macr

	call full_program(Nc_init,radius,rcMM,rcm,rmMM,rmm,mass_c,mass_m,mass,positions,cell_macr,sum_adh_force&
&,velocity_previous,ntotal,num_agg,proba_cell,proba_macr,dt,g,amp_noise,lproteins,&
&adh,nf0,lf,critical_size_data,values,radius_well,friction_cell,friction_substrat,friction_macr,tmax,&
&pi,fcrit,Din,Dext,lpenetration,finfty&
&,tau_growth,ncritical,nc_cells,nc_macr,distancefusion,timeadh,timeadhm,tdelayc,tdelaym,variable_name0,&
&variable_name1,variable_name2,radius_cell,factorg,amp_noise2,rep_macr,probacellcr,deltaG0,phimm,tauagg,&
&radius_macr,zmax_init,layer_macr,steps_persistence,init_rad_c,init_rad_m)

51 continue
	 END
	 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine full_program(Nc_init,radius,rcMM,rcm,rmMM,rmm,mass_c,mass_m,mass,positions,cell_macr,sum_adh_force&
&,velocity_previous,ntotal,num_agg,proba_cell,proba_macr,dt,g,amp_noise,lproteins,&
&adh,nf0,lf,critical_size_data,values,radius_well,friction0,friction1,friction2,tmax,&
&pi,fcrit,Din,Dext,d,finfty&
&,tau_growth,ncritical,nc_cells,nc_macr,distancefusion,timeadh,timeadhm,tdelayc,tdelaym,variable_name0,&
&variable_name1,variable_name2,radius_cell,factor,amp_noise2,rep_macr,probacellcr,deltaG0,phimm,tauagg,&
&radius_macr,zmax,layer_macr,steps_persistence,init_rad_c,init_rad_m)

integer Nc_init,tmax,critical_size_data,values(8),ncritical,nc_cells,nc_macr,&
&variable_name0,variable_name1,variable_name2
real, dimension(Nc_init) :: radius,mass
real, dimension(Nc_init,3) :: positions,sum_adh_force,velocity_previous
real, dimension(Nc_init,2) :: adh_strength
integer, dimension(Nc_init,2) :: cell_macr
integer, dimension(Nc_init) :: ntotal,num_agg
real :: mass_c,radius_well,proba_cell,proba_macr,dt,g,amp_noise,lproteins,adh,nf0,lf,probacellcr,&
&number_corr,friction0,friction1,friction2,pi,fcrit,Din,Dext,d,finfty,tau_growth,distancefusion,phimm&
&,timeadh,timeadhm,tdelayc,tdelaym,radius_cell,amp_noise2,deltaG0&
&,tauagg,radius_macr,zmax,layer_macr,mass_m,steps_persistence,init_rad_c,init_rad_m
character*1000 :: datadat
character(len=8) :: fmt ! format descriptor
character*5 :: x111,x222,x333

fmt = '(I5.5)' ! an integer of width 5 with zeros at the left

write (x111,fmt) variable_name0 ! converting integer to string using a 'internal file'
write (x222,fmt) variable_name1
write (x333,fmt) variable_name2 

datadat='data'//x111//''//x222//''//x333//'.dat'

call initialization(Nc_init,radius,rcMM,mass_c,mass_m,mass,positions,cell_macr,sum_adh_force,velocity_previous,ntotal&
&,radius_well,nc_cells,nc_macr,adh_strength,proba_macr,rmm,init_rad_c,init_rad_m)

call time_process(Nc_init,tmax,mass,radius,positions,cell_macr,num_agg,ntotal,proba_cell,proba_macr,dt,g,amp_noise,&
&lproteins,adh,nf0,lf,critical_size_data,number_corr,friction0,friction1,friction2,rcMM,rcm,rmMM,rmm,radius_well,&
&pi,fcrit,Din,Dext,d,finfty&
&,tau_growth,ncritical,nc_cells,nc_macr,distancefusion,datadat,adh_strength,timeadh,timeadhm,tdelayc,tdelaym,mass_c,mass_m&
&,radius_cell,factor,amp_noise2,rep_macr,probacellcr,deltaG0,phimm,tauagg,radius_macr,layer_macr,steps_persistence)

call final_routine(Nc_init,positions,radius,ntotal,cell_macr,number_corr,values,nc_cells,nc_macr,datadat)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine initialization(Nc_init,radius,rc,mass_c,mass_m,mass,positions,cell_macr,sum_adh_force,velocity_previous,ntotal,&
&radius_well,nc_cells,nc_macr,adh_strength,proba_macr,rmm,init_rad_c,init_rad_m)

integer Nc_init,nc_cells,nc_macr
real, dimension(Nc_init) :: radius,mass
real, dimension(Nc_init,3) :: positions,sum_adh_force,velocity_previous
real, dimension(Nc_init,2) :: adh_strength
integer, dimension(Nc_init,2) :: cell_macr
integer, dimension(Nc_init) :: ntotal
real :: rc,mass_c,radius_well,proba_macr,rmm,&
&mass_m,init_rad_c,init_rad_m


DO 2 i2=1,nc_cells+nc_macr	
		cell_macr(i2,1)=1;cell_macr(i2,2)=0
		sum_adh_force(i2,1)=0;sum_adh_force(i2,2)=0;sum_adh_force(i2,3)=0
		velocity_previous(i2,1)=0;velocity_previous(i2,2)=0;velocity_previous(i2,3)=0
		radius(i2)=rc
		mass(i2)=mass_c
		adh_strength(i2,1)=0
		if (i2.gt.nc_cells) then
			cell_macr(i2,1)=0;cell_macr(i2,2)=1
			adh_strength(i2,2)=proba_macr
			radius(i2)=rmm
			mass(i2)=mass_m
		end if
		
		ntotal(i2)=1

2	 CONTINUE

	 !!!!!!!!!! Initial position
	 
	 CALL initial_pos(Nc_init,radius_well,positions,nc_cells,nc_macr,radius,init_rad_c,init_rad_m)
	 
!print*, nc_cells,nc_macr
	
return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	 SUBROUTINE initial_pos(Nc_init,radius_well,positions,nc_cells,nc_macr,radius,init_rad_c,init_rad_m)
	 integer 

	 integer :: Nc_init,idum,nc_cells,nc_macr
	 REAL, DIMENSION(Nc_init,3) :: positions
	 real :: radius_well,xinit,yinit,zinit,radius(Nc_init),init_rad_c,init_rad_m,radius_init,radiusagg
	 
	 character*8 :: date
	 character*10 :: time
	 character*5 :: zone
	 integer,dimension(8) :: values
	 call date_and_time(date,time,zone,values)
	 idum=(values(1)+values(2)+values(3)+values(4)+values(5)+values(6)+values(7)+values(8))
	 DO 1000 i1000=1,nc_cells+nc_macr

			idum=idum+100

10001		if (i1000.gt.nc_cells) then
			radius_init=min(init_rad_m,2.9e3)
			elseif (i1000.le.nc_cells) then
			radius_init=min(init_rad_c,2.9e3)
			end if
			
			xinit=(2*ran2(idum)-1)*radius_init
			idum=idum+100
			yinit=(2*ran2(idum)-1)*radius_init
			idum=idum+100
			zinit=3e3+(2*ran2(idum)-1)*radius_init			
			idum=idum+100		
			radiusagg=radius(i1000)
			
			if (sqrt(xinit**2+yinit**2+(zinit-3e3)**2).gt.radius_init) then
				goto 10001
			end if
			
		positions(i1000,1)=xinit
		positions(i1000,2)=yinit
		positions(i1000,3)=zinit
		
1000 CONTINUE

	 RETURN
	 END
	 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	  FUNCTION ran2(idum)
      INTEGER :: idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,&
     &IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,&
     &NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
!  (C) Copr. 1986-92 Numerical Recipes Software #$!5,5.){2p491&&k"15.  page 272
      END	
	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	subroutine friction_function_new(Nc_init,agg1,friction,friction0,friction1,friction2,cell_macr,radius_well,positions,lf,radius,rc&
	&,adh_strength,probacellcr,rcMM,rcm,radiusmacr,timeadh,proba_cell,tdelayc,dt,phimm,layer_macr,radius_macr)

	integer :: Nc_init,cell_macr(Nc_init,2),agg1
	real :: friction1,friction2,radius_well,positions(Nc_init,3),friction0,lf,radius(Nc_init),adh_strength(Nc_init,2),pi,phimm,&
	&friction(3),radius0,radiusagg,rc,Surfaceagg1,Vmax1,phi1m,Vcellsagg1,probacellcr,rcMM,rcm,timeadh,f0,proba_cell,tdelayc,dt,&
	&layer_macr,radiusmacr,surface,radius_macr,phi1m_final,phim,normal_vectorx,normal_vectory,normal_vectorz,friction_normal,&
	&friction_parallel,fminimum,friction_normal_opposite

	pi=3.1415926
	radiusagg=radius(agg1)
		
		if (timeadh.gt.0) then
				f0=proba_cell/timeadh*tdelayc 
		end if
	
		if ((cell_macr(agg1,1).eq.0).and.(cell_macr(agg1,2).gt.0)) then
		phi1m=1
		end if
		
		if ((cell_macr(agg1,1).gt.0).and.(cell_macr(agg1,2).gt.0)) then
			radius01c3=r3(adh_strength(agg1,1)-f0,probacellcr,rcMM,rcm)
			Vcellsagg1=real(cell_macr(agg1,1))*4*pi/3*radius01c3
			Vmax1=4.*pi/3.*((Vcellsagg1/(4*pi/3))**(1./3.)+layer_macr)**3.-Vcellsagg1
			phi1m=real(cell_macr(agg1,2))*(radius_macr**3.)*4*pi/3/(Vmax1)
			radiusagg=((Vcellsagg1+Vmax1)/(4*pi/3))**(1./3.)
		end if
		
		if ((cell_macr(agg1,1).gt.0).and.(cell_macr(agg1,2).eq.0)) then
		phi1m=0
		end if
	rc=10
	radius0=sqrt(positions(agg1,1)**2+positions(agg1,2)**2+positions(agg1,3)**2)+radiusagg
	radius1=sqrt(positions(agg1,1)**2+positions(agg1,2)**2)+radiusagg
	phim=cell_macr(agg1,2)*radius_macr**3/radius(agg1)**3
	phi1m_final=max(phi1m,phim)
	surface=(phi1m_final)**(5./3.)*radius(agg1)**(4./3.)
	
	
	do 6547 i6547=1,3
	friction(i6547)=6*pi*radiusagg*friction0!+&
6547	continue
	
	fminimum=6*pi*radiusagg*friction0
	
	if (positions(agg1,3).lt.0) then
	
	normal_vectorx=positions(agg1,1)/(radius0-radiusagg)
	normal_vectory=positions(agg1,2)/(radius0-radiusagg)
	normal_vectorz=positions(agg1,3)/(radius0-radiusagg)
	
		distance=max((radius_well-radius0),0.5)
		friction_normal=6*pi*friction0*radiusagg*(1+0*radiusagg/distance)!+&
		
		friction_normal_opposite=6*pi*friction0*radiusagg+&
		&friction2*surface*min(exp(-(radius_well-radius0)/lf),1.0) 
		
		friction_parallel=6*pi*friction0*radiusagg*(1-0*8/15*pi*min(log(distance/radiusagg),0.))+&
		&friction2*surface*min(exp(-(radius_well-radius0)/lf),1.0)

	radius0=sqrt(positions(agg1,1)**2+positions(agg1,2)**2+positions(agg1,3)**2)+radiusagg
		
	friction(1)=friction_normal;friction(3)=friction_normal_opposite;friction(2)=friction_parallel
		
	end if
	
return
end
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		SUBROUTINE gravity_force(Nc_init,g,m,gf,mass_kpm)
		REAL, DIMENSION(3) :: gf
		real m,g,mass_c,mass_kpm
		integer nn
		
		gf(1)=0
		gf(2)=0
		gf(3)=-mass_kpm*g
		
		return
		end
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		SUBROUTINE noise_force(Nc_init,nf,amp_noise,idum,positions,friction,dt,rad,agg1,cell_macr,radius_well,amp_noise2,radius,lf&
		&,adh_strength,radius_macr,probacellcr,rcMM,rcm,proba_cell,timeadh,tdelayc,layer_macr,gf,theta_last,amp_noise21_last,&
		&steps_persistence)
		REAL, DIMENSION(3) :: nf
		
		real :: amp_noise,amp_noise0,x1,y1,z1,amp_noise1,radius_well,radius(Nc_init),lf,adh_strength(Nc_init,2)
		integer :: idum,agg1,cell_macr(Nc_init,2),radius0,idum1,j,idum0,iter
		real :: positions(Nc_init,3),rad,theta0,phi0,dx,dy,dz,phim,amp_noise3,amp_noise2,x0,y0,z0,phi1m,radius01c3,Vcellsagg1,&
		&Surfaceagg1,Vmax1,amp_noise11,radius_macr,probacellcr,rcMM,rcm,proba_cell,timeadh,dt,tdelayc,layer_macr,amp_noise22,&
		&amp_noise20,theta,vectx0,vectx1,vecty0,vecty1,vectz0,vectz1,vectx,vecty,vectz,nf1x0,nf1y0,nf1z0,nf1x,nf1y,nf1z,&
		&vectx11,vecty11,vectz11,norm11,surface,phim_v,friction(3),friction_matrix(3,3),friction_table(3,3),fn,fp,enx,eny,enz,&
		&det,forcen,gf(3),forces(3),theta_last,amp_noise21_last,amp_noise21,steps_persistence
		
		character*8 :: date
		character*10 :: time
		character*5 :: zone
		integer,dimension(8) :: values
	
		PI=3.14159265359
		lf=10
	radiusagg=radius(agg1)
	radius0=sqrt(positions(agg1,1)**2+positions(agg1,2)**2+positions(agg1,3)**2)
	radius1=sqrt(positions(agg1,1)**2+positions(agg1,2)**2)
	
	
	if ((cell_macr(agg1,2).eq.1).and.(cell_macr(agg1,1).eq.0)) then
	phim=1
	end if
	
	if (cell_macr(agg1,2).eq.0) then
	phim=0
	end if
	
	f0=0
	
	if ((cell_macr(agg1,1).gt.0).and.(cell_macr(agg1,2).gt.0)) then
		if (timeadh.gt.0) then
			f0=proba_cell/timeadh*tdelayc
		end if
		radius01c3=r3(adh_strength(agg1,1)-f0,probacellcr,rcMM,rcm)
		Vcellsagg1=real(cell_macr(agg1,1))*4*pi/3*radius01c3
		Vmax1=4.*pi/3.*((Vcellsagg1/(4*pi/3))**(1./3.)+layer_macr)**3.-Vcellsagg1
		phim=real(cell_macr(agg1,2))*(radius_macr**3.)*4*pi/3/(Vmax1)
		radiusagg=((Vcellsagg1+Vmax1)/(4*pi/3))**(1./3.)
	end if
	
	!Simple noise
	
		amp_noise11=(amp_noise)
		amp_noise1=amp_noise**2.
		
		
		j=0
		dx=0
		dy=0
		dz=0
		idum0=int((positions(agg1,1)*10)**2+(positions(agg1,2)*10)**2)
		idum1=int((positions(agg1,1)*10)**2+(positions(agg1,2)*10)**2)	
		idum0=agg1+idum1+abs(positions(agg1,1))+abs(positions(agg1,2))+abs(positions(agg1,3))
		iter=0
		idum=idum0
		amp_noise0=sqrt(-2*log(ran2(idum)))*cos(2*PI*ran2(idum*11))*sqrt(amp_noise1)
		idum=int(idum*ran2(idum))
		theta0=ran2(idum)*PI
		idum=int(idum*ran2(idum))
		phi0=ran2(idum)*PI*2
		
963		amp_noise0=sqrt(-2*log(ran2(idum)))*cos(2*PI*ran2(idum*11))*sqrt(amp_noise1)
		x0=positions(agg1,1)
		y0=positions(agg1,2)
		z0=positions(agg1,3)
		r0=sqrt(x0**2+y0**2+z0**2)
		nf1x0=amp_noise0*cos(phi0)*sin(theta0)+dx
		nf1y0=amp_noise0*sin(phi0)*sin(theta0)+dy
		nf1z0=amp_noise0*cos(theta0)+dz
		
		nf(1)=nf1x0
		nf(2)=nf1y0
		nf(3)=nf1z0
		idum=int(idum*ran2(idum))
		
		fn=friction(1) !normal contribution
		fp=friction(2) !parallel contribution
		
	    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		radius0=sqrt(positions(agg1,1)**2+positions(agg1,2)**2+positions(agg1,3)**2)
		
		enx=positions(agg1,1)/radius0;eny=positions(agg1,2)/radius0;enz=positions(agg1,3)/radius0		
		
		forces(1)=nf(1);forces(2)=nf(2);forces(3)=nf(3)!+gf(3)
		
		if (radius_well-(radius0+radius(agg1)).gt.50) then
			forces(3)=nf(3)
		end if
		
		forcen=forces(1)*enx+forces(2)*eny+forces(3)*enz !normal force towards the substrate or not.
		
		
		if (forcen.lt.0) then
			fn=friction(3)
		end if
		
friction_matrix(1,1)=fp+enx*enx*(fn-fp);friction_matrix(1,2)=   eny*enx*(fn-fp);friction_matrix(1,3)=   enz*enx*(fn-fp);
friction_matrix(2,1)=   enx*eny*(fn-fp);friction_matrix(2,2)=fp+eny*eny*(fn-fp);friction_matrix(2,3)=   enz*eny*(fn-fp);		
friction_matrix(3,1)=   enx*enz*(fn-fp);friction_matrix(3,2)=   eny*enz*(fn-fp);friction_matrix(3,3)=fp+enz*enz*(fn-fp);

det=friction_matrix(1,1)*friction_matrix(2,2)*friction_matrix(3,3)-friction_matrix(1,1)*friction_matrix(2,3)**2&
&-friction_matrix(2,2)*friction_matrix(1,3)**2-friction_matrix(3,3)*friction_matrix(1,2)**2&
&+2*friction_matrix(1,3)*friction_matrix(2,3)*friction_matrix(1,2)
			
friction_table(1,1)=friction_matrix(2,2)*friction_matrix(3,3)-friction_matrix(2,3)**2
friction_table(2,2)=friction_matrix(1,1)*friction_matrix(3,3)-friction_matrix(1,3)**2
friction_table(3,3)=friction_matrix(1,1)*friction_matrix(2,2)-friction_matrix(1,2)**2

friction_table(1,2)=friction_matrix(1,3)*friction_matrix(2,3)-friction_matrix(1,2)*friction_matrix(3,3)
friction_table(2,1)=friction_table(1,2)

friction_table(1,3)=friction_matrix(1,2)*friction_matrix(2,3)-friction_matrix(1,3)*friction_matrix(2,2)
friction_table(3,1)=friction_table(1,3)

friction_table(2,3)=friction_matrix(1,2)*friction_matrix(1,3)-friction_matrix(2,3)*friction_matrix(1,1)
friction_table(3,2)=friction_table(2,3)
		
		
		
		x1=positions(agg1,1)+((nf(1))*friction_table(1,1)+(nf(2))*friction_table(1,2)+&
		&(forces(3))*friction_table(1,3))*dt/det
		y1=positions(agg1,2)+((nf(1))*friction_table(2,1)+(nf(2))*friction_table(2,2)+&
		&(forces(3))*friction_table(2,3))*dt/det
		z1=positions(agg1,3)+((nf(1))*friction_table(3,1)+(nf(2))*friction_table(3,2)+&
		&(forces(3))*friction_table(3,3))*dt/det
		
		
	radius00=sqrt(positions(agg1,1)**2+positions(agg1,2)**2+positions(agg1,3)**2)
	radius11=sqrt(positions(agg1,1)**2+positions(agg1,2)**2)
	
			if ((positions(agg1,3).lt.0).and.((radius00+radius(agg1)).gt.radius_well)) then
				print*, 'problem initial position noise down 0',radius00+radius(agg1),radius(agg1)
				stop
			elseif ((positions(agg1,3).ge.0).and.((radius11+radius(agg1)).gt.radius_well)) then
				print*, 'problem initial position noise up 0',radius11+radius(agg1),radius(agg1)
				stop
			end if


		idum=idum+1000
		
			
		if (phim.gt.0) then
			
			
		idum=idum+1000
		phim_v=cell_macr(agg1,2)*radius_macr**3/radius(agg1)**3
		phim_final=max(phim,phim_v)
		surface=min(phim_final**(5./3.)*radius(agg1)**(4./3.)/radius_macr**(4./3.),1.)
		amp_noise22=amp_noise2*surface*min(exp(-(radius_well-radius(agg1)-radius0)/lf),1.0)!*sqrt(phim**(5./3.)*radius(agg1)**(4./3.))
		amp_noise21=sqrt(-2*log(ran2(idum)))*cos(2*PI*ran2(idum+1000))
		amp_noise20=amp_noise21*abs(amp_noise22)
			!direction perpendiculaire au rayon
			
		if ((z0.lt.-1e-5).and.(radius_well-(r0+radiusagg).lt.lf)) then !https://math.stackexchange.com/questions/2347215/how-to-find-a-random-unit-vector-orthogonal-to-a-random-unit-vector-in-3d
			idum=idum+1000
			vectx0=y0/sqrt(y0**2+x0**2)
			vecty0=-x0/sqrt(y0**2+x0**2)
			vectz0=0
			vectx11=vecty0*z0
			vecty11=-vectx0*z0
			vectz11=vectx0*y0-vecty0*x0
			norm11=sqrt(vectx11**2+vecty11**2+vectz11**2)
			vectx1=vectx11/norm11
			vecty1=vecty11/norm11
			vectz1=vectz11/norm11
						
			if (theta_last.lt.1e4) then
			theta=theta_last
			amp_noise20=amp_noise21_last*abs(amp_noise22)
			elseif (theta_last.gt.1e4) then
			idum=idum+100
			theta=ran2(idum)*2*pi
			amp_noise20=amp_noise21*abs(amp_noise22)
			theta_last=theta
			amp_noise21_last=amp_noise21
			end if

			vectx=cos(theta)*vectx0+sin(theta)*vectx1
			vecty=cos(theta)*vecty0+sin(theta)*vecty1
			vectz=cos(theta)*vectz0+sin(theta)*vectz1
			nf1x=vectx*amp_noise20
			nf1y=vecty*amp_noise20
			nf1z=vectz*amp_noise20
			nf(1)=nf1x+nf1x0
			nf(2)=nf1y+nf1y0
			nf(3)=nf1z+nf1z0
			
		end if
		end if
		

		
		return
		end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		subroutine new_velocity_force(agg1,nf,gf,normf,af,new_v,Nc_init,friction,friction0,friction1,friction2,dt,&
		&cell_macr,radius_well,positions,lf,radius)
		
		
		real, DIMENSION(3) :: new_v
		real, DIMENSION(3) :: gf,nf,af,normf
		REAL, DIMENSION(Nc_init,3) :: positions
		integer :: agg1,Nc_init,cell_macr(Nc_init,2),j
		real :: displacement,dt,friction0,friction1,friction2,friction(3),radius_well,lf,friction_table(3,3)&
		&,radius(Nc_init),friction_matrix(3,3),fn,fp,det,forces(3),newf_pos(3),frictionp(3),forcen,enx,eny,enz
		
		fn=friction(1) !normal contribution
		fp=friction(2) !parallel contribution
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		forces(1)=nf(1)+gf(1)
		forces(2)=nf(2)+gf(2)
		forces(3)=nf(3)+gf(3)
		
		radius0=sqrt(positions(agg1,1)**2+positions(agg1,2)**2+positions(agg1,3)**2)
		
		enx=positions(agg1,1)/radius0;eny=positions(agg1,2)/radius0;enz=positions(agg1,3)/radius0
		
		forcen=forces(1)*enx+forces(2)*eny+forces(3)*enz !force normal to the substrate
		
		
		
		if (forcen.lt.0) then 
			fn=friction(3)
		end if
		
		
friction_matrix(1,1)=fp+enx*enx*(fn-fp);friction_matrix(1,2)=   eny*enx*(fn-fp);friction_matrix(1,3)=   enz*enx*(fn-fp);
friction_matrix(2,1)=   enx*eny*(fn-fp);friction_matrix(2,2)=fp+eny*eny*(fn-fp);friction_matrix(2,3)=   enz*eny*(fn-fp);		
friction_matrix(3,1)=   enx*enz*(fn-fp);friction_matrix(3,2)=   eny*enz*(fn-fp);friction_matrix(3,3)=fp+enz*enz*(fn-fp);


det=friction_matrix(1,1)*friction_matrix(2,2)*friction_matrix(3,3)-friction_matrix(1,1)*friction_matrix(2,3)**2&
&-friction_matrix(2,2)*friction_matrix(1,3)**2-friction_matrix(3,3)*friction_matrix(1,2)**2&
&+2*friction_matrix(1,3)*friction_matrix(2,3)*friction_matrix(1,2)
			
friction_table(1,1)=friction_matrix(2,2)*friction_matrix(3,3)-friction_matrix(2,3)**2
friction_table(2,2)=friction_matrix(1,1)*friction_matrix(3,3)-friction_matrix(1,3)**2
friction_table(3,3)=friction_matrix(1,1)*friction_matrix(2,2)-friction_matrix(1,2)**2

friction_table(1,2)=friction_matrix(1,3)*friction_matrix(2,3)-friction_matrix(1,2)*friction_matrix(3,3)
friction_table(2,1)=friction_table(1,2)

friction_table(1,3)=friction_matrix(1,2)*friction_matrix(2,3)-friction_matrix(1,3)*friction_matrix(2,2)
friction_table(3,1)=friction_table(1,3)

friction_table(2,3)=friction_matrix(1,2)*friction_matrix(1,3)-friction_matrix(2,3)*friction_matrix(1,1)
friction_table(3,2)=friction_table(2,3)

		do 4000 i4000=1,3

			new_v(i4000)=(friction_table(i4000,1)*forces(1)+friction_table(i4000,2)*forces(2)&
			&+friction_table(i4000,3)*forces(3))/det
	
			displacement=abs(dt*new_v(i4000))
						
4000	continue

		return
		end
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		subroutine new_position_force(agg1,positions,newf_pos,new_velocity,Nc_init,dt,nf,gf,normf,sum_adh_force,friction,radius,&
		&radius_well)
		
		real, DIMENSION(3) :: gf,nf,normf,sum_adh_force
		real, DIMENSION(3) :: newf_pos
		real, DIMENSION(3) :: new_velocity
		real :: dt,displacement,friction(2)
		real :: radius(Nc_init),radius_well,radius0,radius1,rapport0,rapport1
		REAL, DIMENSION(Nc_init,3)  ::	positions
		integer :: agg1
		do 5000 i5000=1,3
		    displacement=dt*new_velocity(i5000)
			newf_pos(i5000)=positions(agg1,i5000)+displacement
			
			
			
			if (newf_pos(3).gt.7e3) then
			newf_pos(3)=7e3-0.1
				
				end if
			
5000	continue

		radius0=sqrt(newf_pos(1)**2+newf_pos(2)**2+newf_pos(3)**2)+radius(agg1)
		rapport0=(radius_well-(radius(agg1)+2))/sqrt(newf_pos(1)**2+newf_pos(2)**2+newf_pos(3)**2)
		radius1=sqrt(newf_pos(1)**2+newf_pos(2)**2)+radius(agg1)
		rapport1=(radius_well-(radius(agg1)+2))/sqrt(newf_pos(1)**2+newf_pos(2)**2)
		if (newf_pos(3).gt.7e3) then
			newf_pos(3)=6e3
		elseif ((newf_pos(3).ge.0).and.(newf_pos(3).lt.7e3).and.(radius1.gt.radius_well)) then
			newf_pos(1)=newf_pos(1)*rapport1
			newf_pos(2)=newf_pos(2)*rapport1
		elseif ((newf_pos(3).lt.0).and.(radius0.gt.radius_well)) then
			newf_pos(1)=newf_pos(1)*rapport0
			newf_pos(2)=newf_pos(2)*rapport0
			newf_pos(3)=newf_pos(3)*rapport0
		end if
		
		radius0=sqrt(newf_pos(1)**2+newf_pos(2)**2+newf_pos(3)**2)+radius(agg1)
		radius1=sqrt(newf_pos(1)**2+newf_pos(2)**2)+radius(agg1)
		return
		end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	subroutine fusion_function(agg1,agg2,mass,radius,positions,cell_macr,Nc_init,num_agg,ntotal,proba_cell,proba_macr&
	&,nc_cells,nc_macr,distancefusion,adh_strength,dt,timeadh,timeadhm,tdelayc,tdelaym,rcMM,rcm,rmMM,rmm,rep_macr,&
	&probacellcr,B,tauagg,radius_cell,radius_macr,phimm,mass_c,mass_m,layer_macr,theta_table)
	
	integer :: agg1,agg2,idum,nc_cells,nc_macr,n,m0,m00,n00
	real :: positions(1:Nc_init,1:3),mass(1:Nc_init),radius(1:Nc_init),rep_macr,distance0,B,&
	&mass0,radius0,positions0(3),distance,distancefusion,adh_strength(1:Nc_init,1:2),dt,timeadh,timeadhm
	integer :: cell_macr0(2),cell_macr(1:Nc_init,1:2),ntotal(1:Nc_init),m,mext
	real :: proportion_cell,proba_adhesion1,proba_adhesion2,proba_adhesion,ran2i,adh_strength0,&
	&proba_cell,proba_macr,adhs2,adhs1,npc1,npc2,npc22,npc11,rc,adh_strength1,proba_cellm,phimm,&
	&f0,f1,adhs11,adhs12,adhs21,adhs22,tdelayc,tdelaym,effectiveradius,r12m3,probacellcr,pi,&
	&deltaEjkr,deltaEjkrsubstrat,Vag,Vagg1,Vagg2,fag,fagg1,fagg2,FFag,Fdesag,deltaF,difference,tauagg,&
	&phi1c,Vcellsagg1,fcellsagg1,Surfaceagg1,Vmax1,phi1m,fsurfaceagg1,Eagg1,&
	&phi2c,Vcellsagg2,fcellsagg2,Surfaceagg2,Vmax2,phi2m,fsurfaceagg2,Eagg2,&
	&phi3c,Vcellsagg3,fcellsagg3,Surfaceagg3,Vmax3,phi3m,fsurfaceagg3,Eagg3,radius_macr,radius_cell,&
	&mmax_int,mmax_real,mass_c,mass_m,alpha,beta,s,s0,beta0,energy_barrier,layer_macr,theta_table(Nc_init,3)

	integer :: new_agg,new_zero,num_agg(1:Nc_init)
	character*8 :: date
	character*10 :: time
	character*1000 :: message_erreur
	character*5 :: zone
	integer,dimension(8) :: values
	radius_well=3000.
	
	if ((cell_macr(agg1,1)+cell_macr(agg1,2).ge.1).and.(cell_macr(agg2,1)+cell_macr(agg2,2).ge.1)) then
	
	distance=sqrt((positions(agg1,1)-positions(agg2,1))**2.+(positions(agg1,2)-positions(agg2,2))**2+&
	&(positions(agg1,3)-positions(agg2,3))**2)-radius(agg1)-radius(agg2)&
	&-distancefusion
	rc=10
	
	adh_strength0=0
	adh_strength1=0
	
	!!!!!!!!!!!!!!!!!!!!!!
	proba_adhesion=0
	pi=3.14159265359
	if (distance.lt.0) then
		do 7003 i7003=1,3
				positions0(i7003)=0
	 7003	continue
		do 7004 i7004=1,2
			cell_macr0(i7004)=0
	 7004	continue
		mass0=0
		radius0=0
		f0=0
		f1=0
			if (timeadh.gt.0) then
			f0=proba_cell/timeadh*tdelayc
			end if
			if (timeadhm.gt.0) then
			f1=proba_macr/timeadhm*tdelaym
			end if
			
		adhs11=adh_strength(agg1,1)-f0
		adhs21=adh_strength(agg1,2)-f1
		
		adhs12=adh_strength(agg2,1)-f0
		adhs22=adh_strength(agg2,2)-f1
		
		effectiveradius=1/(1/radius(agg1)+1/radius(agg2))
		
		if (cell_macr(agg1,1)+cell_macr(agg2,1).gt.0) then
					adh_strength0=(real(cell_macr(agg1,1))*adh_strength(agg1,1)+&
					&real(cell_macr(agg2,1))*adh_strength(agg2,1))/real(cell_macr(agg1,1)+cell_macr(agg2,1))
		end if
		if (cell_macr(agg1,2)+cell_macr(agg2,2).gt.0) then
					adh_strength1=(real(cell_macr(agg1,2))*adh_strength(agg1,2)+&
					&real(cell_macr(agg2,2))*adh_strength(agg2,2))/real(cell_macr(agg1,2)+cell_macr(agg2,2))
		end if
		
		
			radius01c3=0;radius02c3=0;radius01m3=0;radius02m3=0;
			
			radius01c3=r3(adh_strength(agg1,1)-f0,probacellcr,rcMM,rcm)
			radius02c3=r3(adh_strength(agg2,1)-f0,probacellcr,rcMM,rcm)
			radius01m3=rmm**3
			radius02m3=rmm**3
			
			radius012c3=r3(adh_strength0-f0,probacellcr,rcMM,rcm)
			radius012m3=rmm**3
			
				if (cell_macr(agg1,1)+cell_macr(agg2,1).gt.0) then
				r12c3=real(cell_macr(agg1,1))*radius01c3+real(cell_macr(agg2,1))*radius02c3
				end if
				
				r12m3=0
				if (cell_macr(agg1,2)+cell_macr(agg2,2).gt.0) then
				r12m3=real(cell_macr(agg1,2))*radius01m3+real(cell_macr(agg2,2))*radius02m3
				end if
				
				
			radius0=(r12c3+r12m3)**(1./3.)
		
		!!!!!!!!!!!!!!!!!!!!!!
		
		idum=agg1+agg2+cell_macr(agg1,1)+cell_macr(agg1,2)+cell_macr(agg2,1)+cell_macr(agg2,2)+&
		&int(positions(agg1,1)+positions(agg1,2)+positions(agg1,3))*10
		
			posi0=((radius(agg1)**3)*positions(agg1,1)+&
			&(radius(agg2)**3)*positions(agg2,1))/(radius(agg1)**3+radius(agg2)**3)
			posi1=((radius(agg1)**3)*positions(agg1,2)+&
			&(radius(agg2)**3)*positions(agg2,2))/(radius(agg1)**3+radius(agg2)**3)
			posi2=((radius(agg1)**3)*positions(agg1,3)+&
			&(radius(agg2)**3)*positions(agg2,3))/(radius(agg1)**3+radius(agg2)**3)
			
			pp1=1.
			pp2=1.
			pp=1.
			
			if (positions(agg1,3).gt.0) then
				pp1=0.
			end if
			if (positions(agg2,3).gt.0) then
				pp2=0.
			end if
			if (posi.gt.0) then
				pp=0.
			end if
			
		rad=sqrt(posi0**2+posi1**2+pp*posi2**2)+radius0
		rad1=sqrt(positions(agg1,1)**2+positions(agg1,2)**2+pp1*positions(agg1,3)**2)+radius(agg1)
		rad2=sqrt(positions(agg2,1)**2+positions(agg2,2)**2+pp2*positions(agg2,3)**2)+radius(agg2)
		
		p1=min(1.,exp((rad1-radius_well)/(2*rc)))
		p2=min(1.,exp((rad2-radius_well)/(2*rc)))
		p3=min(1.,exp((rad-radius_well)/(2*rc)))
		
		ran2i=ran2(idum)
			if (cell_macr(agg1,1)+cell_macr(agg2,1).gt.0) then
					adh_strength0=(cell_macr(agg1,1)*adh_strength(agg1,1)+&
					&cell_macr(agg2,1)*adh_strength(agg2,1))/(cell_macr(agg1,1)+cell_macr(agg2,1))
			end if
			if (cell_macr(agg1,2)+cell_macr(agg2,2).gt.0) then
					adh_strength1=(cell_macr(agg1,2)*adh_strength(agg1,2)+&
					&cell_macr(agg2,2)*adh_strength(agg2,2))/(cell_macr(agg1,2)+cell_macr(agg2,2))
			end if
					
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		phi1c=radius_cell**3/radius01c3
		Vcellsagg1=real(cell_macr(agg1,1))*4*pi/3*radius01c3
		Vmacragg2=real(cell_macr(agg2,2))*4*pi/3*radius02m3
		l1barrier=0
			if ((cell_macr(agg1,1).gt.0).and.(cell_macr(agg1,2).gt.0)) then
			l1barrier=(real(cell_macr(agg1,1))*radius01c3+real(cell_macr(agg1,2))*radius01m3)**(1./3.)-&
			&(real(cell_macr(agg1,1))*radius01c3)**(1./3.)
			end if
		fcellsagg1=0
			if (cell_macr(agg1,1).gt.1) then 
			fcellsagg1=-phi1c**2*max(adhs11,0.)
			end if
		Surfaceagg1=4.*pi*(Vcellsagg1/(4*pi/3))**(2./3.)
		Vmax1=4.*pi/3.*((Vcellsagg1/(4*pi/3))**(1./3.)+layer_macr)**3.-Vcellsagg1
		phi1m=0
		phi1m=real(cell_macr(agg1,2))*(radius_macr**3.)*4*pi/3/(Vmax1)
		fsurfaceagg1=-phi1c*phi1m*max(adhs21,0.)
		Eagg1=Vcellsagg1*fcellsagg1+Surfaceagg1*fsurfaceagg1
		
		phi2c=radius_cell**3/radius02c3
		Vcellsagg2=real(cell_macr(agg2,1))*4*pi/3*radius02c3
		Vmacragg2=real(cell_macr(agg2,2))*4*pi/3*radius02m3
		l2barrier=0
			if ((cell_macr(agg2,1).gt.0).and.(cell_macr(agg2,2).gt.0)) then
			l2barrier=(real(cell_macr(agg2,1))*radius02c3+real(cell_macr(agg2,2))*radius02m3)**(1./3.)-&
			&(real(cell_macr(agg2,1))*radius02c3)**(1./3.)
			end if
		fcellsagg2=0
			if (cell_macr(agg2,1).gt.1) then 
			fcellsagg2=-phi2c**2*max(adhs12,0.)
			end if
		Surfaceagg2=4.*pi*(Vcellsagg2/(4*pi/3))**(2./3.)
		Vmax2=4.*pi/3.*((Vcellsagg2/(4*pi/3))**(1./3.)+layer_macr)**3.-Vcellsagg2
		phi2m=0
		phi2m=real(cell_macr(agg2,2))*(radius_macr**3.)*4*pi/3/(Vmax2)
		fsurfaceagg2=-phi2c*phi2m*max(adhs22,0.)
		Eagg2=Vcellsagg2*fcellsagg2+Surfaceagg2*fsurfaceagg2
		
		phi3c=radius_cell**3/radius012c3
		Vcellsagg3=real(cell_macr(agg1,1)+cell_macr(agg2,1))*4*pi/3*radius012c3
		fcellsagg3=0
			if ((cell_macr(agg1,1)+cell_macr(agg2,1)).gt.1) then 
			fcellsagg3=-phi3c**2*max(adh_strength0-f0,0.)
			end if
		Surfaceagg3=4.*pi*(Vcellsagg3/(4*pi/3))**(2./3.)
		Vmax3=4.*pi/3.*((Vcellsagg3/(4*pi/3))**(1./3.)+layer_macr)**3.-Vcellsagg3
		phi3m=real(cell_macr(agg1,2)+cell_macr(agg2,2))*(radius_macr**3.)*4*pi/3/(Vmax3)
		fsurfaceagg3=-phi3c*max(phi3m,1.)*max(adh_strength1-f1,0.)
		Eagg3=Vcellsagg3*fcellsagg3+Surfaceagg3*fsurfaceagg3
		
		deltaF=Eagg3-(Eagg1+Eagg2)
		energy_barrier=0
		
		energy_barrier=10*phi1m*phi2m*rep_macr
		energy_barrier_anticd11c=(phi1m*(1-phi2m)+(1-phi1m)*phi2m)*B
		
		barrier=energy_barrier
		proba_adhesion=min(exp(-barrier),1.)*min(exp(-energy_barrier_anticd11c),1.)/tauagg*dt
		
			if ((cell_macr(agg1,2)+cell_macr(agg2,2).eq.0).and.((adhs11.le.0).or.(adhs12.le.0))) then
			proba_adhesion=0
			elseif (cell_macr(agg1,1)+cell_macr(agg2,1).eq.0) then
			proba_adhesion=0
			elseif ((cell_macr(agg1,2).eq.1).and.(cell_macr(agg2,1).eq.1).and.(adhs21.le.0).and.(adhs12.le.0)) then
			proba_adhesion=0
			elseif ((cell_macr(agg2,2).eq.1).and.(cell_macr(agg1,1).eq.1).and.(adhs22.le.0).and.(adhs11.le.0)) then
			proba_adhesion=0
			end if
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!!! if fusion
		
		if (ran2i.lt.proba_adhesion) then
			mass0=mass(agg1)+mass(agg2)
			
			 
			new_agg=agg1
			new_zero=agg2
				if (agg1.gt.agg2) then
				new_agg=agg2
				new_zero=agg1
				end if 
			

			adh_strength(new_agg,1)=adh_strength0
			adh_strength(new_agg,2)=adh_strength1

			cell_macr(new_agg,2)=cell_macr(agg1,2)+cell_macr(agg2,2)
			cell_macr(new_agg,1)=cell_macr(agg1,1)+cell_macr(agg2,1)
			
			cell_macr(new_zero,1)=0
			cell_macr(new_zero,2)=0

			ntotal(new_agg)=cell_macr(new_agg,1)+cell_macr(new_agg,2)
			ntotal(new_zero)=0
				do 7002 i7002=1,3
				posi=((radius(agg1)**3)*positions(agg1,i7002)+&
				&(radius(agg2)**3)*positions(agg2,i7002))/(radius(agg1)**3+radius(agg2)**3)
				
				
				positions(new_agg,i7002)=posi
				7002	continue	

			mass(new_agg)=mass0
			mass(new_zero)=0
			 
			radius03=r3(adh_strength(new_agg,1)-f0,probacellcr,rcMM,rcm)
			
			rc3=cell_macr(new_agg,1)*radius03+cell_macr(new_agg,2)*rmm**3
			
			radius(new_agg)=rc3**(1./3.)
			radius(new_zero)=0
		
			theta_table(new_agg,1)=1.
			theta_table(new_zero,1)=1.
		
	radius_well=3000
	
		if ((cell_macr(new_agg,1).gt.0).and.(cell_macr(new_agg,2).gt.0)) then
			call expulsion_macr(new_agg,cell_macr,phimm,ntotal,radius,adh_strength,mass,mass_c,mass_m,radius_macr,proba_macr,positions&
			&,Nc_init,radius012c3,layer_macr,radius_well)
				end if	

			radius0=sqrt(positions(new_agg,1)**2+positions(new_agg,2)**2+positions(new_agg,3)**2)
			radius1=sqrt(positions(new_agg,1)**2+positions(new_agg,2)**2)
	
					if ((positions(new_agg,3).lt.0).and.((radius0+radius(new_agg)).gt.radius_well)) then
						pos1=positions(new_agg,1)*(radius_well-1-radius(new_agg))/(radius0)
						pos2=positions(new_agg,2)*(radius_well-1-radius(new_agg))/(radius0)
						pos3=positions(new_agg,3)*(radius_well-1-radius(new_agg))/(radius0)
						positions(new_agg,1)=pos1
						positions(new_agg,2)=pos2
						positions(new_agg,3)=pos3
					elseif ((positions(new_agg,3).gt.0).and.((radius1+radius(new_agg)).gt.radius_well)) then
						pos1=positions(new_agg,1)*(radius_well-1-radius(new_agg))/(radius1)
						pos2=positions(new_agg,2)*(radius_well-1-radius(new_agg))/(radius1)
						positions(new_agg,1)=pos1
						positions(new_agg,2)=pos2
					end if	
		end if
			
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!!!!!!! if no fusion
		
		
		if (ran2i.gt.proba_adhesion) then
	
			if ((cell_macr(agg1,1).eq.1).and.(adh_strength(agg1,1).lt.f0)) then 
				adhs1=adh_strength(agg1,1)
				npc1=adhs1+proba_cell/timeadh*dt
				
				npc11=min(proba_cell+f0,npc1)
				adh_strength(agg1,1)=npc11
			end if 
			
			if ((cell_macr(agg2,1).eq.1).and.(adh_strength(agg2,1).lt.f0)) then 
				adhs2=adh_strength(agg2,1)
				npc2=adhs2+proba_cell/timeadh*dt
				npc22=min(proba_cell+f0,npc2)
				adh_strength(agg2,1)=npc22
				
			end if
			
			if ((cell_macr(agg1,2).eq.1).and.(adh_strength(agg1,2).lt.f1)) then 
				npc1=adh_strength(agg1,2)
					if (timeadhm.gt.0) then 
					npc1=adh_strength(agg1,2)+proba_macr/timeadhm*dt
					end if 
				npc11=min(proba_macr+f1,npc1)
				adh_strength(agg1,2)=npc11
			end if 
			
			if ((cell_macr(agg2,2).eq.1).and.(adh_strength(agg2,2).lt.f1)) then 
				npc2=adh_strength(agg2,2)
					if (timeadhm.gt.0) then 
					npc2=adh_strength(agg2,2)+proba_macr/timeadhm*dt
					end if
				npc22=min(proba_macr+f1,npc2)
				adh_strength(agg2,2)=npc22
			end if
			
			
			radius01c3=r3(adh_strength(agg1,1)-f0,probacellcr,rcMM,rcm);radius02c3=r3(adh_strength(agg2,1)-f0,probacellcr,rcMM,rcm)
			radius01m3=rmm**3;radius02m3=rmm**3
			
			radius(agg1)=(real(cell_macr(agg1,1))*radius01c3+real(cell_macr(agg1,2))*radius01m3)**(1./3.)
			radius(agg2)=(real(cell_macr(agg2,1))*radius02c3+real(cell_macr(agg2,2))*radius02m3)**(1./3.)
			
			if ((cell_macr(agg1,1).gt.0).and.(cell_macr(agg1,2).gt.0)) then
			call expulsion_macr(agg1,cell_macr,phimm,ntotal,radius,adh_strength,mass,mass_c,mass_m,radius_macr,&
			&proba_macr,positions,Nc_init,radius01c3,layer_macr,radius_well)
			end if	
			if ((cell_macr(agg2,1).gt.0).and.(cell_macr(agg2,2).gt.0)) then
			call expulsion_macr(agg2,cell_macr,phimm,ntotal,radius,adh_strength,mass,mass_c,mass_m,radius_macr,&
			&proba_macr,positions,Nc_init,radius02c3,layer_macr,radius_well)
			end if	
		
		if ((radius(agg1).ne.radius(agg1)).or.(radius(agg2).ne.radius(agg2)).or.(cell_macr(agg1,1).lt.0).or.(cell_macr(agg2,1).lt.0)&
			&.or.(cell_macr(agg1,2).lt.0).or.(cell_macr(agg2,2).lt.0).or.(radius(agg1).gt.1e4).or.(radius(agg2).gt.1e4)&
			&.or.(adh_strength(agg1,1).ne.adh_strength(agg1,1)).or.(adh_strength(agg1,2).ne.adh_strength(agg1,2))&
			&.or.(adh_strength(agg2,1).ne.adh_strength(agg2,1)).or.(adh_strength(agg2,2).ne.adh_strength(agg2,2))) then
		
			print*, 'problem non fusion !','rad',radius(agg1),radius(agg2),'N',cell_macr(agg1,1),cell_macr(agg2,1),cell_macr(agg1,2)&
			&,cell_macr(agg2,2),&
			&'r12c/m3',r12c3,r12m3,'radiuscell',radius01c3,radius02c3,radius01m3,radius02m3,&
			&'adh',adh_strength(agg1,1),adh_strength(agg1,2),adh_strength(agg2,1),adh_strength(agg2,2),&
			&'fusion ?',ran2i,proba_adhesion
			stop
			end if
	
			distance0=sqrt((positions(agg1,1)-positions(agg2,1))**2+(positions(agg1,2)-positions(agg2,2))**2&
			&+(positions(agg1,3)-positions(agg2,3))**2)-radius(agg1)-radius(agg2)
			
			if (distance0.lt.0) then
				radius_well=3000.
				
		call position_function_new(agg1,agg2,cell_macr(agg1,1),cell_macr(agg2,1),cell_macr(agg1,2),cell_macr(agg2,2)&
		&,positions(agg1,1),positions(agg1,2),positions(agg1,3),positions(agg2,1),positions(agg2,2),positions(agg2,3),&
		&radius_well,radius(agg1),radius(agg2))
				
			end if
			
			new_agg=agg1
			radius0=sqrt(positions(new_agg,1)**2+positions(new_agg,2)**2+positions(new_agg,3)**2)
			radius1=sqrt(positions(new_agg,1)**2+positions(new_agg,2)**2)
					if ((positions(new_agg,3).lt.0).and.((radius0+radius(new_agg)).gt.radius_well)) then
						pos1=positions(new_agg,1)*(radius_well-1-radius(new_agg))/(radius0)
						pos2=positions(new_agg,2)*(radius_well-1-radius(new_agg))/(radius0)
						pos3=positions(new_agg,3)*(radius_well-1-radius(new_agg))/(radius0)
						positions(new_agg,1)=pos1
						positions(new_agg,2)=pos2
						positions(new_agg,3)=pos3
					elseif ((positions(new_agg,3).gt.0).and.((radius1+radius(new_agg)).gt.radius_well)) then
						pos1=positions(new_agg,1)*(radius_well-1-radius(new_agg))/(radius1)
						pos2=positions(new_agg,2)*(radius_well-1-radius(new_agg))/(radius1)
						positions(new_agg,1)=pos1
						positions(new_agg,2)=pos2
					end if	
			new_agg=agg2
			radius0=sqrt(positions(new_agg,1)**2+positions(new_agg,2)**2+positions(new_agg,3)**2)
			radius1=sqrt(positions(new_agg,1)**2+positions(new_agg,2)**2)
	
					if ((positions(new_agg,3).lt.0).and.((radius0+radius(new_agg)).gt.radius_well)) then
						pos1=positions(new_agg,1)*(radius_well-1-radius(new_agg))/(radius0)
						pos2=positions(new_agg,2)*(radius_well-1-radius(new_agg))/(radius0)
						pos3=positions(new_agg,3)*(radius_well-1-radius(new_agg))/(radius0)
						positions(new_agg,1)=pos1
						positions(new_agg,2)=pos2
						positions(new_agg,3)=pos3
					elseif ((positions(new_agg,3).gt.0).and.((radius1+radius(new_agg)).gt.radius_well)) then
						pos1=positions(new_agg,1)*(radius_well-1-radius(new_agg))/(radius1)
						pos2=positions(new_agg,2)*(radius_well-1-radius(new_agg))/(radius1)
						positions(new_agg,1)=pos1
						positions(new_agg,2)=pos2
					end if	
					
	
		end if
	end if
	end if
	
	return
	
	end
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine fragmentation(agg1,mass,radius,positions,cell_macr,Nc_init,num_agg,ntotal,proba_cell,proba_macr&
	&,nc_cells,nc_macr,distancefusion,adh_strength,dt,timeadh,timeadhm,tdelayc,tdelaym,rcMM,rcm,rmMM,rmm,rep_macr,mass_c,mass_m,&
	&probacellcr,B,tauagg,rKP)
	
	integer :: agg1,agg2,idum,nc_cells,nc_macr
	real :: positions(1:Nc_init,1:3),mass(1:Nc_init),radius(1:Nc_init),rep_macr,distance0,probacellcr,B,&
	&mass0,radius0,positions0(3),distance,distancefusion,adh_strength(1:Nc_init,1:2),dt,timeadh,timeadhm
	integer :: cell_macr0(2),cell_macr(1:Nc_init,1:2),ntotal(1:Nc_init),cell_macr_intermediaire,m,IR,rnbin
	real :: proportion_cell,proba_adhesion1,proba_adhesion2,proba_adhesion,ran2i,adh_strength0,lambda,proportion_cell1&
	&proba_cell,proba_macr,adhs2,adhs1,npc1,npc2,npc22,npc11,rc,adh_strength1,mass_c,mass_m,radius01m3,normalization
	real :: f0,f1,adhs11,adhs12,adhs21,adhs22,tdelayc,tdelaym,effectiveradius,r12m3,DeltaG,radius01c3,deltaG0,&
	&proportion_celli,proportion_cellf,Vi,Vf,deltaEjkr,deltaEjkrsubstrat,Fdesag,Fag,deltaF,difference,proba_fragmentation,&
	&energy_barrier,p0,radius1,cc,energy_barrier0,mm
	integer :: new_agg,new_zero,num_agg(1:Nc_init),n,m0,c,c0,c00,m00
	character*8 :: date
	character*10 :: time
	character*5 :: zone
	integer,dimension(8) :: values
	
	rc=10
	f0=0
	
	if ((cell_macr(agg1,2).ge.1).or.(cell_macr(agg1,1).gt.1)) then
		
			if (timeadh.gt.0) then
			f0=proba_cell/timeadh*tdelayc
			end if
	
		adhs11=adh_strength(agg1,1)-f0
		adhs12=adh_strength(agg1,2)
	
		radius01c3=r3(adhs11,probacellcr,rcMM,rcm)
		radius01m3=rmm**3
		radius0=(real(cell_macr(agg1,1))*radius01c3+real(cell_macr(agg1,2)-1)*radius01m3)**(1./3.)
		radius1=(real(cell_macr(agg1,1)-1)*radius01c3+real(cell_macr(agg1,2))*radius01m3)**(1./3.)
	        
			p=ran2(agg1+cell_macr(agg1,2)+cell_macr(agg1,1)+int(radius0)*5)
			p0=ran2(agg1*2+cell_macr(agg1,2)+cell_macr(agg1,1)+int(radius1)*5)
			
			
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
		m=0
		if ((cell_macr(agg1,2).ge.1).and.(cell_macr(agg1,1).ge.1)) then
		effectiveradius=1/(1/radius0+1/rmMM)
		proportion_celli=(real(cell_macr(agg1,1))*radius01c3)/(real(cell_macr(agg1,1))*radius01c3+real(cell_macr(agg1,2))*radius01m3)
		proportion_macri=(real(cell_macr(agg1,2)-1)*radius01m3)/(real(cell_macr(agg1,1))*radius01c3+real(cell_macr(agg1,2))*radius01m3)		
		proportion_cellf=(real(cell_macr(agg1,1))*radius01c3)/(real(cell_macr(agg1,1))*radius01c3+real(cell_macr(agg1,2)-1)*radius01m3)	
		proportion_macrf=(real(cell_macr(agg1,2)-1)*radius01m3)/(real(cell_macr(agg1,1))*radius01c3+real(cell_macr(agg1,2)-1)*radius01m3)
		Vi=(real(cell_macr(agg1,1))*radius01c3+real(cell_macr(agg1,2))*radius01m3)*4*pi/3
		Vf=(real(cell_macr(agg1,1))*radius01c3+real(cell_macr(agg1,2)-1)*radius01m3)*4*pi/3
		deltaEjkr=B*(proportion_cellf*max(adhs12,0.))**(5./3.)*effectiveradius**(4./3.)
		!difficile de prendre la différence d'énergie libre avec le substrat en compte car on ne sait pas encore où le macrophage va sortir.
		!!!!deltaEjkrsubstrat=rep_macr**(5./3.)*((1-proportion_cell3)**(5./3.)*radius0**(4./3.)*p3&
		!!! &-((1-proportion_cellf)**(5./3.)*radius(agg1)**(4./3.)*p1+(1-proportion_cell2)**(5./3.)*radius(agg2)**(4./3.)*p2))
		Fdesag=Vf*(-proportion_cellf*proportion_macrf*max(adhs12,0.))
		Fag=Vi*(-proportion_celli*proportion_macri*max(adhs12,0.))
		if (cell_macr(agg1,1).gt.1) then
		Fag=Vi*(-(proportion_celli)**2*max(adhs12,0.)-proportion_celli*proportion_macri*max(adhs12,0.))
		Fdesag=Vf*(-(proportion_cellf)**2*max(adhs12,0.)-proportion_cellf*proportion_macrf*max(adhs12,0.))
		end if
		deltaF=Fag-Fdesag
		energy_barrier=max(-deltaF,deltaEjkr)
		proba_fragmentation=min(exp(-energy_barrier),1.)
		mm=proba_fragmentation/(tauagg*(radius(agg1)**3+rmm**3)**(1./3.))*dt*cell_macr(agg1,2)
		m=int(mm)!RNBIN(cell_macr(agg1,2),proba_fragmentation)
			if (p.lt.(mm-real(m))) then
			m=m+1
			end if
		m=min(m,cell_macr(agg1,2))
		end if 
		
		c=0
		if (cell_macr(agg1,1).gt.1) then
		!!!! desagregation KPs
		effectiveradius=1/(1/radius0+1/rKP)
	proportion_celli=(real(cell_macr(agg1,1))*radius01c3)/(real(cell_macr(agg1,1))*radius01c3+real(cell_macr(agg1,2))*radius01m3)
	proportion_macri=(real(cell_macr(agg1,2))*radius01m3)/(real(cell_macr(agg1,1))*radius01c3+real(cell_macr(agg1,2))*radius01m3)		
	proportion_cellf=(real(cell_macr(agg1,1)-1)*radius01c3)/(real(cell_macr(agg1,1)-1)*radius01c3+real(cell_macr(agg1,2))*radius01m3)
	proportion_macrf=(real(cell_macr(agg1,2))*radius01m3)/(real(cell_macr(agg1,1)-1)*radius01c3+real(cell_macr(agg1,2))*radius01m3)				
		Vi=(real(cell_macr(agg1,1))*radius01c3+real(cell_macr(agg1,2))*radius01m3)*4*pi/3
		Vf=(real(cell_macr(agg1,1)-1)*radius01c3+real(cell_macr(agg1,2))*radius01m3)*4*pi/3
		deltaEjkr=B*(proportion_cellf*max(adhs11,0.))**(5./3.)*effectiveradius**(4./3.)
		!difficile de prendre la différence d'énergie libre avec le substrat en compte car on ne sait pas encore où le macrophage va sortir.
		!!!!deltaEjkrsubstrat=rep_macr**(5./3.)*((1-proportion_cell3)**(5./3.)*radius0**(4./3.)*p3&
		!!! &-((1-proportion_cellf)**(5./3.)*radius(agg1)**(4./3.)*p1+(1-proportion_cell2)**(5./3.)*radius(agg2)**(4./3.)*p2))
		Fdesag=Vf*(-proportion_cellf*proportion_macrf*max(adhs12,0.))
		Fag=Vi*(-proportion_celli*proportion_macri*max(adhs12,0.))
		if (cell_macr(agg1,1).gt.1) then
		Fag=Vi*(-(proportion_celli)**2*max(adhs12,0.)-proportion_celli*proportion_macri*max(adhs12,0.))
		end if
		if (cell_macr(agg1,1).gt.2) then
		Fdesag=Vf*(-(proportion_cellf)**2*max(adhs12,0.)-proportion_cellf*proportion_macrf*max(adhs12,0.))
		end if
		deltaF=Fag-Fdesag
		energy_barrier0=max(-deltaF,deltaEjkr)
		proba_fragmentation=min(exp(-energy_barrier0),1.)
		cc=proba_fragmentation/(tauagg*(radius(agg1)**3+rKP**3)**(1./3.))*dt*cell_macr(agg1,1)
		c=int(mm)
			if (p0.lt.(cc-real(c))) then
			c=c+1
			end if
		c=min(c,cell_macr(agg1,1))
		end if
													
		
		if (c.eq.cell_macr(agg1,1)) then
			m=cell_macr(agg1,2)
		end if
		
			if ((m.ne.m).or.(c.ne.c)) then
			print*, 'bug fragmentation'
			stop
			end if
			
			
	       cell_macr_intermediaire=cell_macr(agg1,2)
	       cell_macr(agg1,2)=max(cell_macr_intermediaire-m,0)
		   cell_macr_intermediaire=cell_macr(agg1,1)
	       cell_macr(agg1,1)=max(cell_macr_intermediaire-c,0)
		   ntotal(agg1)=cell_macr(agg1,1)+cell_macr(agg1,2)
	       radius(agg1)=(real(cell_macr(agg1,1))*radius01c3+real(cell_macr(agg1,2))*radius01m3)**(1./3.)
	       mass(agg1)=real(cell_macr(agg1,1))*mass_c+real(cell_macr(agg1,2))*mass_m
	
	    n=0
	    m0=0
	    c0=0
		if (m.ge.1) then
	    do while (m0.lt.m)
			n00=n
            n=n00+1
            if (cell_macr(n,1)+cell_macr(n,2).eq.0) then
				cell_macr(n,1)=0
                cell_macr(n,2)=1
				ntotal(n)=1
                radius(n)=rmm
                adh_strength(n,2)=proba_macr;adh_strength(n,1)=0
                mass(n)=mass_m
				beta0=0
			idum=n+agg1+int(radius(agg1)+10*abs(alpha)+beta0)
			alpha=-1+2*ran2(idum)
			idum=idum+1000
			beta=-1+2*ran2(idum)
				s=(alpha**2+beta**2)
				s0=s
				beta0=0
				if (s.ge.1) then
				s=s0-1
				beta=s*alpha
				alpha=sqrt(s-beta**2)
				end if
                positions(n,1)=positions(agg1,1)+(radius(agg1)+radius(n))*2*alpha*sqrt(1-s);
				positions(n,2)=positions(agg1,2)+(radius(agg1)+radius(n))*2*beta*sqrt(1-s);
				positions(n,3)=positions(agg1,3)+(radius(agg1)+radius(n))*(1-2*s)
				if ((positions(n,1).ne.positions(n,1)).or.(positions(n,2).ne.positions(n,2)).or.&
				&(positions(n,3).ne.positions(n,3)).or.(radius(n).ne.radius(n)).or.(radius(agg1).ne.radius(agg1))) then
				print*, 'problem position', s,positions(n,1),positions(n,2),positions(n,3),radius(agg1),radius(n),radius01c3,&
				&radius01m3,cell_macr(agg1,1),cell_macr(agg1,2),m
				stop
				end if 
				m00=m0
				m0=m00+1
            end if
	    end do
		end if
		if (c.ge.1) then
		do while (c0.lt.c)
			n00=n
            n=n00+1
            if (cell_macr(n,1)+cell_macr(n,2).eq.0) then
				cell_macr(n,1)=1
                cell_macr(n,2)=0
				ntotal(n)=1
                radius(n)=rKP
                adh_strength(n,2)=0;adh_strength(n,1)=adh_strength(agg1,1)
                mass(n)=mass_c
				idum=idum+1000
				beta0=0
			alpha=-1+2*ran2(idum)
			idum=idum+1000
			beta=-1+2*ran2(idum)
				s=(alpha**2+beta**2)
				s0=s
				beta0=0
				if (s.ge.1) then
				s=s0-1
				beta=s*alpha
				alpha=sqrt(s-beta**2)
				end if
                positions(n,1)=positions(agg1,1)+(radius(agg1)+radius(n))*2*alpha*sqrt(1-s);
				positions(n,2)=positions(agg1,2)+(radius(agg1)+radius(n))*2*beta*sqrt(1-s);
				positions(n,3)=positions(agg1,3)+(radius(agg1)+radius(n))*(1-2*s)
				if ((positions(n,1).ne.positions(n,1)).or.(positions(n,2).ne.positions(n,2)).or.&
				&(positions(n,3).ne.positions(n,3)).or.(radius(n).ne.radius(n)).or.(radius(agg1).ne.radius(agg1))) then
				print*, 'problem position', s,positions(n,1),positions(n,2),positions(n,3),radius(agg1),radius(n),radius01c3,&
				&radius01m3,cell_macr(agg1,1),cell_macr(agg1,2),m
				stop
				end if 
				c00=c0
				c0=c00+1
            end if
	    end do
	    end if

	    
	    
	end if
	
	return 
	end
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	subroutine boucle2(Nc_init,radius,cell_macr,rcMM,rcm,rmMM,rmm,dt,pi,fcrit,Din,Dext,d,finfty,tau_growth,&
			&ncritical,nc_cells,nc_macr,adh_strength,proba_cell,timeadh,tdelayc,ntotal,mass_c,mass_m,radius_cell,factor&
			&,positions,sum_adh_force,g,mass,amp_noise,lproteins,rc,adh,radius_well,friction,friction0,friction1,friction2&
			&,velocity_previous,i10,nf0,number_corr,tmax,lf,amp_noise2,probacellcr,phimm,radius_macr,layer_macr,theta_table&
			&,steps_persistence)
			
	integer :: Nc_init,ncritical,nc_cells,nc_macr,ntotal(Nc_init),new_agg
	integer,dimension(Nc_init,2) :: cell_macr
	real :: dt,rc,pi,fcrit,Din,Dext,d,finfty,tau_growth,mass_c,radius_cell,factor,phimm,radius_macr,adh_strength(Nc_init,2)
	real,dimension(Nc_init) :: radius,mass
	integer :: idum,time,tmax,agg1
	real, DIMENSION(3) :: gf,nf,normf,newf_pos,new_v,af,sum_adh_force
	real, DIMENSION(Nc_init) :: m
	real :: radius_well,friction(3),amp_noise,lproteins,mm,radius_check,nf0,friction0,friction1,friction2,lf&
	&,amp_noise2
	REAL, DIMENSION(Nc_init,3)  ::	positions,velocity_previous
	real :: number_corr,probacellcr,layer_macr,mass_m,mass_kpm,theta_table(Nc_init,3),steps_persistence,distance_surf
	
do 270 i270=1,nc_cells+nc_macr
	
	agg1=i270
	new_agg=i270
	radius0=sqrt(positions(agg1,1)**2+positions(agg1,2)**2+positions(agg1,3)**2)
	radius1=sqrt(positions(agg1,1)**2+positions(agg1,2)**2)

		call growth(Nc_init,i270,radius,cell_macr,rcMM,rcm,rmMM,rmm,dt,pi,fcrit,Din,Dext,d,finfty,&
		&tau_growth,ncritical,nc_cells,nc_macr,adh_strength,proba_cell,timeadh,tdelayc,ntotal,radius_cell,factor,&
		&probacellcr,phimm,radius_macr,mass,mass_c,mass_m,positions,layer_macr)
	
	
	radius0=sqrt(positions(new_agg,1)**2+positions(new_agg,2)**2+positions(new_agg,3)**2)
	radius1=sqrt(positions(new_agg,1)**2+positions(new_agg,2)**2)
	
	if ((positions(new_agg,3).lt.0).and.((radius0+radius(new_agg)).gt.radius_well)) then
		pos1=positions(new_agg,1)*(radius_well-1-radius(new_agg))/(radius0)
		pos2=positions(new_agg,2)*(radius_well-1-radius(new_agg))/(radius0)
		pos3=positions(new_agg,3)*(radius_well-1-radius(new_agg))/(radius0)
		positions(new_agg,1)=pos1
		positions(new_agg,2)=pos2
		positions(new_agg,3)=pos3
	elseif ((positions(new_agg,3).gt.0).and.((radius1+radius(new_agg)).gt.radius_well)) then
		pos1=positions(new_agg,1)*(radius_well-1-radius(new_agg))/(radius1)
		pos2=positions(new_agg,2)*(radius_well-1-radius(new_agg))/(radius1)
		positions(new_agg,1)=pos1
		positions(new_agg,2)=pos2
	elseif ((positions(new_agg,3).gt.7e3)) then
		pos3=6e3
		positions(new_agg,3)=pos3
	end if	
	
	radius0=sqrt(positions(agg1,1)**2+positions(agg1,2)**2+positions(agg1,3)**2)
	radius1=sqrt(positions(agg1,1)**2+positions(agg1,2)**2)
	if ((positions(agg1,3).lt.0).and.(radius0+radius(agg1).gt.radius_well)) then
	print*, 'problem initial position boucle2growth 1 down',radius0+radius(agg1),radius(agg1)
	stop
	elseif ((positions(agg1,3).gt.0).and.(radius1+radius(agg1).gt.radius_well)) then
	print*, 'problem initial position boucle2growth 1 up',radius1+radius(agg1),radius(agg1)
	stop
	end if

270	continue

		idum=1e5
		positions=positions

		do 271 i271=1,nc_cells+nc_macr

			if (cell_macr(i271,1)+cell_macr(i271,2).gt.0) then

			gf(1)=0;gf(2)=0;gf(3)=0			
			af(1)=0;af(2)=0;af(3)=0
			nf(1)=0;nf(2)=0;nf(3)=0
			normf(1)=0;normf(2)=0;normf(3)=0
			new_v(1)=0;new_v(2)=0;new_v(3)=0
			newf_pos(1)=0;newf_pos(2)=0;newf_pos(3)=0
			sum_adh_force(1)=0;sum_adh_force(2)=0;sum_adh_force(3)=0
			
			mm=m(i271)
			rc=10
			
			call friction_function_new(Nc_init,i271,friction,friction0,friction1,friction2,cell_macr,radius_well,positions,&
			&lf,radius,rc,adh_strength,probacellcr,rcMM,rcm,rmm,timeadh,proba_cell,tdelayc,dt,phimm,layer_macr,&
			&radius_macr)
			
			nn=cell_macr(i271,1)+cell_macr(i271,2)
			mass_kpm=real(cell_macr(i271,1))*mass_c+real(cell_macr(i271,2))*mass_m
			
			call gravity_force(Nc_init,g,mm,gf,mass_kpm) 
			
			if (positions(i271,3).le.0) then
			distance_surf=radius_well-sqrt(positions(i271,1)**2+positions(i271,2)**2+positions(i271,3)**2)-radius(i271)
			elseif (positions(i271,3).gt.0) then
			distance_surf=radius_well-sqrt(positions(i271,1)**2+positions(i271,2)**2)-radius(i271)
			endif
			
			theta_last=theta_table(i271,2)
			amp_noise21_last=theta_table(i271,3)
			idum=i10*i271+1000
			if ((theta_table(i271,1).gt.0.5).or.(distance_surf.gt.10).or.(i10.lt.2).or.(ran2(idum).le.steps_persistence*dt)) then
			theta_last=1e5
			amp_noise21_last=0
			end if
			
			call noise_force(Nc_init,nf,amp_noise,idum,positions,friction,dt,radius_well,i271,cell_macr,&
			&radius_well,amp_noise2,radius,lf,adh_strength,radius_macr,probacellcr,rcMM,rcm,proba_cell,timeadh,tdelayc&
			&,layer_macr,gf,theta_last,amp_noise21_last,steps_persistence)
			theta_table(i271,1)=0
			theta_table(i271,2)=theta_last
			theta_table(i271,3)=amp_noise21_last
		
			call new_velocity_force(i271,nf,gf,normf,sum_adh_force,new_v,Nc_init,friction,friction0,friction1,friction2,dt,&
			&cell_macr,radius_well,positions,lf,radius)

			do 272 i272=1,3
				velocity_previous(i271,i272)=new_v(i272)
272				continue

			call new_position_force(i271,positions,newf_pos,new_v,Nc_init,dt,nf,gf,normf,sum_adh_force,friction,radius,&
			&radius_well)
			
			positions(i271,1)=newf_pos(1)          
			positions(i271,2)=newf_pos(2)
			positions(i271,3)=newf_pos(3)
			
			end if

271 	 continue


	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	 subroutine boucle_fusion(mass,radius,positions,cell_macr,Nc_init,num_agg,ntotal,proba_cell,proba_macr,nc_cells,nc_macr&
	 &,distancefusion,adh_strength,dt,timeadh,timeadhm,tdelayc,tdelaym,rcMM,rcm,rmMM,rmm,rep_macr,mass_c,mass_m,probacellcr,deltaG0,&
	 &tauagg,rKP,radius_macr,phimm,layer_macr,theta_table)

	 real :: positions(1:Nc_init,1:3),mass(1:Nc_init),radius(1:Nc_init),adh_strength(1:Nc_init,1:2),dt,timeadh,timeadhm,&
	 &tdelayc,tdelaym,f0,f1,rc,ntot,proba_cell0,rep_macr,mass_c,mass_m,proba_cellm,probacellcr,deltaG0,tauagg,rKP,radius_macr,phimm,&
	 &pi,Vcellsagg3,Vmax3,phi3m,ran2i,mmax_int,mmax_real,layer_macr,theta_table(Nc_init,3)
	 integer :: num_agg(1:Nc_init),cell_macr(1:Nc_init,1:2),ntotal(1:Nc_init),nc_cells,nc_macr,dejapasse(1:Nc_init,1:Nc_init),agg1,&
	 &ntt0,ntt1,m,mext,m0,n,n0,n00,m00,idum,agg2

	 pi=3.14159

	 f0=0
	 f1=0
		 if (timeadh.gt.0) then
		 f0=proba_cell/timeadh*tdelayc
		 end if 
		 if (timeadhm.gt.0) then
		 f1=proba_macr/timeadhm*tdelaym
		 end if 

	 do 8002 i8002=1,nc_cells+nc_macr
	 if ((adh_strength(i8002,1).gt.(proba_cell+f0)*1.1)) then
				print*, 'problem adh 8002 before',adh_strength(i8002,1),proba_cell+f0,cell_macr(i8002,1),i8002
				stop
				end if
	 if (cell_macr(i8002,1)+cell_macr(i8002,2).gt.1) then          !maturation inside the aggregates
			adh0=max(adh_strength(i8002,1),f0)
			adh1=max(adh_strength(i8002,2),f1)

			ntot=real(cell_macr(i8002,1)+cell_macr(i8002,2))/(1e2)
			
			adh_strength(i8002,1)=min(adh0+(proba_cell+f0-adh0)/timeadh*dt,proba_cell+f0)
			
			adh_strength(i8002,2)=proba_macr
					if (timeadhm.gt.0) then
					adh_strength(i8002,2)=min(adh1+(proba_macr+f1-adh1)/timeadhm*dt,proba_macr+f1)
					end if
					
			radius01c3=r3(adh_strength(i8002,1)-f0,probacellcr,rcMM,rcm)
			radius01m3=rmm**3
			radius(i8002)=(real(cell_macr(i8002,1))*radius01c3+real(cell_macr(i8002,2))*radius01m3)**(1./3.)
	
	 end if
	 
	 call expulsion_macr(i8002,cell_macr,phimm,ntotal,radius,adh_strength,mass,mass_c,mass_m,radius_macr,proba_macr&
	 &,positions,Nc_init,radius01c3,layer_macr,radius_well)

			
			
8002		continue

	 do 8000 i8000=1,nc_cells+nc_macr
	 ntt0=cell_macr(i8000,1)+cell_macr(i8000,2)
	 if (ntt0.ge.1) then 
		do 8001 i8001=1,nc_cells+nc_macr
		
			if ((adh_strength(i8000,1).gt.(proba_cell+f0)*1.1).or.((adh_strength(i8001,1).gt.(proba_cell+f0)*1.1))) then
				print*, 'problem adh before',adh_strength(i8000,1),adh_strength(i8001,1),proba_cell+f0,cell_macr(i8000,1)&
				&,i8000,cell_macr(i8001,1),i8001
				stop
				end if
			ntt1=cell_macr(i8001,1)+cell_macr(i8001,2)
			if ((ntt1.ge.1).and.(i8001.ne.i8000).and.(dejapasse(i8000,i8001).eq.0)) then 
			
			
			agg1=i8000;
			agg2=i8001
					
				call fusion_function(i8000,i8001,mass,radius,positions,cell_macr,Nc_init,&
				&num_agg,ntotal,proba_cell,proba_macr,nc_cells,nc_macr,distancefusion,adh_strength,dt,timeadh,timeadhm,&
				&tdelayc,tdelaym,rcMM,rcm,rmMM,rmm,rep_macr,probacellcr,deltaG0,tauagg,rKP,radius_macr,phimm,mass_c,mass_m,layer_macr&
				&,theta_table)
				
				dejapasse(i8000,i8001)=1
				if ((adh_strength(i8000,1).gt.(proba_cell+f0)*1.1).or.((adh_strength(i8001,1).gt.(proba_cell+f0)*1.1))) then
				print*, 'problem adh',adh_strength(i8000,1),adh_strength(i8001,1),proba_cell+f0,cell_macr(i8000,1)&
				&,i8000,cell_macr(i8001,1),i8001
				stop
				end if
			 end if
8001 	continue 
	end if
8000 continue


	 return
	 end
	 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine growth(Nc_init,agg1,radius,cell_macr,rcMM,rcm,rmMM,rmm,dt,pi,fcrit,Din,Dext,lpenetration,finfty,tau_growth&
&,ncritical,nc_cells,nc_macr,adh_strength,proba_cell,timeadh,tdelayc,ntotal,radius_cell,factor,probacellcr,phimm,radius_macr&
&,mass,mass_c,mass_m,positions,layer_macr)

integer Nc_init,agg1,dcell,idum,ncritical,nc_cells,nc_macr,ntotal(Nc_init)
integer,dimension(Nc_init,2) :: cell_macr
real :: growth_vol,rcritical,dt,rc,factor,probacellcr,l,phi&
&tau_growth,fcrit,Din,Dext,d,finfty,fcrit0,radius_cell,proba_cellm,mass_c,positions(Nc_init,3),mass(Nc_init)
real :: rinit,ddcell,new_random,proportion_cell,adh_strength(Nc_init,2),lpenetration,agg2,rmax,phimm,&
&phic,phicmax,l0,R0,L1,R1,D1,D0,nR0,nR1,radius_macr,leff,layer_macr,mass_m
integer :: cell_init
real,dimension(Nc_init) :: radius
CHARACTER*16 filename
character*8 :: date
character*10 :: time
character*5 :: zone
integer,dimension(8) :: values

if (radius(agg1).ne.radius(agg1)) then
	print*,'problem radius', radius(agg1),cell_macr(agg1,1),cell_macr(agg1,2)
	STOP
end if 

if (cell_macr(agg1,1).gt.0) then
		agg2=agg1
		rinit=radius(agg1)
		cell_init=cell_macr(agg1,1)
		proportion_cell=real(cell_macr(agg1,1))/real(cell_macr(agg1,1)+cell_macr(agg1,2))
		f0=proba_cell/timeadh*tdelayc
		rc3=r3(adh_strength(agg1,1)-f0,probacellcr,rcMM,rcm)
		rmax=(real(cell_macr(agg1,1))*rcm**3+real(cell_macr(agg1,2))*rmm**3*phimm)**(1./3.)
		rad_th=(real(cell_macr(agg1,1))*rc3+real(cell_macr(agg1,2))*rmm**3)**(1./3.)
		phi=(real(cell_macr(agg1,1))*radius_cell**3+real(cell_macr(agg1,2))*rmm**3*phimm)/(radius(agg1)+radius_cell)**3!min(radius_cell**3/rc3,(cell_macr(agg1,1)+cell_macr(agg1,2))*radius_cell**3/(radius(agg1)+radius_cell)**3)
		phimax=(real(cell_macr(agg1,1))*radius_cell**3+real(cell_macr(agg1,2))*rmm**3*phimm)/(rmax+radius_cell)**3
		

			
		fcrit0=fcrit
		
		
		if (cell_macr(agg1,2).gt.0) then
			R0=(real(cell_macr(agg1,1))*rc3)**(1./3.)                          !R0 rayon ext de la zone des KPs
			R1=(real(cell_macr(agg1,1))*rc3+cell_macr(agg1,2)*rmm**3.)**(1./3.)				   !R1 rayon ext de la zone des macr
			R1=R0+layer_macr
			Vmax=4*pi/3*(R1**3-R0**3)
			leff=radius_cell-(R1-R0)
			phic=real(cell_macr(agg1,1))*radius_cell**3./(real(cell_macr(agg1,1))*rc3+max(leff,0.)**3)
			phim=real(cell_macr(agg1,2))*radius_macr**3./(real(cell_macr(agg1,2))*rmm**3.+radius_cell**3)
			phim=real(cell_macr(agg1,2))*radius_macr**3./Vmax

			phicmax=radius_cell**3./min(rcm,rcmm)**3.
			taugdepend=(phicmax-phic)*factor  
			l0=lpenetration*sqrt((1-phic)/phic)					       !penetration length of the KP zone
			l1=lpenetration*sqrt((1-phim)/phim)					   !penetration length of the macr zone
			D1=1/(1-phim)				!D1=Dmacr/Dext
			D0=(1-phim)/(1-phic)					!D0=Dkp/Dmacr
			nR1=1/(1-D1+D1*(R1/l1)*(-1/Tanh((R0-R1)/l1)+(R0/l1)/Sinh((R0-R1)/l1)**2./&
			&(-1+D0-D0*(R0/l0)/Tanh(R0/l0)+(R0/l1)/Tanh((R0-R1)/l1))))		!nR1 concentration de nutriments à la limite agrégat/exterieur
			nR0=nR1*(R1/l1)/Sinh((R0-R1)/l1)/(-1+D0-D0*(R0/l0)/Tanh(R0/l0)+(R0/l1)/Tanh((R0-R1)/l1))					   !nR0 concentration de nutriments à la limite KPs/macrs
		end if
		
		if (cell_macr(agg1,2).eq.0) then
			phic=real(cell_macr(agg1,1))*radius_cell**3./(real(cell_macr(agg1,1))*rc3+radius_cell**3)

			phicmax=real(cell_macr(agg1,1))*radius_cell**3./(real(cell_macr(agg1,1))*min(rcm,rcmm)**3.+radius_cell**3)
			taugdepend=(phicmax-phic)*factor !taugdepend=(phimax-phic)*factor 
			D0=1/(1-phic)
			l0=lpenetration*sqrt((1-phic)/phic)					       !penetration length of the KP zone
			R0=(real(cell_macr(agg1,1))*rc3)**(1./3.)    
			nR0=1/(1-D0+D0*(R0/l0)/Tanh(R0/l0))					   !nR0 concentration de nutriments à la limite KPs/macrs
		end if
		
		growth_vol=tau_growth*(-4./3.*fcrit0*pi*R0**3.-4.*nR0*pi*R0*l0**2*(1-R0/l0/Tanh(R0/l0))+taugdepend*4./3.*pi*(R0**3))
		
		
		rc=rc3**(1.0/3.0)
		dcell=int(dt*growth_vol/((4*pi/3)*rc3))
		ddcell=abs(dt*growth_vol/((4*pi/3)*rc3))-real(dcell)
		
		
		call date_and_time(date,time,zone,values)
		
		idum=int(values(1)+values(2)+values(3)+values(4)+values(5)+values(6)+values(7)+values(8))
		new_random=ran2(idum)
		
			if (new_random.lt.abs(ddcell))  then
			dcell=dcell+1
			end if
		
		cell_macr(agg1,1)=max(dcell+cell_init,0)
		
		if (cell_macr(agg1,1)+cell_macr(agg1,2).eq.0) then
			print*, 'dies in situ',cell_init,phic,phim,l0,l1,rinit,growth_vol,dcell,&
			&ddcell,dcell+cell_init,cell_macr(agg1,1),max(dcell+cell_init,0)
			stop
		end if
		
		radius(agg1)=max((rinit**3+dcell*(rc3))**(1./3.),0.)
		
			if (growth_vol.lt.0) then
				cell_macr(agg1,1)=max(-dcell+cell_init,0)
				radius(agg1)=max((rinit**3-dcell*(rc3))**(1./3.),0.)
			end if
		
		call expulsion_macr(agg1,cell_macr,phimm,ntotal,radius,adh_strength,mass,mass_c,mass_m,radius_macr,proba_macr,positions,Nc_init,&
		&rc3,layer_macr,radius_well)
		

end if


		
return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine time_process(Nc_init,tmax,mass,radius,positions,cell_macr,num_agg,ntotal,proba_cell,proba_macr,dt,g,amp_noise,&
&lproteins,adh,nf0,lf,number_corr,friction0,friction1,friction2,rcMM,rcm,rmMM,rmm,radius_well,&
&pi,fcrit,Din,Dext,d,finfty,tau_growth,ncritical&
&,nc_cells,nc_macr,distancefusion,datadat,adh_strength,timeadh,timeadhm,tdelayc,tdelaym,mass_c,mass_m,radius_cell,factor,&
&amp_noise2,rep_macr,probacellcr,deltaG0,phimm,tauagg,radius_macr,layer_macr,steps_persistence)

integer Nc_init,tmax,ncritical,nc_cells,nc_macr
real, dimension(Nc_init) :: radius,mass
real, dimension(Nc_init,3) :: positions,sum_adh_force,velocity_previous,theta_table
real, dimension(Nc_init,2) :: adh_strength
integer, dimension(Nc_init,2) :: cell_macr
integer, dimension(Nc_init) :: ntotal,num_agg
real :: rc,mass_c,mass_m,radius_well,dt,g,amp_noise,lproteins,adh,friction0,friction1,friction2,nf0,pi,fcrit,Din,Dext,d,finfty,&    
&tau_growth,distancefusion,timeadh,timeadhm,tdelayc,tdelaym,radius_cell,factor,amp_noise2,rep_macr!!!!!!!!!!rc pas sur nécessaire
real :: friction(3),number_corr,lf,probacellcr,deltaG0,phimm,tauagg,radius_macr,layer_macr,steps_persistence
character*1000 :: datadat
number_corr=0.0
pi=4.D0*DATAN(1.D0)

open(111111, file = datadat, status='unknown')  
	 
	 do 10 i10=1,tmax
	 print*, i10
			
			call boucle2(Nc_init,radius,cell_macr,rcMM,rcm,rmMM,rmm,dt,pi,fcrit,Din,Dext,d,finfty,tau_growth,&
			&ncritical,nc_cells,nc_macr,adh_strength,proba_cell,timeadh,tdelayc,ntotal,mass_c,mass_m,radius_cell,factor&
			&,positions,sum_adh_force,g,mass,amp_noise,lproteins,rc,adh,radius_well,friction,friction0,friction1,friction2&
			&,velocity_previous,i10,nf0,number_corr,tmax,lf,amp_noise2,probacellcr,phimm,radius_macr,layer_macr,theta_table&
			&,steps_persistence)
			
			call boucle_fusion(mass,radius,positions,cell_macr,Nc_init,num_agg,ntotal,proba_cell,proba_macr,nc_cells,nc_macr&
			&,distancefusion,adh_strength,dt,timeadh,timeadhm,tdelayc,tdelaym,rcMM,rcm,rmMM,rmm,rep_macr,mass_c,mass_m,probacellcr,deltaG0,&
			&tauagg,radius_cell,radius_macr,phimm,layer_macr,theta_table)
			
	if (mod(int(i10+99),int(1/dt)).eq.0) then
	 write(111111,*) (i10-1)*dt
	 do 351 i351=1,nc_cells+nc_macr
	 f0=proba_cell/timeadh*tdelayc
	 rc3=r3(adh_strength(i351,1)-f0,probacellcr,rcMM,rcm)
	 radius_cells=(real(cell_macr(i351,1))*rc3)**(1./3.)
	 color_number=i351
	 if (cell_macr(i351,1)+cell_macr(i351,2).gt.0) then
	 write(111111,*) positions(i351,1),'  ',positions(i351,2),'  ',positions(i351,3),'  ',radius(i351),'  ',color_number,'  ',&
	 &real(cell_macr(i351,1))/real(cell_macr(i351,1)+cell_macr(i351,2)),'  ',cell_macr(i351,1),'  ',cell_macr(i351,2)&
	 &,'  ',adh_strength(i351,1),'  ',adh_strength(i351,2),'  ',radius_cells
	 end if 
351	 continue
	 write(111111,*) '  '
			end if 
					
10   continue
close(111111)
return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine final_routine(Nc_init,positions,radius,ntotal,cell_macr,number_corr,values,nc_cells,nc_macr,datadat)

integer Nc_init,cell_macr(Nc_init,2),ntotal(Nc_init),nc_cells,nc_macr
real positions(Nc_init,3),radius(Nc_init),number_corr
real :: total,min_relative_distance,max_relative_distance,radiusi2,relative_distance(Nc_init,Nc_init),&
&color_number
CHARACTER*16 filename
character*8 :: date
character*10 :: time
character*5 :: zone
integer,dimension(8) :: values,values_fin
character*1000 :: datadat

print*,'finale positions'
	 do 21 i21=1,nc_cells+nc_macr
	 radiusi21=sqrt(positions(i21,1)**2+positions(i21,2)**2+positions(i21,3)**2)
	 if (radius(i21).gt.1.0) then
	 print*, 'cell', i21, ', position:', '; x:', positions(i21,1), '; y:', positions(i21,2), '; z:', positions(i21,3),&
	 &';radius: ',radiusi21,'radius cells', radius(i21)
	 end if
21	 continue
	 total=0
	 do 32 i32=1,nc_cells+nc_macr
	 total=total+cell_macr(i32,1)
32	 continue

	 
	 min_relative_distance=1e10
	 max_relative_distance=0
	 do 30 i30=1,nc_cells+nc_macr
		 do 31 i31=1,nc_cells+nc_macr
		 
		 if ((radius(i30).gt.11).and.(radius(i31).gt.11).and.(i30.ne.i31)) then
		 relative_distance(i30,i31)=sqrt((positions(i30,1)-positions(i31,1))**2+(positions(i30,2)-positions(i31,2))**2&
		 &+(positions(i30,3)-positions(i31,3))**2)-radius(i30)-radius(i31)
		 
		 if (relative_distance(i30,i31).gt.max_relative_distance) then
		 max_relative_distance=relative_distance(i30,i31)
		 end if 
		 
		 if (relative_distance(i30,i31).lt.min_relative_distance) then
		 min_relative_distance=relative_distance(i30,i31)
		 end if 
		 
		 end if
 31		continue
 30	 continue
	 print*, 'min',min_relative_distance,'max',max_relative_distance
	 	 print*, 'number of corrections',number_corr
call date_and_time(date,time,zone,values_fin)
	 print*, 'heure dÃ©but: ',values(5),'h',values(6),'min',values(7),'s'
	 print*, 'heure fin: ',values_fin(5),'h',values_fin(6),'min',values_fin(7),'s'
	
	close(11111)
	
return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	function r3(adhesion_value,proba,rlarge,rsmall)
	
	Real :: r3,adhesion_value,proba,rsmall,rlarge
	r3=rlarge**3+(rsmall**3-rlarge**3)*(1-exp(-max(adhesion_value,0.)/proba))/(1-exp(-1/proba))
	r3=rlarge**3+(rsmall**3-rlarge**3)*(max(adhesion_value,0.)/proba)
	
	return
	end 
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	function rnbin(n1,p1)
	
	integer :: n1,rnbin,idum
	real :: p1,random0
	
	idum=int(n1*p1)
	
	rnbin=0
	do 9000 i9000=1,n1 
		random0=ran2(idum)
		if (random0.le.p1) then
			rnbin=rnbin+1
		end if
		idum=idum+int(ran2(idum)*n1+1)
9000 continue
	return 
	end
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine expulsion_macr(agg2,cell_macr,phimm,ntotal,radius,adh_strength,mass,mass_c,mass_m,radius_macr,proba_macr,positions&
&,Nc_init,radius02c3,layer_macr,radius_well)

integer :: Nc_init,cell_macr(Nc_init,2),mext,agg2,m0,m00,n,n00,n0,ntotal(Nc_init),j,new_agg
real :: Vcellsagg3,Vmax3,pi,radius_macr,phimm,phi3m,mmax_real,mmax_int,m,radius(Nc_init),adh_strength(Nc_init,2),&
&mass(Nc_init),mass_c,mass_m,radius_well,&
&alpha,beta,beta0,s,s0,positions(Nc_init,3),layer_macr
character*8 :: date
character*1000 :: time,message_erreur
character*5 :: zone
integer,dimension(8) :: values
	radius_well=3000.
	j=0
		if ((cell_macr(agg2,1).gt.0).and.(cell_macr(agg2,2).gt.0)) then
		
				pi=3.14159 
				
				Vcellsagg3=real(cell_macr(agg2,1))*4*pi/3*radius02c3
				Vmax3=4.*pi/3.*((Vcellsagg3/(4*pi/3))**(1./3.)+layer_macr)**3.-Vcellsagg3
				phi3m=real(cell_macr(agg2,2))*(radius_macr**3.)*4*pi/3/(Vmax3)
				
				mmax_real=phimm*Vmax3/(4.*pi/3.)/(radius_macr**3.)
				mmax_int=real(int(mmax_real))
				m=int(mmax_int)
				mext=0
				
				ran2i=ran2(int(Vmax3))
				
				if (ran2i.lt.(mmax_real-mmax_int)) then
					m=m+1
				end if	
				
		if (m.lt.(cell_macr(agg2,2))) then
		
					mext=cell_macr(agg2,2)-m
					cell_macr(agg2,2)=m
					ntotal(agg2)=cell_macr(agg2,1)+cell_macr(agg2,2)
					radius(agg2)=(real(cell_macr(agg2,1))*radius02c3+real(cell_macr(agg2,2))*radius_macr**3./phimm)**(1./3.)
					mass(agg2)=mass_c*real(cell_macr(agg2,1))+mass_m*real(cell_macr(agg2,2))
					
					new_agg=agg2
	radius0=sqrt(positions(new_agg,1)**2+positions(new_agg,2)**2+positions(new_agg,3)**2)
	radius1=sqrt(positions(new_agg,1)**2+positions(new_agg,2)**2)
	
					if ((positions(new_agg,3).lt.0).and.((radius0+radius(new_agg)).gt.radius_well)) then
						pos1=positions(new_agg,1)*(radius_well-1-radius(new_agg))/(radius0)
						pos2=positions(new_agg,2)*(radius_well-1-radius(new_agg))/(radius0)
						pos3=positions(new_agg,3)*(radius_well-1-radius(new_agg))/(radius0)
						positions(new_agg,1)=pos1
						positions(new_agg,2)=pos2
						positions(new_agg,3)=pos3
					elseif ((positions(new_agg,3).gt.0).and.((radius1+radius(new_agg)).gt.radius_well)) then
						pos1=positions(new_agg,1)*(radius_well-1-radius(new_agg))/(radius1)
						pos2=positions(new_agg,2)*(radius_well-1-radius(new_agg))/(radius1)
						positions(new_agg,1)=pos1
						positions(new_agg,2)=pos2
					end if	
					m0=0
					n00=0
					n=0
				
					do while (m0.lt.mext)
						n00=n
						n=n00+1
						if (cell_macr(n,1)+cell_macr(n,2).eq.0) then
						j=0
							cell_macr(n,1)=0
							cell_macr(n,2)=1
							ntotal(n)=1
							radius(n)=radius_macr
							adh_strength(n,2)=proba_macr;adh_strength(n,1)=0
							mass(n)=mass_m
9653						beta0=0
							alpha=-1+2*ran2(int(Vmax3)+n+j)
							beta=-1+2*ran2(int(Vmax3)-n+j)
							s=(alpha**2+beta**2)
							s0=s
							beta0=0  
								if (s.ge.1) then
									s=s0-1
									beta=s*alpha
									alpha=sqrt(s-beta**2)
								end if
								
							positions(n,1)=positions(agg2,1)+(radius(agg2)+radius(n))*2*alpha*sqrt(1-s);
							positions(n,2)=positions(agg2,2)+(radius(agg2)+radius(n))*2*beta*sqrt(1-s);
							positions(n,3)=positions(agg2,3)+(radius(agg2)+radius(n))*(1-2*s)
							rad00=sqrt(positions(agg2,1)**2+positions(agg2,2)**2+positions(agg2,3)**2)+radius(agg2)
							rad0=sqrt(positions(n,1)**2+positions(n,2)**2+positions(n,3)**2)+radius(n)
							rad1=sqrt(positions(n,1)**2+positions(n,2)**2)+radius(n)
							if (j.gt.1e4) then
							print*, 'problem position exp macr',positions(agg2,1),positions(agg2,2),positions(agg2,3),rad00
							stop
							end if
							if ((positions(n,3).lt.0).and.(rad0.gt.radius_well)) then
							j=j+1000
							goto 9653
							elseif ((positions(n,3).gt.0).and.(rad1.gt.radius_well)) then
							j=j+1000
							GOTO 9653
							end if
							
							
							if ((positions(n,1).ne.positions(n,1)).or.(positions(n,2).ne.positions(n,2)).or.&
								&(positions(n,3).ne.positions(n,3)).or.(radius(n).ne.radius(n)).or.(radius(agg2).ne.radius(agg2))) then
								print*, 'problem position', s,positions(n,1),positions(n,2),positions(n,3),radius(agg2),radius(n),&
								&cell_macr(agg2,1),cell_macr(agg2,2),m
								stop
							end if 
							
							m00=m0
							m0=m00+1
						end if
					end do
		end if	
		end if		
return 
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rotate(nposx0,nposy0,nposz0,nposx1,nposy1,nposz1,radiusagg0,radiusagg1,radius_well)

real:: nposx0,nposy0,nposz0,nposx1,nposy1,nposz1,radiusagg0,radiusagg1,alpha,beta,eps,s,&
&dis,disx,disy,disz,theta,phi,deltax,deltay,deltaz,dx,dy,dz,E,radius_well,new_radius,pi,&
&nnposx1,nnposy1,nnposz1,E1,E2,E3,new_dis
integer :: idum0,idum1,j

		idum=0
		theta=0.5
		beta=0.5
		eps=pi/1e3
		j=0
		pi=3.14159
		jj=0
5789	phi=ran2(idum)*2*pi
			
	dis=max(sqrt((nposx0-nposx1)**2.+(nposx0-nposx1)**2.+(nposx0-nposx1)**2.),radiusagg0+radiusagg1+0.5)
			
			deltax=dis*cos(phi)*sin(theta)
			deltay=dis*sin(phi)*sin(theta)
			deltaz=dis*cos(theta)
			
			dx=(nposx1-nposx0);dy=(nposy1-nposy0);dz=(nposz1-nposz0);
			E1=sqrt((dx**2*dz**2+dy**2*dz**2+(-dy**2-dx**2)**2))
			E2=sqrt((dx**2+dy**2))
			E3=sqrt((dx**2+dy**2+dz**2))
			nnposx1=nposx0+(deltax*(dx*dz)/E1-deltay*dy/E2+deltaz*dx/E3)
			nnposy1=nposy0+(deltax*dy*dz/E1+deltay*dx/E2+deltaz*dy/E3)
			nnposz1=nposz0+(deltax*(-dy**2-dx**2)/E1+deltaz*dz/E3)
	
	new_radius=sqrt(nnposx1**2+nnposy1**2+nnposz1**2)
	new_dis=sqrt((nposx0-nnposx1)**2+(nposy0-nnposy1)**2+(nposz0-nnposz1)**2)
	
	
if (new_radius.gt.radius_well) then

	if (jj.lt.3) then
	jj=jj+1
	idum=idum+1000
	goto 5789
	end if
	
	jj=0
	idum=idum+1000
	theta=theta+eps
	j=j+1
	if (theta.gt.pi) then
	theta=pi-theta
	end if
	goto 5789
end if

nposx1=nnposx1;nposy1=nnposy1;nposz1=nnposz1
new_dis=sqrt((nposx0-nposx1)**2+(nposy0-nposy1)**2+(nposz0-nposz1)**2)
	if ((abs(new_dis-dis).gt.1e-3).or.(sqrt((nposx0)**2+(nposy0)**2+(nposz0)**2).gt.radius_well)&
	&.or.(sqrt((nposx1)**2+(nposy1)**2+(nposz1)**2).gt.radius_well)) then
	print*, "problem rotate",new_dis,dis,sqrt((nposx0)**2+(nposy0)**2+(nposz0)**2),&
	&sqrt((nposx1)**2+(nposy1)**2+(nposz1)**2)
	stop
	end if


return 
end
	

	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	
subroutine position_function_new(agg0,agg1,cell0,cell1,macr0,macr1,posx0,posy0,posz0,posx1,posy1,posz1,radius_well&
&,radiusagg0,radiusagg1)

real :: posx,posy,posz,beta0,alpha,beta,s,s0,distance,radius0,radius_well,radius1,radiusagg&
&,variation,radiusagg2,eps,dis,posx0,posy0,posz0,posx1,posy1,posz1,radiusagg0,radiusagg1,&
&weight0,weight1,distance0,distance1,barycentrex,barycentrey,barycentrez,disx00,disy00,disz00,disx10,disy10,disz10,&
&disx0,disy0,disz0,disx1,disy1,disz1,dis10,dis00,rapport,radius0_up,radius0_do,radius1_up,radius1_do,dist_init,&
&nposx0,nposy0,nposz0,nposx1,nposy1,nposz1,rapport0,rapport1,translationx,translationy,translationz,new_dis,a,&
&dist_corr,v10x,v10xy,v10z,displacement,w0,w1
integer :: idum,idum0,idum1,j,cell0,cell1,macr0,macr1,agg0,agg1
		j=0		

		distance=radiusagg0+radiusagg1
		
		idum=abs(int(radiusagg0)+int(distance)+int(posx)+int(posy)+int(posz))
		
		w0=1-(radiusagg0**3)/(radiusagg0**3+radiusagg1**3);w1=1-w0
		dist_init=sqrt((posx0-posx1)**2+(posy0-posy1)**2+(posz0-posz1)**2)
		displacement=max(distance-dist_init,0.)
		v10x=(posx0-posx1)/dist_init;v10y=(posy0-posy1)/dist_init;v10z=(posz0-posz1)/dist_init !vector de agg1 à agg0
		
		
		
	nposx0=posx0+displacement*v10x*w0;nposy0=posy0+displacement*v10y*w0;nposz0=posz0+displacement*v10z*w0
	nposx1=posx1-displacement*v10x*w1;nposy1=posy1-displacement*v10y*w1;nposz1=posz1-displacement*v10z*w1
		
		
		dist_corr=sqrt((nposx0-nposx1)**2+(nposy0-nposy1)**2+(nposz0-nposz1)**2)
		
		radius0_up=sqrt(nposx0**2+nposy0**2)+radiusagg0;radius1_up=sqrt(nposx1**2+nposy1**2)+radiusagg1
		radius0_do=sqrt(nposx0**2+nposy0**2+nposz0**2)+radiusagg0;radius1_do=sqrt(nposx1**2+nposy1**2+nposz1**2)+radiusagg1
		
		a=0
			if ((nposz0.gt.0).and.(radius0_up.gt.radius_well)) then
			a=1
				rapport0=(radius_well-radiusagg0-0.5)/(radius0_up-radiusagg0)
				nposx0=nposx0*rapport0;nposy0=nposy0*rapport0
			elseif ((nposz1.gt.0).and.(radius1_up.gt.radius_well)) then
			a=2
				rapport1=(radius_well-radiusagg1-0.5)/(radius1_up-radiusagg1)
				nposx1=nposx1*rapport1;nposy1=nposy1*rapport1
			elseif ((nposz0.lt.0).and.(radius0_do.gt.radius_well)) then
			a=3
				rapport0=(radius_well-radiusagg0-0.5)/(radius0_do-radiusagg0)
				nposx0=nposx0*rapport0;nposy0=nposy0*rapport0;nposz0=nposz0*rapport0
			elseif ((nposz1.lt.0).and.(radius1_do.gt.radius_well)) then
			a=4
				rapport1=(radius_well-radiusagg1-0.5)/(radius1_do-radiusagg1)
				nposx1=nposx1*rapport1;nposy1=nposy1*rapport1;nposz1=nposz1*rapport1
			elseif (nposz0.gt.7e3) then
			a=5
				nposz0=6e3
			elseif (nposz1.gt.7e3) then
			a=6
				nposz1=6e3
			end if
		
	posx0=nposx0;posy0=nposy0;posz0=nposz0;
	posx1=nposx1;posy1=nposy1;posz1=nposz1
	new_dis=sqrt((posx0-posx1)**2+(posy0-posy1)**2+(posz0-posz1)**2)
	
return 
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	


