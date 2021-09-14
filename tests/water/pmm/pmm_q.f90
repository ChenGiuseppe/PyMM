program calctotale
!gfortran pmm_q.f90 -llapack ./libxtcf.so ./libxdrf.a -o pmm_q

        implicit none 
	integer :: at3, at1, at2, ninte, n, i, j, k,r,nat, tot 
	integer :: asolu3, asolu1, asolu2 
	integer :: id1, ret, frame, nframe, step, g, m, l
	integer :: nsolv, nsolu,nrad, nroot, fframe, lframe, diagok
	integer :: jarr  !stato di arrivo
	integer :: ipart !stato di partenza
	integer, allocatable :: a(:), asolu(:), num(:), delta(:,:)
	real :: xperint, yperint, zperint, verxperpro, veryperpro, verzperpro
	double precision :: proj, proj2, proj_perp
        real :: verxat1qm, veryat1qm, verzat1qm, verxat2qm 
	real :: veryat2qm, verzat2qm, verxperproqm, veryperproqm, verzperproqm, paraqm, parabqm, gammaqm
        real :: verx3qm, very3qm, verz3qm, modat1qm, modat2qm, modperproqm
	real ::  time, box(9), prec, verxat1, veryat1, verzat1, Ar(3,3)
	real :: xcdm, ycdm, zcdm, xperpro, yperpro, zperpro, verxat2, veryat2, verzat2
	real :: xcdmqm, ycdmqm, zcdmqm
	real, allocatable :: p1(:,:), p2(:,:), p3(:,:)
	real :: versx, versy, versz, versx1, versy1, versz1
	integer :: atom1, atom2

	real :: x1p, y1p, z1p, x2p, y2p, z2p, x3p, y3p, z3p, verx3 
	real ::very3, verz3, nverx3, nvery3, nverz3, x1pd, y1pd, z1pd , nx1pd, ny1pd, nz1pd
	real ::  Enerx, Enery, Enerz, para, parab, gamma, modat1, modat2, modperpro
	real :: distx, disty, distz, dist, mu0xqm, mu0yqm, mu0zqm, Enerxqm, Eneryqm, Enerzqm
	real, allocatable :: xint(:), yint(:), zint(:), ch(:), x(:), xintcm(:), yintcm(:), zintcm(:), coorx(:), coory(:), coorz(:) 
	real, allocatable :: mx(:,:), my(:,:), mz(:,:)
	real*8, allocatable:: mass(:), massqm(:)
	real, allocatable :: E0(:), ExMux(:,:), EyMuy(:,:), EzMuz(:,:)
	real*8, allocatable :: Eivec(:,:), Eprim(:), Himp(:,:), Hdiag(:,:), Hprim(:,:)
        real*8, allocatable :: Tvec(:,:), d0x(:), d0y(:), d0z(:) 
	real, allocatable :: xsolu(:), ysolu(:), zsolu(:),VRi(:),Vll(:)
	real*8, allocatable :: Qdiag(:,:)
	real*8 ::  mu0x, mu0y, mu0z
	real, allocatable :: mxp(:,:), myp(:,:), mzp(:,:)
	character(len=50) :: interne, solvente, soluto, topol
	character(len=60) :: traj
	character(len=50) :: file_masse, mdipoli, energie,file_cq_charges
	character(len=50) :: output,nucleo
	character(len=4), allocatable :: atnom(:)
	real :: Potential
	real*8 :: Carica_cq, Mass_cq, mass_cqqm !carica e massa del centro quantistico
        character(len=5) :: sdum
	integer :: indice_1, indice_2
	!variabili legate alle strutture del fit e matrici di rotazione
	real*8, allocatable :: r0_1(:), r0_2(:), rfin(:)
        real*8 the,phi,psi
	real*8 AA(3),BB(3),CC(3)
        real*8 D(3,3)
	!per valutare il CPU time
	real :: tini
	real :: tfin



!Per la subroutine di diagonalizzazione
        CHARACTER          JOBZ, UPLO
        INTEGER            INFO, LDA, LWORK  !, N
        INTEGER LWMAX
        DOUBLE PRECISION   WORK(1000)
        double precision, allocatable :: W(:)

!chiamo cpu_time
call cpu_time(tini)
!!! computational part ... !!!


!leggo i parametri e file di input

!leggo le coordinate interne
	read(*,'(A50)') interne
!leggo la matrice dei dipoli xyz
	read(*,'(A50)') mdipoli 
!leggo il numero di radici 
	read(*,300) nroot
!leggo l'indice del solvente
	read(*,'(A50)') solvente
!leggo l'indice del soluto (contiene anche le masse)
	read(*,'(A50)') soluto
!leggo la topologia
	read(*,'(A50)') topol 
!leggo il numero di atomi
	read(*,300) nat
!leggo le energie 
	read(*,'(A50)') energie
!leggo la traiettoria
	read(*,'(A60)') traj
!leggo i frame
	read(*,300) fframe
	read(*,300) lframe 
!leggo il numero di radici in uscita
	read(*,300) nrad
!leggo lo stato di partenza
	read*, ipart    
!leggo lo stato di arrivo
	read*, jarr     
!leggo gli indici dei dipoli perturbati
	read*, indice_1
	read*, indice_2
!leggo il file dove stampare l'output
	read(*,'(A50)'), output
!leggo la carica del centro quantistico
	read*, Carica_cq
!leggo il file di cariche del cq
	read(*,'(A50)') file_cq_charges
!se maggiore di 0 agli elementi diagonali si aggiunge \sum_{i}q_{i}V(r_{i})
        read*, diagok


!faccio alcune operazioni
	tot=nat*3
	nframe=lframe-fframe+1 

	
!leggo le coordinate interne  
        open (unit=15, file=interne, action='READ')

	read(15,300) ninte
        nsolu=ninte
	allocate(xint(ninte))
        allocate(yint(ninte))
        allocate(zint(ninte))    
        allocate(xintcm(ninte))
        allocate(yintcm(ninte))
        allocate(zintcm(ninte))
	allocate(massqm(ninte))
        allocate(coorx(ninte))
        allocate(coory(ninte))
        allocate(coorz(ninte))
	allocate(r0_1(3*nsolu))

	do n=1,ninte
		read(15,*) massqm(n), xint(n), yint(n), zint(n)
                print*, xint(n),yint(n),zint(n),'coor'   
	end do

	close(15)

	!leggo le cariche dal file delle cariche per ogni stato
        allocate(Qdiag(nrad,ninte))
        open (unit=25, file=file_cq_charges, action='READ')
        !ciclo sugli stati elettronici
        do m=1,nrad
        !ciclo sul numero degli atomi
	read(25,*) (Qdiag(m,n),n=1,ninte)
	do n=1,ninte
	!print*, Qdiag(m,n)
	end do 
        end do
        close(25)


	mass_cqqm=0.0
	     do i=1,ninte
             mass_cqqm=mass_cqqm+massqm(i)
        end do

	 xcdmqm=0.0                            
	 ycdmqm=0.0
	 zcdmqm=0.0
                                  
	do i=1,ninte
        	xcdmqm=xcdmqm+xint(i)*massqm(i)
        	ycdmqm=ycdmqm+yint(i)*massqm(i)
        	zcdmqm=zcdmqm+zint(i)*massqm(i)
	end do
	xcdmqm=xcdmqm/mass_cqqm
	ycdmqm=ycdmqm/mass_cqqm
	zcdmqm=zcdmqm/mass_cqqm

	!calcolo della struttura di riferimento r0_1 su cui fittare tutto, ottenuta dalla geometria qm
	j=1
        do i=1,ninte
        r0_1(j)=xint(i)-xcdmqm
        r0_1(j+1)=yint(i)-ycdmqm
        r0_1(j+2)=zint(i)-zcdmqm
        j=j+3
        enddo

	open(unit=60, action='write')
        write(60,*)'initial structure'
        do i=1, ninte
	!print*, r0_1(3*i-2), r0_1(3*i-1), r0_1(3*i)
        
	!!! scrivo il file con la struttura di riferimento r0_1 >>> fort.60 !!!
	write(60,*) massqm(i),r0_1(3*i-2),r0_1(3*i-1),r0_1(3*i)
        
	end do
	close(60)


	!leggo la matrice dei dipoli in au, del momento angolare e del contributo nucleare per un nucleo 
	allocate(mx(nroot,nroot))
	allocate(my(nroot,nroot))
	allocate(mz(nroot,nroot))

        open(unit=18, file=mdipoli, action='READ')

        do i=1,nroot
        do j=1,nroot
        read(18,*) sdum, sdum, mx(i,j), my(i,j), mz(i,j)
        !print*, i, j, mx(i,j), my(i,j), mz(i,j), 'dipoli' 
        end do
        end do
        close(18)


	!leggo l'indice degli atomi del sistema che genera il campo
	open(unit=20, file=solvente, action='READ')
	read(20,*) nsolv
	allocate(a(nsolv))

	!controllo su nsolv=nat-ninte ovvero se gli atomi totali meno quelli del qc sono uguali a quelli del campo
        if (nsolv.ne.(nat-ninte)) print*,'ERRORE NUMERO ATOMI SOLVENTE'   
	do i=1,nsolv
	read(20,*) a(i)
	end do 
	close(20)
	!print*, 'ho letto indice del solvente'
	
	!leggo l'indice degli atomi su cui devo calcolare il campo
	open(unit=25, file=soluto, action='READ') 
	read(25,*) nsolu
	allocate(asolu(nsolu))
	allocate(mass(nsolu))
	
	!controllo su nsolu==ninte
        if (nsolu.ne.ninte) print*,'ERRORE NUMERO ATOMI SOLUTO'   
        
	Mass_cq=0.0
        do i=1,nsolu
        read(25,*) asolu(i), mass(i) 
        Mass_cq=Mass_cq+mass(i)
        end do
	close(25)

	allocate(r0_2(3*nsolu))
	allocate(rfin(3*nsolu))
	

	!leggo la topologia con le cariche
	allocate(num(nat))
	allocate(atnom(nat))
	allocate(ch(nat))

	open(unit=30, file=topol, action='READ')
	
	do j=1,nat
	read(30,*) sdum, ch(j)
	end do 
	close(30)
	!print*, 'ho letto la topologia'
	
	!leggo il file delle energie
	allocate(E0(nroot))
	allocate(Vll(nroot))

        !Variabile per la subroutine di diagonalizzazione
        allocate(W(nroot))

	open(unit=35, file=energie, action='READ')
	do i=1,nroot
	read(35,*) E0(i)
        !print*, E0(i),i,'energie'
	end do
	!print*, 'ho letto le energie'

	!costruisco una delta di Kroneker che mi servirà poi
	allocate(delta(nroot,nroot))

	do i=1,nroot
	do j=1,nroot
		if (i==j) then
			delta(i,j)=1
		else
			delta(i,j)=0
		end if
	end do
	end do

	!apro la traiettoria 
	!alloco le variabili che uso nella traiettoria
	allocate(x(tot))
	allocate(mxp(nroot,nroot))
	allocate(myp(nroot,nroot))
	allocate(mzp(nroot,nroot))
	allocate(ExMux(nroot,nroot))
        allocate(EyMuy(nroot,nroot))
        allocate(EzMuz(nroot,nroot))
	allocate(Hprim(nroot,nroot))       
	allocate(Hdiag(nrad,nrad))
	allocate(Himp(nrad,nrad))
        allocate(Eprim(nrad))
        allocate(Eivec(nrad,nrad))        
        allocate(Tvec(nrad,nrad))   
        allocate(d0x(nrad))
        allocate(d0y(nrad))
        allocate(d0z(nrad))   
	allocate(xsolu(ninte))
        allocate(ysolu(ninte))
        allocate(zsolu(ninte))
        allocate(VRi(ninte))

	open(unit=11, file=output)
	!open(unit=1111,file='out_dipoli_pert')
	!open(unit=1112,file='out_coordinate')

	call xdrfopen (id1,traj,'r',ret)
	if (ret.eq.0) then
                print*, '> error opening xtc file:'
                print*, '>', traj
                stop
        end if 

	print*, 'ho aperto la traiettoria'

	do frame=1,fframe-1
		call readxtc(id1,nat,step,time,box,x,prec,ret) 
	end do

	do frame=1,nframe
		call readxtc(id1,nat,step,time,box,x,prec,ret) 
	do j=1,tot
	!Trasformazione in unità atomiche
	x(j)=x(j)*10      !/0.52918 
	end do

	!calcolo il centro di massa del centro quantistico nel riferimento GROMACS

                xcdm=0.0
                ycdm=0.0
                zcdm=0.0

                do m=1,nsolu
                        g=asolu(m)
                        xcdm=xcdm+x(3*g-2)*mass(m)
                        ycdm=ycdm+x(3*g-1)*mass(m)
                        zcdm=zcdm+x(3*g)*mass(m)
                end do
                
		xcdm=xcdm/Mass_cq
                ycdm=ycdm/Mass_cq
                zcdm=zcdm/Mass_cq

	!mi scrivo la struttura da confrontare con la struttura di riferimento (r0_1)
	do m=1,nsolu
        g=asolu(m)
        end do

	j=1
        do m=1,nsolu
	g=asolu(m)
        r0_2(j)=x(3*g-2)-xcdm
        r0_2(j+1)=x(3*g-1)-ycdm
        r0_2(j+2)=x(3*g)-zcdm
	j=j+3
        enddo

	!ora fa il fit si fa rispetto alla struttura r0_1
	!per adesso ruota la struttura Gromacs, r0_2, sul riferimento r0_1

        call mwsfit(nsolu,r0_2,r0_1,rfin,AA,BB,CC,mass,mass_cq,the,phi,psi)
	
	write(600,*)'frame'
        do i=1,nsolu
        !!! scrivo la struttura finale ottenuta dal fit >>> fort.600 !!!
	write(600,*) mass(i),rfin(3*i-2),rfin(3*i-1),rfin(3*i)
        enddo
	
	!!! scrivo la matrice di rotazione >>> fort.67 !!!
	write(67,*)'rotation matrix',frame, the, phi, psi
        write(67,*)AA(1),BB(1),CC(1)
        write(67,*)AA(2),BB(2),CC(2)
        write(67,*)AA(3),BB(3),CC(3)

!       do i=1,3
!       do j=1,3
!       D(i,j)=AA(i)*AA(j)+BB(i)*BB(j)+CC(i)*CC(j)
!       enddo
!       enddo

!       write(67,*)'***********'
!       write(67,*)D(1,1),D(1,2),D(1,3)
!       write(67,*)D(2,1),D(2,2),D(2,3)
!       write(67,*)D(3,1),D(3,2),D(3,3)
!       write(67,*)'***********'
	
	
	!Ho ottenuto la matrice AABBCC che fa la rotazione da MD a QM, ora applico la trasposta per ottenere i dipoli da QM a MD
	!la matrice di rotazione scritta in fort.67 permette di passare da r0_2 a r_01, la trasposta fa il contrario
	!viene applicata la trasposta per portare i dipoli dal riferimento QM a MD per perturbarli correttamente
	!calcolo le proiezioni di mu sui tre vettori

	do i=1,nroot
	do j=1,nroot
	mxp(i,j)=mx(i,j)*AA(1)+my(i,j)*AA(2)+mz(i,j)*AA(3)
	myp(i,j)=mx(i,j)*BB(1)+my(i,j)*BB(2)+mz(i,j)*BB(3) 
	mzp(i,j)=mx(i,j)*CC(1)+my(i,j)*CC(2)+mz(i,j)*CC(3)
        end do
	end do
	!sono i dipoli ruotati


	!calcolo il campo elettrico nel sistema della proteina
	!inizializzo il campo elettrico
		Enerx=0.0
		Enery=0.0
		Enerz=0.0
		Potential=0.0
	!print*, 'ho inizializzato il campo elettrico'
	
	!giro sul solvente per calcolare il campo
	!converto le coordinate in au per ottenere il campo in queste unità
	!\frac{1}{4\pi\epsilon_{0}}=1
	xcdm=xcdm/0.52918
	ycdm=ycdm/0.52918
	zcdm=zcdm/0.52918
	
	do m=1,nsolv
		k=a(m)
		distx=((x(3*k-2)/0.52918)-xcdm)**2
                disty=((x(3*k-1)/0.52918)-ycdm)**2
                distz=((x(3*k)/0.52918)-zcdm)**2 

		dist=sqrt(distx+disty+distz)

		Enerx=Enerx+(ch(k)*(xcdm-(x(3*k-2)/0.52918))/(dist)**3)
                Enery=Enery+(ch(k)*(ycdm-(x(3*k-1)/0.52918))/(dist)**3)
                Enerz=Enerz+(ch(k)*(zcdm-(x(3*k)/0.52918))/(dist)**3)
	
		!Calcolo anche il potenziale elettrico nel centro di massa
		Potential=Potential+ch(k)/dist
	end do

	!print*, 'ho calcolato il campo generato dal soluto'
	!print*, Enerx, Enery, Enerz

	!Calcolo il potenziale su ogni atomo V(R_i)
                do m=1,nsolu
                VRi(m)=0.0
                    g=asolu(m)
                    xsolu(m)=x(3*g-2)/0.52918
                    ysolu(m)=x(3*g-1)/0.52918
                    zsolu(m)=x(3*g)/0.52918
                        do n=1,nsolv
                        k=a(n)
                        distx=(x(3*k-2)/0.52918-xsolu(m))**2
                        disty=(x(3*k-1)/0.52918-ysolu(m))**2
                        distz=(x(3*k)/0.52918-zsolu(m))**2
                        dist=sqrt(distx+disty+distz)
                        VRi(m)=VRi(m)+ch(k)/dist
                        end do
                end do


		!\sum_{i}q_{i}V(r_{i}) to obtain the diagonal elements of Hprim
                do l=1,nroot
                Vll(l)=0.0
                        do m=1,nsolu
                        Vll(l)=Vll(l)+Qdiag(l,m)*VRi(m)
                        end do
                end do

		!calcolo il prodotto scalare fra dipolo e campo elettrico
		do i=1,nroot
		do j=1,nroot
			ExMux(i,j)=-Enerx*mxp(i,j)
			EyMuy(i,j)=-Enery*myp(i,j)
			EzMuz(i,j)=-Enerz*mzp(i,j)
		end do
		end do
		!print*, 'ho scritto la matrice dei prodotti scalari'	
		
		!costruisco la matrice da diagonalizzare
		do i=1,nroot
		do j=1,nroot
			Hprim(i,j)=E0(i)*delta(i,j)+ExMux(i,j)+EyMuy(i,j)+EzMuz(i,j)
		end do
		end do

		!se diagok > 0 vengono aggiunti i termini Vll alla diagonale
		do i=1,nroot
                if (diagok > 0) then
                Hprim(i,i)=E0(i)+Vll(i)
		!print*, Hprim(i,i)  
                endif
                end do
		!print*, 'ho costruito la matrice da diagonalizzare'
		!riempio la matrice che diagonalizzo


		do i=1,nrad
		do j=1,nrad
			Hdiag(i,j)=Hprim(i,j)
		end do
		end do

                !chiamo la subroutine che diagonalizza
                
                LWORK = -1
                JOBZ='V'
                UPLO='U'
                LDA=nrad
                LWMAX=1000

                CALL DSYEV( 'V', 'U', nrad, Hprim, LDA, W, WORK, LWORK, INFO )
                LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
                !Solve eigenproblem
                CALL DSYEV( 'V', 'U', nrad, Hprim, LDA, W, WORK, LWORK, INFO )

                !print*, 'ho diagonalizzato'
                !Eprim sono gli autovalori e Eivec gli autovettori
		
		!se diagok > 0 allora l'autovalore i-esimo viene lasciato così
		!se diagok < 0 l'autovalore viene addizionato di Carica_cq*Potential che è il termine di ordine zero
		do i=1,nrad
                if (diagok > 0) then
                Eprim(i)=W(i)
                else
                Eprim(i)=W(i)+Carica_cq*Potential
                endif
                end do
                !do i=1,nrad
                !Eprim(i)=W(i)+Carica_cq*Potential
                !end do

                do i=1,nrad
                do j=1,nrad
                Eivec(i,j)=Hprim(i,j)
                end do
                end do
		
		!print*, 'ho diagonalizzato'


		!calcolo la matrice trasposta di Hdiag 
		!decommentare per la stampa degli autovettori
		write(251,*) 'frame',frame 
		write(251,*) 'autovettori'
		do i=1,nrad
                do j=1,nrad  
                    Tvec(i,j)=Eivec(j,i)
		!!! scrivo gli autovettori >>> fort.251 !!!
		write(251,*) i,j, Eivec(i,j)
                end do
                end do

		!calcolo le componenti del vettore d0 
 
        	do i=1,nrad
                d0x(i)=0.0 
                d0y(i)=0.0 
                d0z(i)=0.0 
                do j=1,nrad
                	d0x(i)=d0x(i)+mxp(i,j)*Eivec(j,indice_1)           
                        d0y(i)=d0y(i)+myp(i,j)*Eivec(j,indice_1) 
                        d0z(i)=d0z(i)+mzp(i,j)*Eivec(j,indice_1)      
                end do
        	end do

        mu0x=0.0
        mu0y=0.0
        mu0z=0.0
        do i=1,nrad 
        	mu0x=mu0x+d0x(i)*Tvec(indice_2,i)
                mu0y=mu0y+d0y(i)*Tvec(indice_2,i)  
                mu0z=mu0z+d0z(i)*Tvec(indice_2,i)   
        end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Da stampare per il calcolo delle cariche perturbate
!       write(1111,*) frame,  mu0x, mu0y, mu0z, Carica_cq
!       do m=1,nsolu
!       write(1112,*) frame, (xsolu(m)-xcdm),(ysolu(m)-ycdm),(zsolu(m)-zcdm)
!       enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
	!!!STAMPA LE FREQUENZE IN UNITA ATOMICHE E I DIPOLI DI TRANSIZIONI 
	!!! PER LO SPETTRO DI ASSORBIMENTO
	!write(11,*) frame, (Eprim(jarr)-Eprim(ipart))/(2*3.14), (mu0x**2+mu0y**2+mu0z**2)
	!!! PER STAMPARE SOLO ENERGIA GROUND Eprim(ipart)          
        write(11,*) frame, Eprim(ipart)
	!!! STAMPA ENERGIA STATO ECCITATO E DIPOLO STATO ECCITATO (mettere 2 2 2 2 in script)
	!write(11,*) frame, Eprim(ipart),(mu0x**2+mu0y**2+mu0z**2)
	!SCRITTURA DIFFERENZA DIPOLI TRANSIZIONE
	!write(11,*) frame, 100*(sqrt(mu0x**2+mu0y**2+mu0z**2)-sqrt(mx(indice_1,indice_2)**2&
	!&+my(indice_1,indice_2)**2+mz(indice_1,indice_2)**2))/1.677                        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	end do
	close(11)

call cpu_time(tfin)
print*, 'CPU time ', tfin - tini 
write(33,*) 'CPU time ', tfin - tini
close(33)

100     format(i6,f12.3,f12.3)
202     format(f12.3,f12.3,f12.3)
320     format(i6)
300	format(i10)
350  	format(3i6)
400	format(i3,3f6.3)
420     format(i6,2x,e12.6,2x,e12.6,2x,e12.6)
421     format(i6,i3,i3,2x,e12.6,2x,e12.6,2x,e12.6)
450	format(i6,f8.3)
1010	format(15i5)
end program calctotale


        subroutine angles(the,phi,psi,A,B,C,Athe,Aphi,Apsi,Bthe,Bphi,Bpsi,Cthe,Cphi,Cpsi)
        implicit none
        real*8 A(3),B(3),C(3),Athe(3),Aphi(3),Apsi(3), Bthe(3),Bphi(3),Bpsi(3),Cthe(3),Cphi(3),Cpsi(3)
	real*8 :: cosphi, sinphi, costhe, sinthe, cospsi, sinpsi
	real*8 :: the, phi, psi

        cosphi=cos(phi)
        sinphi=sin(phi)
        costhe=cos(the)
        sinthe=sin(the)
        cospsi=cos(psi)
        sinpsi=sin(psi)

        A(1)=cosphi*cospsi-sinphi*costhe*sinpsi
        A(2)=-(cosphi*sinpsi+sinphi*costhe*cospsi)
        A(3)=sinphi*sinthe

        B(1)=sinphi*cospsi+cosphi*costhe*sinpsi
        B(2)=cosphi*costhe*cospsi-sinphi*sinpsi
        B(3)=-cosphi*sinthe

        C(1)=sinthe*sinpsi
        C(2)=sinthe*cospsi
        C(3)=costhe

        Athe(1)=sinphi*sinthe*sinpsi
        Aphi(1)=-sinphi*cospsi-cosphi*costhe*sinpsi
        Apsi(1)=-cosphi*sinpsi-sinphi*costhe*cospsi
        Bthe(1)=-cosphi*sinthe*sinpsi
        Bphi(1)=cosphi*cospsi-sinphi*costhe*sinpsi
        Bpsi(1)=-sinphi*sinpsi+cosphi*costhe*cospsi
        Cthe(1)=costhe*sinpsi
        Cphi(1)=0.
        Cpsi(1)=sinthe*cospsi

        Athe(2)=sinphi*sinthe*cospsi
        Aphi(2)=sinphi*sinpsi-cosphi*costhe*cospsi
        Apsi(2)=-cosphi*cospsi+sinphi*costhe*sinpsi
        Bthe(2)=-cosphi*sinthe*cospsi
        Bphi(2)=-sinphi*costhe*cospsi-cosphi*sinpsi
        Bpsi(2)=-cosphi*costhe*sinpsi-sinphi*cospsi
        Cthe(2)=costhe*cospsi
        Cphi(2)=0.
        Cpsi(2)=-sinthe*sinpsi

        Athe(3)=sinphi*costhe
        Aphi(3)=cosphi*sinthe
        Apsi(3)=0.
        Bthe(3)=-cosphi*costhe
        Bphi(3)=sinphi*sinthe
        Bpsi(3)=0.
        Cthe(3)=-sinthe
        Cphi(3)=0.
        Cpsi(3)=0.

        end
      
        subroutine mwsfit(NATOMS,r,rref,rfin,A,B,C,ATmass,ATmasstot,the,phi,psi)

        implicit none !!double precision(a-h,o-z)

	integer :: i, j, lk
        real*8 A(3),B(3),C(3),Athe(3),Aphi(3),Apsi(3)
  	real*8 Bthe(3),Bphi(3),Bpsi(3),Cthe(3),Cphi(3),Cpsi(3)
  	real*8 r(3*natoms),rref(3*natoms), rfin(3*natoms)
	reaL*8 V(natoms,3)
  	real*8 Wthe(natoms,3),Wphi(natoms,3),Wpsi(natoms,3)
  	real*8 r1(3*natoms),ATmass(natoms)
  	real*8 :: ATmasstot
	integer :: natoms, kll
	real*8 :: sum, the, phi, psi, the1, phi1, psi1, DtheL,DphiL,DpsiL
	real*8 :: DtheL1, DphiL1, DpsiL1, DDtheL, DDphiL, DDpsiL, dthe, dphi, dpsi
	real*8 :: ulambda, pot, pot1, deltapot, DGDr, DGDG 


! Il potenziale da minimizzare con una procedura steepest descent è L=sum_i mass_i * (\vec r_i - \vec rref_i)^ 2/ mass_total
! V(i,1) = \partial_x_i L   V(i,2) = \partial_y_i L   V(i,3) = \partial_z_i L (i = indice dell'atomo) 
! La matrice di rotazione con rispetto gli angoli di Euler ROT = (A, B, C), A,B,C sono 3-vettori 

! angoli iniziali 


           the=0.5
           phi=0.5
           psi=0.5


! calcolo matrice rotazione (A,B,C). Athe è la derivata di A con rispetto l'angolo theta, Aphi è
! la derivata di A con rispetto l'angolo phi, ecc.

          call angles(the,phi,psi,A,B,C,Athe,Aphi,Apsi,Bthe,Bphi,Bpsi,Cthe,Cphi,Cpsi)


        do i=1,NATOMS
         do j=1,3
          if(j.eq.1)lk=2
          if(j.eq.2)lk=1
          if(j.eq.3)lk=0
          V(i,j)=2.d0*ATmass(i)*(r(3*i-2)*A(j)+r(3*i-1)*B(j)+r(3*i)*C(j)-rref(3*i-lk))/ATmasstot

          Wthe(i,j)=r(3*i-2)*Athe(j)+r(3*i-1)*Bthe(j)+r(3*i)*Cthe(j)
          Wphi(i,j)=r(3*i-2)*Aphi(j)+r(3*i-1)*Bphi(j)+r(3*i)*Cphi(j)
          Wpsi(i,j)=r(3*i-2)*Apsi(j)+r(3*i-1)*Bpsi(j)+r(3*i)*Cpsi(j)
         enddo
        enddo


        DtheL=0.0
        do i=1,NATOMS
          sum=0.
          do j=1,3
           sum=sum+V(i,j)*Wthe(i,j)
          enddo
          DtheL=DtheL+sum
        enddo
        DphiL=0.0
        do i=1,NATOMS
          sum=0.
          do j=1,3
           sum=sum+V(i,j)*Wphi(i,j)
          enddo
          DphiL=DphiL+sum
        enddo
        DpsiL=0.0
        do i=1,NATOMS
          sum=0.
          do j=1,3
           sum=sum+V(i,j)*Wpsi(i,j)
          enddo
          DpsiL=DpsiL+sum
        enddo

! fisso lo step iniziale

         ulambda=1.d-10

! calcolo del potenziale L

         pot=0.d0
         do i=1,natoms
          do j=1,3
          if(j.eq.1)lk=2
          if(j.eq.2)lk=1
          if(j.eq.3)lk=0
          pot=ATmass(i)*(r(3*i-2)*A(j)+r(3*i-1)*B(j)+r(3*i)*C(j)-rref(3*i-lk))**2/ATmasstot+pot
          enddo
         enddo

         deltapot=1.
         KLL=0
        DO WHILE(dabs(deltapot).gt.1.d-15)
         KLL=KLL+1

!         write(6,*)pot1,deltapot

         if(KLL.ne.1)pot=pot1
! nuovi angoli

        the1=the-ulambda*DtheL
        phi1=phi-ulambda*DphiL
        psi1=psi-ulambda*DpsiL

          call angles(the1,phi1,psi1,A,B,C,Athe,Aphi,Apsi,Bthe,Bphi,Bpsi,Cthe,Cphi,Cpsi)

         pot1=0.

        do i=1,NATOMS
         do j=1,3
          if(j.eq.1)lk=2
          if(j.eq.2)lk=1
          if(j.eq.3)lk=0
          V(i,j)=2.d0*ATmass(i)*(r(3*i-2)*A(j)+r(3*i-1)*B(j)+r(3*i)*C(j)-rref(3*i-lk))/ATmasstot
          Wthe(i,j)=r(3*i-2)*Athe(j)+r(3*i-1)*Bthe(j)+r(3*i)*Cthe(j)
          Wphi(i,j)=r(3*i-2)*Aphi(j)+r(3*i-1)*Bphi(j)+r(3*i)*Cphi(j)
          Wpsi(i,j)=r(3*i-2)*Apsi(j)+r(3*i-1)*Bpsi(j)+r(3*i)*Cpsi(j)

          pot1=ATmass(i)*(r(3*i-2)*A(j)+r(3*i-1)*B(j)+r(3*i)*C(j)-rref(3*i-lk))**2/ATmasstot+pot1

         enddo
        enddo

         deltapot=pot1-pot

        DtheL1=0.0
        do i=1,NATOMS
          sum=0.
          do j=1,3
           sum=sum+V(i,j)*Wthe(i,j)
          enddo
          DtheL1=DtheL1+sum
        enddo

        DphiL1=0.0
        do i=1,NATOMS
          sum=0.
          do j=1,3
           sum=sum+V(i,j)*Wphi(i,j)
          enddo
          DphiL1=DphiL1+sum
        enddo
        DpsiL1=0.0
        do i=1,NATOMS
          sum=0.
          do j=1,3
           sum=sum+V(i,j)*Wpsi(i,j)
          enddo
          DpsiL1=DpsiL1+sum
        enddo

        do i=1,3
         DDtheL=DtheL1-DtheL
         DDphiL=DphiL1-DphiL
         DDpsiL=DpsiL1-DpsiL
         dthe=the1-the
         dphi=phi1-phi
         dpsi=psi1-psi
        enddo

         DGDr=0.d0
         DGDG=0.d0

         DGDr=DDtheL*dthe+DDphiL*dphi+DDpsiL*dpsi
         DGDG=DDtheL*DDtheL+DDphiL*DDphiL+DDpsiL*DDpsiL
! nuovo step

         ulambda=dabs(DGDr/DGDG)

         the=the1
         phi=phi1
         psi=psi1

         DtheL=DtheL1
         DphiL=DphiL1
         DpsiL=DpsiL1

         ENDDO

!         write(6,*)'rotation angles',the,phi,psi
!         write(6,*)'rotation matrix'
!         write(6,*)A(1),B(1),C(1)
!         write(6,*)A(2),B(2),C(2)
!         write(6,*)A(3),B(3),C(3)

         write(6,*)pot1,deltapot

        do i=1,NATOMS
         do j=1,3
          if(j.eq.1)lk=2
          if(j.eq.2)lk=1
          if(j.eq.3)lk=0
          rfin(3*i-lk)=r(3*i-2)*A(j)+r(3*i-1)*B(j)+r(3*i)*C(j)
         enddo
        enddo

         end
