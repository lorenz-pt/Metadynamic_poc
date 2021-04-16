        program metadi

        parameter(nlattice = 500)
        parameter(length = 55 )
        real*8 dq, hgt, ext, metad, charge, qtrh  !! Metadynamics parameters
        real*8 den,cc,arate,dt, par,dist(nlattice)!! HMC parameters
        real*8 store(length),td_pot,vin


        integer in_time, measures,step            !! Omelyan/ metropolis integer parameters
        common/metain/ qtrh, dq, hgt, ext,arate
        common/potential/td_pot(length), vin

        open(1,file = 'meta_input',status = 'old')
        open(2,file = 'meta_misure_2_500v7',status ='unknown')
        open(3,file = 'tdp_2_500v7',status = 'unknown')

c============================================================================
C Parameters to be passed
c============================================================================
        
        read(1,*) measures                      !! number of measures
        read(1,*) in_time                       !! integration time
        read(1,*) dt                            !! time step of the integration
        read(1,*) par                           !! beta * hbar * chi
        read(1,*) qtrh                          !! treshold value of the charge
        read(1,*) dq                            !! charge step
        read(1,*) hgt                           !! heigth of the potential
        read(1,*) ext                           !! strength of the potential outside the barrier 
        read(1,*) start                         !! 1 for random start, 0 for cold start

        pi = 3.141592653589793

        den = float(nlattice)/(2.*par*pi*pi)    !!potential coefficient
        cc = den*pi                             !!strength coefficient

        step = 1
        arate = 0.0                             !! acceptance rate
        

        call ranstart                           !! initialize random number generator
        call initialize_lattice(start)          !! initial value of the field
        call lattice_grid()                     !! prepare the lattice grid
        call momentum()                         !!initialize the momenta
        call metropolis_update(den,charge,step) !!initialize charge, distance and fields initial valuez

        do i = 1,length
            store(i) = 0.0
        enddo

c----------------------------------------------------------------------------
C START MEASURING------------------------------------------------------------
c----------------------------------------------------------------------------

        do i_measures = 1,measures              

            call momentum()                             !! extraction of the new momenta

            call hmc(dt,in_time,cc)                     !! evolution of the field via Omelyan integrator

            call metropolis_update(den,charge,step)     !! metropolis step which returns the charge value

            k = mod(i_measures,15)                      !! every 20 steps update the time dependent potential

            if (k.EQ.0) then
                call time_dependent_p(charge) 

                if (i_measures.GT.100000) then
                    do ip = 1,length
                        store(ip) = store(ip) + td_pot(ip)
                    enddo
                endif
                        
            endif

            write(2,*) charge, i_measures               !! write the values of the charge on file: meta_misure
        enddo

        store = store/float(10000)
        print *, 'Acceptance rate:', arate/float(measures)
        
        do i = 2,length-1
            write(3,*) store(i), -qtrh + (i-2.)*dq
        enddo

        call ranfinish                                  !! save the seed

        close(1)                                        !! close opened files
        close(2)
        close(3)
        
        end
c============================================================================
C MAIN CODE ENDS HERE =======================================================
c============================================================================

C****************************************************************************

c============================================================================
C Initialize lattice
c============================================================================
            subroutine initialize_lattice(start)
            
            parameter (nlattice = 500)

            common/lattice/y(nlattice)
            

            if (start.eq.0) then
            do i = 1,nlattice             !! Cold start
                y(i) = 0.0
            enddo

            elseif (start.eq.1) then
            do i = 1,nlattice                
                x = 1.0 - 2.*ran2()       !! frand() random between -1 e 1
                y(i) = x
            enddo
            endif                      


            return
            end
c============================================================================
C Lattice grid
c============================================================================
            subroutine lattice_grid()

            parameter (nlattice = 500)
            
            common/move/np(nlattice),nm(nlattice)
              
            do i = 1,nlattice
                np(i) = i + 1
                nm(i) = i - 1
            enddo
            np(nlattice) = 1        !! periodic boundary condition      
            nm(1) = nlattice        

            return
            end
c============================================================================
C Generate momenta (BOX MULLER)
c============================================================================
            subroutine momentum()
            parameter(nlattice = 500)
            real*8 e_k, pp
            common/momenta/pp(nlattice), e_k            
            pi = 3.141592653589793
         

            n_half = nlattice/2
            e_k = 0.0

            do i_mom = 1,n_half                            !!two momenta are extracted for every iteration

                ray = sqrt(-2*log(1.-ran2()))              !! ray 
                phi = 2.*pi*ran2()                         !! angle 

                pp1 = ray*cos(phi)                         !! x1
                pp2=  ray*sin(phi)                         !! x2

                pp(i_mom) = pp1
                pp(n_half + i_mom) = pp2
                e_k = e_k + (pp1*pp1)/2. + (pp2*pp2)/2.    !! initial kinetic energy

                
            enddo
            return
            end
            
c============================================================================
C HMC + Meta (Omelyan)
c============================================================================
            subroutine hmc(dt,in_time,cc)
            
            parameter(nlattice = 500)
            real*8 qtrh, dq, hgt, ext
            real*8 arate,pp,e_k

            common/metain/ qtrh, dq, hgt, ext,arate
            common/lattice/y(nlattice)
            common/momenta/pp(nlattice),e_k
            common/move/np(nlattice),nm(nlattice)
            

            real*8 pi,dpi,lam, dist(nlattice),f,fc,force(nlattice),cc
            real*8 metad,dt

            dpi =  6.283185307
            pi =  3.141592653589793
            lam = 0.1931833275                                          !!lambda Omelyan parameter 
            omely = 1.-2.*lam       

            ddt = (2.*ran2()-1.)*0.005 + 0.025                          !! time step chosen casually around 0.025

c*****************************************************************************
C First step
c*****************************************************************************
            do i = 1,nlattice
                y(i) = y(i) + lam*pp(i)*ddt                             !! first update of the field
            enddo
            do i = 1,nlattice
                dist(i) = y(np(i))-y(i)                                 !! distance
            enddo
            
            fc = metad(dist)                                            !! metad coefficient                                  
            
            do i = 1,nlattice
                f = fc*(-cos(dpi*dist(i)) + cos(dpi*dist(nm(i))))       !! time dependent force
                force(i)=-cc*(sin(dpi*dist(nm(i)))-                     !! force from the derivative of the action
     &          sin(dpi*dist(i)))+f

                pp(i) = pp(i) + (ddt/2.)*force(i)                       !! first update of the momentum
                y(i) = y(i) + pp(i)*ddt*omely                           !! second update of the field
            enddo

            do i = 1,nlattice
                dist(i) = y(np(i))-y(i)                                 !! distance 
            enddo

            fc = metad(dist)                                            !! metad coefficient       

            do i = 1,nlattice
                f = fc*(-cos(dpi*dist(i)) + cos(dpi*dist(nm(i))))
                force(i)=-cc*(sin(dpi*dist(nm(i)))-
     &          sin(dpi*dist(i)))+f

                pp(i) = pp(i) + (ddt/2.)*force(i)
                y(i) = y(i) + pp(i)*ddt*2.*lam
            enddo
            do i = 1,nlattice
                dist(i) = y(np(i))-y(i)
            enddo

            fc = metad(dist)
c*****************************************************************************
C Intermediate steps (repeat, jumping the first step by multiplying by two the time step at the end)
c*****************************************************************************
        do k = 2,in_time-1

            do i = 1,nlattice
                f = fc*(-cos(dpi*dist(i)) + cos(dpi*dist(nm(i))))
                force(i)=-cc*(sin(dpi*dist(nm(i)))-
     &          sin(dpi*dist(i)))+f

                pp(i) = pp(i) + (ddt/2.)*force(i)
                y(i) = y(i) + pp(i)*ddt*omely
            enddo

            do i = 1,nlattice
                dist(i) = y(np(i))-y(i)
            enddo

            fc = metad(dist)
            do i = 1,nlattice
                f = fc*(-cos(dpi*dist(i)) + cos(dpi*dist(nm(i))))
                force(i)=-cc*(sin(dpi*dist(nm(i)))-
     &          sin(dpi*dist(i)))+f

                pp(i) = pp(i) + (ddt/2.)*force(i)
                y(i) = y(i) + pp(i)*ddt*2.*lam
            enddo
            do i = 1,nlattice
                dist(i) = y(np(i))-y(i)
            enddo

            fc = metad(dist)
            
        enddo
c*****************************************************************************
C last step
c*****************************************************************************
         do i = 1,nlattice
                f = fc*(-cos(dpi*dist(i)) + cos(dpi*dist(nm(i))))
                force(i)=-cc*(sin(dpi*dist(nm(i)))-
     &          sin(dpi*dist(i)))+f

                pp(i) = pp(i) + (ddt/2.)*force(i)
                y(i) = y(i) + pp(i)*ddt*omely
            enddo

            do i = 1,nlattice
                dist(i) = y(np(i))-y(i)
            enddo

            fc = metad(dist)
            do i = 1,nlattice
                f = fc*(-cos(dpi*dist(i)) + cos(dpi*dist(nm(i))))       
                force(i)=-cc*(sin(dpi*dist(nm(i)))-
     &          sin(dpi*dist(i)))+f

                pp(i) = pp(i) + (ddt/2.)*force(i)
                y(i) = y(i) + pp(i)*ddt*lam
            enddo
            return
            end

c============================================================================
C Metropolis
c============================================================================
            subroutine metropolis_update(den,charge,step)

            parameter(nlattice = 500)
            parameter(length = 55)                              !!length of the potential

            real*8 dist(nlattice),den,pp
            real*8 qtrh,dq,hgt,ext,td_pot(length),q,charge
            real*8 e_k,e_end,vin,vend,lgf,lgi,arate
            real*8 charge_0, y_0(nlattice),dist_0(nlattice)

            integer step

            common/lattice/y(nlattice)
            common/momenta/pp(nlattice), e_k
            common/metain/qtrh, dq, hgt, ext,arate
            common/move/np(nlattice),nm(nlattice)
            common/potential/td_pot,vin
            
            save charge_0,y_0,dist_0                             !! keep the last values of the field, charge and distance
           
            

            pi = 3.141592653589793  
            dpi =  6.283185307

            vend = 0.0
            charge = 0.0
            e_end = 0.0
           



            do i = 1,nlattice                       !! final potential

                dist(i) = y(np(i)) -y(i)
                si = sin(pi*dist(i))

                vend = vend + si*si*den
            enddo
        
            do i = 1,nlattice

                ti = pp(i)
                e_end = e_end +ti*ti               !! final kinetic energy

            enddo

            e_end = e_end/2.

            if (step.EQ.1)then                      !! initialize the 'storage fields' and parameters

                charge_0 = 0.0                     
                vin = 0.0

                do i = 1,nlattice

                    y_0(i) = y(i)
                    dist_0(i) = dist(i)
                    
                enddo
                vin = vend
            endif

            step = step +1
            
            do i = 1,nlattice
                charge = charge + sin(dpi*dist(i))/dpi  !! charge
            enddo

            index = int((charge+qtrh)/dq + 2.)          !! index associated with the charge
            q = -qtrh-dq + (float(index)-1.)*dq                 !! value of the charge on the grid
            if (index.GE.(length-1)) then                   !! potential outside the barrier (rigth)
                arm = (charge - qtrh)
                vend = vend + ext*arm*arm + td_pot(length-1)
                !print *, 'potential:', vend, vin,index
                !print *, 'kinetice:', e_end, e_k,index
            elseif (index.LE.1) then                    !! potential inside the barrier (left)
                arm = (charge + qtrh)
                vend = vend + ext*arm*arm + td_pot(2)
                !print *, 'potential:', vend, vin
               ! print *, 'kinetice:', e_end, e_k,index
            else                                        !! potential inside the barrier (triangular potential)
                vend = vend + td_pot(index) +
     &          (td_pot(index+1) - td_pot(index))*(charge-q)/dq
            endif

            lgi = -e_k-vin                              !! logarithm of the initial probability
            lgf = -vend-e_end                           !!                  final
            
            

            if (log(ran2()).LT.(lgf-lgi)) then          !! accept reject with the replacement of the
                charge_0 = charge                       !! storage variables
                vin = vend
                do i = 1,nlattice
                    dist_0(i) = dist(i)
                    y_0(i) = y(i)
                   
                enddo
                    arate = arate +1.                   !! acceptance rate
                   
            endif

            do i = 1,nlattice
                dist(i)  = dist_0(i)
                y(i) =y_0(i)
            enddo
            charge = charge_0

            return
            end

c============================================================================
C Biased potential construction
c============================================================================  
            subroutine time_dependent_p(charge)
            parameter(length = 55)
            
            real*8 charge,qtrh,dq,ext,hgt,td_pot,q
            real*8 vin, arate,frac
            common/metain/qtrh, dq, hgt, ext,arate
            common/potential/td_pot(length),vin

            index = int((charge+qtrh)/dq + 2.)              !!find the position on the 'charge lattice'
            q = -qtrh-dq + (float(index)-1.)*dq                    !!compute the charge value on the grid
            frac = (charge-q)/dq
         
            
            if (index.LT.(length-1).AND.index.GT.1) then   
                tdpi = td_pot(index)
                tdpip1 = td_pot(index+1)
                !! remove the last value of the time dependent potential from
                !! the potential associated to the last accepted path
                vin = vin-tdpi-(tdpip1-tdpi)*frac 

                !! update the potential
                tdpi=tdpi+hgt*(1.-frac)  
                tdpip1=tdpip1+hgt*frac

                !! add the new value of the time dependent potential from
                !! the potential associated to the last accepted path
                vin = vin+tdpi+(tdpip1-tdpi)*frac

                td_pot(index)= tdpi 
                td_pot(index+1) = tdpip1
            endif

            if (index.EQ.(length-1))then
                tdpi = td_pot(index)
                tdpip1 = td_pot(index+1)
                vin = vin-tdpi                             
                !! update the potential
                tdpi=tdpi+hgt*(1.-frac)  
                tdpip1=tdpip1+hgt*frac
                !! add the new value of the time dependent potential from
                !! the potential associated to the last accepted path
                vin = vin+tdpi
                td_pot(index)= tdpi 
                td_pot(index+1) = tdpip1
            endif

            if (index.EQ.1)then
                tdpi = td_pot(index)
                tdpip1 = td_pot(index+1)
                vin = vin-tdpip1                             
                !! update the potential
                tdpi=tdpi+hgt*(1.-frac)  
                tdpip1=tdpip1+hgt*frac
                !! add the new value of the time dependent potential from
                !! the potential associated to the last accepted path
                vin = vin+tdpip1
                td_pot(index)= tdpi 
                td_pot(index+1) = tdpip1
            endif


            return
            end


c============================================================================
C strength construction
c============================================================================ 
        real*8  function metad(dist)

        parameter(nlattice = 500)
        parameter(length = 55)                  

        real*8 dist(nlattice),vi,arate
        real*8 charge,qtrh,dq,hgt,ext,td_pot

        common/metain/qtrh, dq, hgt, ext,arate
        common/potential/td_pot(length),vin
        
        integer i_c, index

        
        charge = 0.0
        dpi = 6.283185307

        do i_c = 1,nlattice                             !!compute the charge
            charge = charge + sin(dpi*dist(i_c))
        enddo
        charge = charge/dpi

        index = int((charge+qtrh)/dq + 2.)              !!find the position on the 'charge lattice'
        

        if (index.GE.(length-1) )then                       !! return the time dependent strength coefficient
            metad = -2.*ext*(charge - qtrh)
        elseif (index.LE.1) then
            metad = -2.*ext*(charge + qtrh)
          
        else
            metad = (td_pot(index)-td_pot(index+1))/dq
        endif

        return
        end

c============================================================================
C Casual number generator ran2
c============================================================================
      function ran2()
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real ran2,am,eps,rnmx
      parameter(im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1, 
     &         ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,
     &         ir2=3791,ntab=32,ndiv=1+imm1/ntab,eps=1.2e-7,
     &         rnmx=1.-eps)
      integer idum2,j,k,iv,iy
      common /dasav/ idum,idum2,iv(ntab),iy
c      save iv,iy,idum2
c      data idum2/123456789/, iv/NTAB*0/, iy/0/

      if(idum.le.0) then
            idum=max0(-idum,1)
            idum2=idum
            do j=ntab+8,1,-1
                  k=idum/iq1
                  idum=ia1*(idum-k*iq1)-k*ir1
                  if(idum.lt.0) idum=idum+im1
                  if(j.le.ntab) iv(j)=idum
            enddo
      iy=iv(1)
      endif
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if(idum.lt.0) idum=idum+im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if(idum2.lt.0) idum2=idum2+im2
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1) iy=iy+imm1
      ran2=min(am*iy,rnmx)

      return
      end

      subroutine ranstart
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      common /dasav/ idum,idum2,iv(32),iy

      open(unit=23, file='randomseed', status='unknown')
      read(23,*) idum
      read(23,*,end=117) idum2
      do i=1,32
         read(23,*) iv(i)
      enddo
      read(23,*) iy
      close(23)
      goto 118                          !!takes account of the first start
  117 if(idum.ge.0) idum = -idum -1     !!
      close(23)
  118 continue                          !!

      return
      end

      subroutine ranfinish
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      common /dasav/ idum,idum2,iv(32),iy

      open(unit=23, file='randomseed', status='unknown')
      write(23,*) idum
      write(23,*) idum2
      do i=1,32
        write(23,*) iv(i)
      enddo
      write(23,*) iy
      close(23)

      return
      end

