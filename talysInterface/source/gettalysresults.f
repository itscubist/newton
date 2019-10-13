
      subroutine gettalysresults(jz,jn,jpar,jexc,jexcout,nnz,nnn,nnpar,
     + nnexc,nnexcout, excprobs, excenes1, excbinw1, excenes2, excbinw2,
     + sepenes)
c
c +---------------------------------------------------------------------
c | Author: Baran Bodur
c | Date  : 2019-09-30
c | Task  : Get Calculated Decay Probabilities and Energies From TALYS
c +---------------------------------------------------------------------
c ******************** Set defaults and read input *********************
c
c  
c
c
c Declarations and Include talys.cmb
      include "talys.cmb" 

      integer nnz, nnn, nnpar, nnexc, nnexcout
      integer jz, jn, jpar, jexc, jexcout
      real excprobs, excenes1, excbinw1, excenes2, excbinw2, sepenes
c      real excprobs(0:numZchan,0:numNchan,0:numpar,0:numex+1,0:numex+1) 
c      real excbinw(0:numZchan,0:numNchan,0:numex)
c      real excenes(0:numZchan,0:numNchan,0:numex+1)
c      real sepenes(0:numZchan,0:numNchan,0:numpar)
      

c Will need these for getting array sizes in c++
      nnz = numZchan
      nnn = numNchan
      nnpar = numpar+1
      nnexc = numex
      nnexcout = numex

c Probability of decay
      excprobs = feedexcl(jz,jn,jpar,jexc,jexcout)

c To Print if Required for Debugging:

c      if(.not.excprobs.eq.0) then
c        print *,"Z: ", jz, " N: ", jn, " Exc: ", jexc, " via: ",jpar, 
c     + " to: ", jexcout, " probabilitiy: ", excprobs
c      endif

c Initial Energy and Bin Width
      excenes1 = Ex(jz,jn,jexc) 
      if(jexc.le.nnexc) then 
        excbinw1 = deltaEx(jz,jn,jexc)
      else 
        excbinw1 = 0
      endif
c Final Energy and Bin Width
      excenes2 = Ex(jz,jn,jexcout) 
      if(jexcout.le.nnexcout) then 
        excbinw2 = deltaEx(jz,jn,jexcout)
      else 
        excbinw2 = 0
      endif
      sepenes = S(jz,jn,jpar)

c Correct for size = last element + 1, and also counting of original
c state as an excited state
      nnexc = nnexc+2
      nnexcout = nnexcout+1

c Loop to get parameters
c      do jz = 0,nz
c        do jn = 0,nn
c          do jpar = 0,npar
c            sepenes(jz,jn,jpar) = S(jz,jn,jpar)
c            do jexc = 0,nexc+1
c              excenes(jz,nz,jexc) = Ex(jz,nz,jexc)
c              if(jexc.le.nexc) excenes(jz,nz,jexc) = deltaEx(jz,nz,jexc)
c              do jexcout = 0,nexcout+1
c                excprobs(jz,jn,jpar,jexc,jexcout) = feedexcl(jz,jn,jpar,jexc,jexcout)
c              end do
c            end do
c          end do
c        end do
c      end do

      return
      end
