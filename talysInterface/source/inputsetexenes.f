      subroutine inputsetexenes(exEne,exSpin,exParity,popV)
c
c +---------------------------------------------------------------------
c | Author: Baran Bodur
c | Date  : 2019-10-11
c | Task  : Set an excitation energy population from c++
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer exSpin, exParity, exBin
      real exEne
      real popV

      exBin = INT(exEne * 10.0 + 0.5) + 1
      PdistJP(exBin,exSpin,exParity)=PdistJP(exBin,exSpin,exParity)+popV
      PdistE(exBin) = PdistE(exBin) + popV
      
      print *, "Energy: ", exEne,"Energy Bin: ", exBin, " Spin: ", 
     + exSpin, " Parity: ", exParity, " Pop. Value: ",
     + PdistJP(exBin,exSpin,exParity)
      
      end
