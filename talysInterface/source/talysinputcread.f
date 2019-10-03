      subroutine talysinputcread(inputfilename,projectile,tarZ,
     + tarA, exSpin, exParity, exEne)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : June 28, 2004
c | Task  : User input and defaults
c +---------------------------------------------------------------------
c
c ******************** Set defaults and read input *********************
c
c readinput   : subroutine to read input file
c input1      : subroutine to read input for first set of variables
c input2      : subroutine to read input for second set of variables
c input3      : subroutine to read input for third set of variables
c input4      : subroutine to read input for fourth set of variables
c input5      : subroutine to read input for fifth set of variables
c input6      : subroutine to read input for sixth set of variables
c checkkeyword: subroutine to check for errors in keywords
c checkvalue  : subroutine to check for errors in values
c
      character*64 inputfilename
      character*1 projectile
      integer tarZ, tarA, exSpin, exParity
      real exEne

      call readinputcfile(inputfilename)
      call input1cread(projectile, tarZ, tarA, exEne, exSpin, exParity)
      call input2
      call input3
      call input4
      call input5
      call input6
      call checkkeyword
      call checkvalue
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
