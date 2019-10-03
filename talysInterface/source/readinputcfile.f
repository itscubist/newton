       subroutine readinputcfile(inputfilename)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : February 25, 2012
c | Task  : Read input
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      integer i
      character*64 s1
      character*64 s2
      character*64 s3
      character*64 inputfilename
c
c ************************** User Input ********************************
c
c iso     : counter for isotope
c inline  : input line
c numlines: maximum number of input lines
c nlines  : number of input lines
c
c We read the complete input file first as a set of character strings.
c The actual keywords will be read from these later on. For natural
c elements, the input file only needs to be read once.
c

      open(1, file = trim(inputfilename), status = 'old') 
      if (iso.ne.1) return
      if (nlines.gt.0) return
      i=1
   10 read(1,'(a80)',end=100) inline(i)
      i=i+1
      if (i.gt.numlines) then
        write(*,'(" TALYS-error: Number of input lines exceeds ",i5)')
     +    numlines
        write(*,'(" numlines in talys.cmb should be increased")')
        stop
      endif
      goto 10
  100 nlines=i-1
c
c ************** Convert uppercase to lowercase characters *************
c
c For easy handling of all the input parameters, the whole input is
c converted to lowercase characters, with the exception of filenames or
c other character strings.
c
c convert: subroutine to convert input line from upper case to lowercase
c
      do 110 i=1,nlines
        call convert(i)
  110 continue
      nlines0=nlines
      return
      end
Copyright (C)  2013 A.J. Koning, S. Hilaire and S. Goriely
