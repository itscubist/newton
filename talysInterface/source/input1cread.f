      subroutine input1cread(projectile,tarZ,tarA,exEne,exSpin,exParity)
c
c +---------------------------------------------------------------------
c | Author: Arjan Koning
c | Date  : December 12, 2016
c | Task  : Read input for first set of variables
c +---------------------------------------------------------------------
c
c ****************** Declarations and common blocks ********************
c
      include "talys.cmb"
      character*1 projectile
      integer tarZ, tarA, exSpin, exParity, exBin
      real exEne

      logical      projexist,massexist,elemexist,enerexist,lexist,fexist
      character*1  ch
      character*10 afile
      character*15 bestchar
      character*40 bestpath
      character*80 word(40),key,value,bestfile
      integer      i,i2,inull,ipath,type,iz,nbest,J,k,nen,pbeg,parity,
     +             inum,negrid,lenbest,istat
      real         Ein,etmp,enincF,deninc,E
c
c ************ Read first set of variables from input lines ************
c
c energyfile : file with incident energies
c projexist  : logical for existence of projectile
c massexist  : logical for existence of mass
c elemexist  : logical for existence of element
c enerexist  : logical for existence of energy
c flagnatural: flag for calculation of natural element
c flagmicro  : flag for completely microscopic Talys calculation
c flagastro  : flag for calculation of astrophysics reaction rate
c flagbest   : flag to use best set of adjusted parameters
c flagbestbr : flag to use best set of branching ratios
c flagbestend: flag to put best set of parameters at end of input file
c bestpath   : alternative directory for best values
c eninc      : incident energy in MeV
c numJ       : maximal J-value
c EdistE     : energy of population distribution
c PdistE     : population distribution, spin-independent
c PdistJP    : population distribution per spin and parity
c Ztarget    : charge number of target nucleus
c Ltarget    : excited level of target
c Starget    : symbol of target nucleus
c enincF     : final incident energy
c deninc     : incident energy increment
c Emaxtalys  : maximum acceptable energy for TALYS
c Estop      : incident energy above which TALYS stops
c
c 1. Initializations
c
      energyfile='                                                     '
      projexist=.false.
      massexist=.false.
      elemexist=.false.
      enerexist=.false.
      flagnatural=.false.
      flagmicro=.false.
      flagastro=.false.
      flagbest=.false.
      flagbestbr=.true.
      flagbestend=.false.
      bestpath='                                        '
      do 5 i=0,numen6+2
        eninc(i)=0.
    5 continue
      do 7 i=0,numex
        EdistE(i)=0.
        PdistE(i)=0.
        do 8 parity=-1,1,2
          do 9 J=0,numJ
            PdistJP(i,J,parity)=0.
    9     continue
    8   continue
    7 continue
      Ztarget=0
      Ltarget=0
      Starget='  '
      enincF=0.
      deninc=0.
      Estop=Emaxtalys
c
c nlines     : number of input lines
c getkeywords: subroutine to retrieve keywords and values from input
c              line
c inline     : input line
c word       : words on input line
c key        : keyword
c value      : value or string
c ch         : character
c
c The keyword is identified and the corresponding values are read.
c Erroneous input is immediately checked. The keywords and number of
c values on each line are retrieved from the input.
c
      do 10 i=1,nlines
        call getkeywords(inline(i),word)
        key=word(1)
        value=word(2)
        ch=word(2)(1:1)
c
c 2. The projectile is read
c
c ptype0: type of incident particle
c
        if (key.eq.'projectile') then
          projexist=.true.
c c++ overwrite
          ptype0=projectile
          goto 10
        endif
        
c
c 3. The target mass is read
c
c Atarget: mass number of target nucleus
c
        if (key.eq.'mass') then
          massexist=.true.
          read(value,*,end=500,err=500) Atarget
c   c++ overwrite
          Atarget = tarA 
          goto 10
        endif
c
c 4. The nuclear symbol or charge number is read
c
c numelem: number of elements
c
        if (key.eq.'element') then
          elemexist=.true.
          if (ch.ge.'0'.and.ch.le.'9') then
            read(value,*,end=500,err=500) Ztarget
c c++ overwrite
            Ztarget = tarZ
            if (Ztarget.lt.1.or.Ztarget.gt.numelem) goto 500
            goto 10
          else
            read(value,'(a2)',end=500,err=500) Starget
            Starget(1:1)=achar(iachar(Starget(1:1))-32)
            goto 10
          endif
        endif
c
c 5. The level of the target is read
c
        if (key.eq.'ltarget') then
          read(value,*,end=500,err=500) Ltarget
          goto 10
        endif
c
c 6. The incident energy or file with incident energies is read
c
        if (key.eq.'energy') then
          enerexist=.true.
          if ((ch.ge.'0'.and.ch.le.'9').or.ch.eq.'.') then
            read(value,*,end=500,err=500) eninc(1)
            read(word(3),*,iostat=istat) enincF
            if (istat.ne.0) goto 10
            read(word(4),*,iostat=istat) deninc
            if (istat.ne.0) goto 10
            goto 10
          else
            eninc(1)=0.
            energyfile=value
            goto 10
          endif
        endif
        if (key.eq.'estop') then
          read(value,*,end=500,err=500) Estop
          goto 10
        endif
c
c 7. Test for completely microscopic and/or astrophysical calculation
c
        if (key.eq.'micro') then
          if (ch.eq.'n') flagmicro=.false.
          if (ch.eq.'y') flagmicro=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 500
          goto 10
        endif
        if (key.eq.'astro') then
          if (ch.eq.'n') flagastro=.false.
          if (ch.eq.'y') flagastro=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 500
          goto 10
        endif
c
c 8. Possibility to use "best" adjusted parameter sets
c
        if (key.eq.'best') then
          if (ch.eq.'n') flagbest=.false.
          if (ch.eq.'y') flagbest=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 500
          goto 10
        endif
        if (key.eq.'bestbranch') then
          if (ch.eq.'n') flagbestbr=.false.
          if (ch.eq.'y') flagbestbr=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 500
          goto 10
        endif
        if (key.eq.'bestend') then
          if (ch.eq.'n') flagbestend=.false.
          if (ch.eq.'y') flagbestend=.true.
          if (ch.ne.'y'.and.ch.ne.'n') goto 500
          goto 10
        endif
        if (key.eq.'bestpath') then
          bestpath=value
          goto 10
        endif
   10 continue
c
c The four main keywords MUST be present in the input file.
c
      if (.not.projexist) then
        write(*,'(" TALYS-error: projectile must be given")')
        stop
      endif
      if (.not.massexist) then
        write(*,'(" TALYS-error: mass must be given")')
        stop
      endif
      if (.not.elemexist) then
        write(*,'(" TALYS-error: element must be given")')
        stop
      endif
      if (.not.enerexist) then
        write(*,'(" TALYS-error: energy must be given")')
        stop
      endif
c
c Manual input of structure path and null device.
c
c path   : directory containing structure files to be read
c inull  : counter for null device
c nulldev: null device
c
      do 50 i=1,nlines
        call getkeywords(inline(i),word)
        key=word(1)
        value=word(2)
        ch=word(2)(1:1)
        ipath=0
        if (key.eq.'strucpath') then
          do 60 i2=11,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              ipath=ipath+1
              path(ipath:ipath)=ch
            endif
   60     continue
        endif
        inull=0
        if (key.eq.'nulldev') then
          nulldev='             '
          do 70 i2=9,80
            ch=inline(i)(i2:i2)
            if (ch.ne.' ') then
              inull=inull+1
              nulldev(inull:inull)=ch
            endif
   70     continue
        endif
   50 continue
c
c Test to check accessibility of structure files and null device
c
      ipath=len_trim(path)
      if (path(ipath:ipath).ne.'/') then
        ipath=ipath+1
        path(ipath:ipath)='/'
      endif
      inquire (file=trim(path)//'abundance/H.abun',exist=lexist)
      if (.not.lexist) then
        write(*,'(" TALYS-error: Structure database not installed:",
     +    " change path in machine.f or strucpath keyword",
     +    " in input file")')
        stop
      endif
      if (inull.gt.13) then
        write(*,'(" TALYS-warning: null device should contain 13",
     +    " characters or less")')
      endif
c
c ************* Process first set of input variables *******************
c
c 1. Identification of target and initial compound nucleus
c
c nuc: symbol of nucleus
c
      if (Ztarget.eq.0) then
        do 110 iz=1,numelem
          if (nuc(iz).eq.Starget) then
            Ztarget=iz
            goto 200
          endif
  110   continue
      else
        Starget=nuc(Ztarget)
      endif
c
c Special case for loop over fission fragments or residual
c products to be evaporated
c
c flagffruns: flag to denote that run is for fission fragment
c flagrpruns: flag to denote that run is for residual product
c
  200 if (flagffruns) then
        Ztarget=Zff
        Atarget=Aff
        energyfile='ff000000.ex'
      endif
      if (flagrpruns) then
        Ztarget=Zrp
        Atarget=Arp
        energyfile='rp000000.ex'
      endif
      if (flagffruns.or.flagrpruns) then
        Starget=nuc(Ztarget)
        Ltarget=0
        ptype0='0'
        write(energyfile(3:5),'(i3.3)') Ztarget
        write(energyfile(6:8),'(i3.3)') Atarget
      else
        Ztarget0=Ztarget
        Atarget0=Atarget
      endif
c
c 2. Assignment of index k0 to incident particle
c
c flaginitpop: flag for initial population distribution
c parsym     : symbol of particle
c k0         : index of incident particle
c
c Throughout TALYS, the initial particle can always be identified
c by the index k0.
c
      flaginitpop=.false.
      k0=0
      do 210 type=0,6
        if (ptype0.eq.parsym(type)) then
          k0=type
          goto 220
        endif
  210 continue
c
c It is also possible to define a population distribution as the
c initial state, through the keyword projectile 0. In that case,
c we assign a photon projectile.
c
  220 if (ptype0.eq.'0') flaginitpop=.true.
c
c A calculation for a natural element is specified by target mass 0
c
c abundance: subroutine for natural abundances
c iso      : counter for isotope
c nbest    : number of lines in best file
c isotope  : isotope number of residual nucleus
c Ntarget  : neutron number of target nucleus
c Zinit    : charge number of initial compound nucleus
c parZ     : charge number of particle
c Ninit    : neutron number of initial compound nucleus
c parN     : neutron number of particle
c Ainit    : mass number of initial compound nucleus
c
      if (Atarget.eq.0) then
        flagnatural=.true.
        if (iso.eq.1) then
          call abundance
        else
          if (flagbest) then
            nbest=nlines-nlines0
            do 230 i=nbest+1,nlines
              inline(i-nbest)=inline(i)
  230       continue
            nlines=nlines0
          endif
        endif
        Atarget=isotope(iso)
      endif
      Ntarget=Atarget-Ztarget
      Zinit=Ztarget+parZ(k0)
      Ninit=Ntarget+parN(k0)
      Ainit=Zinit+Ninit
c
c 3. Determine incident energies
c
c Ein: incident energy
c nin: counter for incident energies
c fexist: flag for energy grid
c
c 1. Normal case of a projectile with a target nucleus
c
      Ein=0.
      fexist=.false.
      if (.not.flaginitpop) then
c
c A. If no incident energy is given in the input file, incident energies
c    should be read from a file.
c
c incidentgrid: subroutine with predefined incident energy grids
c
        eninc(0)=0.
        if (eninc(1).eq.0.) then
          inquire (file=energyfile,exist=lexist)
          if (.not.lexist) then
            call incidentgrid(energyfile(1:14),fexist)
            if (fexist) then
              goto 300
            else
              write(*,'(" TALYS-error: give a single incident energy ",
     +          "in the input file using the energy keyword ")')
              write(*,'(14x,"or specify a range of incident energies ",
     +          "in a file ")')
              write(*,'(14x,"or give a correct name for a pre-defined",
     +          " energy grid ",a73)') energyfile
              stop
            endif
          endif
          nen=0
          open (unit=2,file=energyfile,status='old')
  310     read(2,*,end=320,err=510) Ein
          if (Ein.ne.0.) then
            nen=nen+1
c
c There is a maximum number of incident energies
c
c numenin : maximal number of incident energies
c
            if (nen.gt.numenin) then
              write(*,'(" TALYS-error: there are more than",i4,
     +          " incident energies in file ",a73)') numenin,energyfile
              write(*,'(" numenin in talys.cmb should be increased")')
              stop
            endif
            eninc(nen)=Ein
          endif
          goto 310
  320     close (unit=2)
          if (nen.eq.0) then
            write(*,'(" TALYS-error: there are no",
     +        " incident energies in file ",a73)') energyfile
            stop
          endif
c
c Sort incident energies in ascending order and remove double points
c
c numinc: number of incident energies
c
          do 322 i=1,nen
            do 324 k=1,i
              if (eninc(i).ge.eninc(k)) goto 324
              etmp=eninc(i)
              eninc(i)=eninc(k)
              eninc(k)=etmp
  324       continue
  322     continue
          numinc=nen
          do 326 i=1,nen-1
            if (eninc(i).eq.eninc(i+1)) then
              do 328 k=i+1,nen
                eninc(k)=eninc(k+1)
  328         continue
              numinc=numinc-1
            endif
  326     continue
c
c The minimum and maximum value of all the incident energies is
c determined.
c
c enincmin: minimum incident energy
c enincmax: maximum incident energy
c
          enincmin=eninc(1)
          enincmax=eninc(1)
          do 330 nen=2,numinc
            enincmin=min(enincmin,eninc(nen))
            enincmax=max(enincmax,eninc(nen))
  330     continue
        else
          if (enincF.eq.0.) then
c
c B1. Single value given in the user input file
c
            numinc=1
            enincmin=eninc(1)
            enincmax=eninc(1)
          else
c
c B2. Energy grid based on input values E1, E2, dE
c
c negrid: number of grid points
c
            if (enincF.le.eninc(1)) then
              write(*,'(" TALYS-error: final incident energy should",
     +          " be larger than the first ",a80)') inline(i)
              stop
            endif
            if (deninc.lt.0.) then
              write(*,'(" TALYS-error: energy increment should",
     +          " be larger than zero ",a80)') inline(i)
              stop
            endif
            if (deninc.eq.0.) then
              negrid=10
              deninc=(enincF-eninc(1))/(negrid-1)
            endif
            nen=1
  335       nen=nen+1
            if (nen.gt.numenin) then
              write(*,'(" TALYS-error: number of incident energies ",
     +          " greater than ",i4)') numenin
              stop
            endif
            E=eninc(nen-1)+deninc
            if (E.lt.enincF-1.e-4) then
              eninc(nen)=E
              goto 335
            else
              eninc(nen)=enincF
              numinc=nen
              enincmin=eninc(1)
              enincmax=eninc(nen)
            endif
          endif
        endif
c
c Remove incident energies above energy given by Estop
c
        nen=numinc
        do 340 i=1,nen
          if (eninc(i).gt.Estop) then
            numinc=i-1
            goto 300
          endif
  340   continue
      else
c
c 2. Population distribution as the initial state
c
c npopE   : number of energies for population distribution
c npopJ   : number of spins for population distribution
c npopP   : number of parities for population distribution
c pbeg    : help variable
c numbins : maximal number of continuum excitation energy bins
c
        inquire (file=energyfile,exist=lexist)
        if (.not.lexist) then
          write(*,'(" TALYS-error: if projectile 0, specify a range",
     +        " of excitation energies in a file ",a73)') energyfile
          stop
        endif
        open (unit=2,file=energyfile,status='old')
        read(2,*,end=510,err=510) npopE,npopJ,npopP
        if (npopE.lt.2.or.npopE.gt.numpop) then
          write(*,'(" TALYS-error: 2 <= bins <=",i4," in population ",
     +      "distribution file")') numpop
          stop
        endif
        if (npopJ.lt.0.or.npopJ.gt.numJ+1) then
          write(*,'(" TALYS-error: 0 <= number of spins <=",i3,
     +      " + 1 in population distribution file")') numJ
          stop
        endif
        if (npopJ.gt.0.and.(npopP.lt.1.or.npopP.gt.2)) then
          write(*,'(" TALYS-error: 1 <= number of parities <= 2",
     +      " in population distribution file")')
          stop
        endif
c
c Only excitation energy distribution (no spins)
c
        exBin = INT(exEne * 10.0 + 0.5) + 1
        print *,exBin, exSpin, exParity
        
        if (npopJ.eq.0) then
          do 350 nen=1,npopE
            read(2,*,end=510,err=510) EdistE(nen),PdistE(nen)
  350     continue
          PdistE(exBin) = 1.0
        else
c
c Spin-dependent excitation energy distribution (no total)
c
          if (npopP.eq.1) then
            pbeg=1
          else
            pbeg=-1
          endif
          do 360 parity=pbeg,1,2
            do 370 nen=1,npopE
                read(2,*,end=510,err=510) EdistE(nen),
     +            (PdistJP(nen,J,parity),J=0,npopJ-1)
  370       continue
  360     continue
          do 375 nen=1,npopE
            if (EdistE(nen).le.EdistE(nen-1)) then
              if (EdistE(1).eq.0.) goto 375
              write(*,'(" TALYS-error: excitation energies must",
     +          " be given in ascending order, or the number",
     +          " of population bins is not correct")')
              stop
            endif
  375     continue
          if (npopP.eq.1) then
            do 380 nen=1,npopE
              do 380 J=0,npopJ-1
                PdistJP(nen,J,-1)=PdistJP(nen,J,1)
  380       continue
          endif
        endif
        close (unit=2)
        numinc=1
        eninc(1)=EdistE(npopE)
        enincmin=eninc(1)
        enincmax=eninc(1)

c Add the Bin that is non-zero manually based on exEne given from c++
        PdistJP(exBin,exSpin,exParity) = 1.0
      endif
c
c In case of built-in energy range, write an explicit 'energies' file
c
  300 if (enincF.gt.0..or.fexist) then
        open (unit=2,file='energies',status='replace')
        do 385 nen=1,numinc
          write(2,'(1p,g12.4)') eninc(nen)
  385   continue
        close (unit=2)
      endif
c
c If requested by input: retrieve best set of adjusted input parameters
c
c afile   : TALYS file with mass number
c bestchar: help variable
c inum    : counter
c ipath   : counter
c lenbest : length of best file
c bestfile: adjusted "best" parameter file
c convert : subroutine to convert input line from upper case to lowercase
c
      if (flagbest) then
        afile='000.talys'
        write(afile(1:3),'(i3.3)') Atarget
        bestchar=ptype0//'-'//trim(nuc(Ztarget))//afile
        if (bestpath(1:1).eq.' ')
     +    bestpath='best/                                   '
        do 390 i=1,40
          if (bestpath(i:i).eq.' ') then
            if (bestpath(i-1:i-1).ne.'/') then
              bestpath(i:i)='/'
              lenbest=i
            else
              lenbest=i-1
            endif
            goto 400
          endif
  390   continue
  400   if (Starget(2:2).eq.' ') then
          if (Ltarget.eq.0) then
            write(bestpath(lenbest+1:lenbest+5),'(a1,i3.3,"/")')
     +        Starget(1:1),Atarget
            write(bestpath(lenbest+6:lenbest+20),'(a15)') bestchar
          else
            write(bestpath(lenbest+1:lenbest+6),'(a1,i3.3,"m/")')
     +        Starget(1:1),Atarget
            write(bestpath(lenbest+7:lenbest+21),'(a15)') bestchar
          endif
        else
          if (Ltarget.eq.0) then
            write(bestpath(lenbest+1:lenbest+6),'(a2,i3.3,"/")')
     +        Starget(1:2),Atarget
            write(bestpath(lenbest+7:lenbest+21),'(a15)') bestchar
          else
            write(bestpath(lenbest+1:lenbest+7),'(a2,i3.3,"m/")')
     +        Starget(1:2),Atarget
            write(bestpath(lenbest+8:lenbest+22),'(a15)') bestchar
          endif
        endif
        bestfile=trim(path)//bestpath
        inquire (file=bestfile,exist=lexist)
        if (.not.lexist) then
          write(*,'(" TALYS-warning: best file does not exist: ",a)')
     +      trim(bestfile)
          return
        endif
        open (unit=3,file=bestfile,status='old')
        inum=0
  410   read(3,'(a80)',end=420) key
        inum=inum+1
        i=numlines-inum
        inline(i)=key
        call convert(i)
        goto 410
  420   close (unit=3)
        if (inum.gt.0) then
          if (flagbestend) then
            do 430 i=1,inum
              inline(nlines+i)=inline(numlines-i)
  430       continue
          else
            do 440 i=nlines,1,-1
              inline(i+inum)=inline(i)
  440       continue
            do 450 i=1,inum
              inline(i)=inline(numlines-i)
  450       continue
          endif
          nlines=nlines+inum
        endif
      endif
      return
  500 write(*,'(" TALYS-error: Wrong input: ",a80)') inline(i)
      stop
  510 write(*,'(" TALYS-error: Problem in file ",a73)') energyfile
      write(*,'(" after E= ",es12.5)') Ein
      stop
      end
Copyright (C)  2016 A.J. Koning, S. Hilaire and S. Goriely
