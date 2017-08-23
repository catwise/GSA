      subroutine Setup(IUnitIn, IUnitOut, ColNam, ColTyp, TblFil,
     +                 RAstr, Decstr, IA1, IA2, ID1, ID2,
     +                 RAmax, RAmin, Decmax, Decmin, nColTyp,
     +                 Field, IFa, IFb, NF, ColTyp3, ColTyp4,
     +                 NLines, LRecL, FilNam, OK)
c-----------------------------------------------------------------------
c
c  Set up an input file for gsa; read file on IUnitIn, write a direct-
c  access file on IUnitOut; strip headers from table files after finding
c  and saving the column names and column types and finding the RA and
c  Dec column numbers; return the logical record length
c
c-----------------------------------------------------------------------
                     Integer*4  MaxFld
                     Parameter (MaxFld = 1000)
c
      character*5000 ColNam, ColTyp, Line, ColTyp3, ColTyp4
      character*150  FilNam
      character*25   RAstr, Decstr, Field(MaxFld)
      character*11   FmtLine
      character*1    Chr
      real*4         RAmax, RAmin, Decmax, Decmin, RTmp
      integer*4      IUnitIn, IUnitOut, NLines, LRecL, IA1, IA2, ID1,
     +               ID2, LNBlnk, MaxLen, NTblHdr, N, IFa(MaxFld),
     +               IFb(MaxFld), NF, nColTyp
      Logical        TblFil, OK, GotCN, GotCT
      Byte           IChr
      Equivalence   (Chr, IChr)
c
c-----------------------------------------------------------------------
c
      OK = .True.
      GotCN = .False.
      GotCT = .False.
      NTblHdr = 0
      MaxLen = 0
      NLines = 0
      nColTyp = 0
c
      if (TblFil) then
10      read (IUnitIn, '(A5000)', end = 3000, err = 3001) Line
        NTblHdr = NTblHdr + 1
c
        if (Line(1:1) .eq. '|') then
          nColTyp = nColTyp + 1
          if (.not.GotCN) then
            ColNam = Line
            GotCN = .True.
            if (ColNam(LNBlnk(ColNam):LNBlnk(ColNam)) .ne. '|') then
              print *,
     +         'WARNING: table-file header does not terminate with "|"'
              print *,ColNam(1:LNBlnk(ColNam))
              N = LNBlnk(ColNam) + 1
              if (N .gt. 5000) N = 5000
              ColNam(N:N) = '|'
              print *,'modified to the following:'
              print *,ColNam(1:LNBlnk(ColNam))
            end if
c
            NF = 0
            do 20 N = 1, LNBlnk(ColNam)-1
              if (ColNam(N:N) .eq. '|') NF = NF + 1
20          continue
            if (NF .gt. MaxFld) then
              print *,'WARNING: too many fields in table file: ',NF
              print *,'         max = ',MaxFld
              print *,'         file: ',FilNam(1:LNBlnk(FilNam))
              print *,'         processing only the max value'
              NF = MaxFld
            end if
c
            call GetFlds(ColNam,Field,IFa,IFb,NF)
c
            do 40 N = 1, NF
              if (Field(N) .eq. RAstr) then
                IA1 = IFa(N)
                IA2 = IFb(N)
                go to 50
              end if
40          continue
            print *,'ERROR: can''t find '//RAstr(1:LNBlnk(RAstr))
     +              //' in header line:'
            print *,Line(1:LNBlnk(Line))
            OK = .False.
            return
c
50          do 60 N = 1, NF
              if (Field(N) .eq. Decstr) then
                ID1 = IFa(N)
                ID2 = IFb(N)
                go to 100
              end if
60          continue
            print *,'ERROR: can''t find '//Decstr(1:LNBlnk(Decstr))
     +              //' in header line:'
            print *,Line(1:LNBlnk(Line))
            OK = .False.
            return
100         continue
c
          else if (.not.GotCT) then
            ColTyp = Line
            GotCT = .True.
          else if (nColTyp .eq. 3) then
            ColTyp3 = Line
          else if (nColTyp .eq. 4) then
            ColTyp4 = Line
          end if
        else
          Chr = Line(1:1)              ! check for "\"
          if (IChr .ne. 92) then
            NLines = 1
            MaxLen = LNBlnk(Line)
            NTblHdr = NTblHdr - 1
            go to 300
          end if
        end if
        go to 10
      end if
c
c-----------------------------------------------------------------------
c
300   read(IUnitIn, '(A5000)', end = 400, err = 3002) Line
      NLines = NLines + 1
      if (LNBlnk(Line) .gt. MaxLen) MaxLen = LNBlnk(Line)
      go to 300
c
400   rewind(IUnitIn)
      if (LNBlnk(Line) .eq. 0) NLines = NLines - 1  ! clip any blank
      if (TblFil) then                              !  leftovers
        do 410 N = 1, NTblHdr
          read(IUnitIn, '(A5000)', end = 3003, err = 3004) Line
410     continue
      end if
c
      LRecL = MaxLen
      if ((LRecL .lt. IA2) .or. (LRecL .lt. ID2) .or. (LRecL .lt. 2))
     + go to 3009
      call MakeFmtD(FmtLine,MaxLen)
      RAmax =  -9999.9
      RAmin =   9999.9
      Decmax = -9999.9
      Decmin =  9999.9
      do 500 N = 1, NLines
         read(IUnitIn, '(A5000)', end = 3005, err = 3006) Line
         write(IUnitOut,FmtLine) Line(1:MaxLen)
         read (Line(IA1:IA2), *, err = 3007) RTmp
         if (RTmp .lt. RAmin) RAmin = RTmp
         if (RTmp .gt. RAmax) RAmax = RTmp
         read (Line(ID1:ID2), *, err = 3007) RTmp
         if (RTmp .lt. Decmin) Decmin = RTmp
         if (RTmp .gt. Decmax) Decmax = RTmp
500   continue
      return
c
c-----------------------------------------------------------------------
c
3000  print *,'ERROR: unexpected EoF while reading header of table file'
      print *,'       ',FilNam(1:LNBlnk(FilNam))
      OK = .False.
      return
c
3001  print *,'ERROR: read error while reading header of table file'
      print *,'       ',FilNam(1:LNBlnk(FilNam))
      OK = .False.
      return
c
3002  print *,'ERROR: read error while reading data in file'
      print *,'       ',FilNam(1:LNBlnk(FilNam))
      OK = .False.
      return
c
3003  print *,'ERROR: unexpected EoF while reading header of table file'
      print *,'       ',FilNam(1:LNBlnk(FilNam))
      print *,'This error occurred on the second pass through the file'
      OK = .False.
      return
c
3004  print *,'ERROR: read error while reading header of table file'
      print *,'       ',FilNam(1:LNBlnk(FilNam))
      print *,'This error occurred on the second pass through the file'
      OK = .False.
      return
c
3005  print *,'ERROR: unexpected EoF while reading data in file'
      print *,'       ',FilNam(1:LNBlnk(FilNam))
      print *,'This error occurred on the second pass through the file'
      OK = .False.
      return
c
3006  print *,'ERROR: read error while data in file'
      print *,'       ',FilNam(1:LNBlnk(FilNam))
      print *,'This error occurred on the second pass through the file'
      OK = .False.
      return
c
3007  print *,'ERROR: read error on RA value '//Line(IA1:IA2)
     +       //' in file'
      print *,'       ',FilNam(1:LNBlnk(FilNam))
      print *,'This error occurred on data line no. ', N
      OK = .False.
      return
c
3008  print *,'ERROR: read error on Dec value '//Line(IA1:IA2)
     +       //' in file'
      print *,'       ',FilNam(1:LNBlnk(FilNam))
      print *,'This error occurred on data line no. ', N
      OK = .False.
      return
c
3009  print *,
     + 'ERROR: data line length must be at least as long as the last'
      print *,
     + '       column in the RA field and the Dec field; found a length'
      print *,
     + '       of ',LRecL,' in file ',FilNam(1:LNBlnk(FilNam))
      print *,
     + '       last column in the RA field: ',IA2
      print *,
     + '       last column in the Dec field: ',ID2
      OK = .False.
      return
c
      end
c
c=======================================================================
c
      Subroutine MakeFmtD(FmtLine,L)
c
      Character*11 FmtLine
      Character*4  Flen
      Integer*4    L
c
c-----------------------------------------------------------------------
c
      write(Flen,'(I4)') L
      FmtLine = '(A'//Flen//'$)'
      return
      end
