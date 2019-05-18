c-----------------------------------------------------------------------
c       GSA  -  General Source Association
c
c  version 1.0  A60721: initial version
c          1.1  A60728: added discrepancy statistics and clipping of
c                       any blank partial lines at end of input files
c                       (saw this in some SWIRE data, e.g., irac2.inp)
c          1.2  A60731: added discrepancy statistics for unit vectors;
c          1.3    "     optional offsets (file2-file1) to be subtracted
c                       from file-2 source position parameters
c          1.4  A60802: added option to output index list with or
c                       without full merged records
c          1.5  A60815: added fdate usage
c          1.6  A60818: added RA2 & Dec2 to null1 output in RA1 & Dec1
c                       location for use by mgsa
c          1.7  A61107: added option to merge two bandmerged files via
c                       common srcid# values in a specified band
c          1.8  A70410: added look-ahead for possible better match when
c                       keeping single best match only; added -bh option
c          2.0  A80720: ported to WISE linux g95; deleted WinMac tags
c          2.1  A90414: fixed MAXFLD & ColNam inconsistency
c          2.2  A90714: added retention of all input "pipe" header
c                       lines; set "dist" to -1 or -2 to indicate first
c                       or second record null in unmatched source output
c          2.3  A90924: added statistics output table file ("-o2")
c          2.4  A90929: added #Assns to statistics output table file
c          2.5  A91001: fixed "pipe" tags for dist
c          2.6  A91113: fixed indexing error for single-source files
c          2.7  B60715: installed on-demand-randomized temporary work
c                       files so that more than one instance of gsa can
c                       run in the same working directory
c          2.8  B80125: installed CatWISE-specific requirement to match
c                       mdetID in addition to meet radial-distance rqmt
c          2.8  B80128: installed option to skip nulls in match fields,
c                       i.e., treat as not matched
c          2.8  B80419: added spec for temporary work directory name
c          2.8  B81002: added "-nm2" option: no multi-matches in file 2
c          2.9  B81009: added OutFNam into scratch file names to reduce
c                       probability of name collisions
c          2.91 B81012: added debug trace for specific mdetID values
c                       (see gsa-dbg261.f and best2cm-dbg261.f)
c          2.92 B81014: added "Best2Come" processing in mdetID mode
c          2.93 B81220: force Dec order when abs(DeltaRA) > 350, i.e.,
c                       indicating an RA 0-360 crossing, which makes
c                       RA order incorrectly assume no matches across
c                       RA boundary, hence need Dec order.
c          2.94 B81224: force Dec order when CatWISE = T; also extend
c                       search from 1.1*Dist to 1.5*Dist; change all
c                       reals to real*8 except in random-number gen.
c          2.95 B81224: restored forcing Dec order when CatWISE = T;
c                       added prints on RA/Dec min/max
c          2.96 B90509: added -ar to allow RA order
c
c NOTE: when redelivering the linux version, be sure to set the path
c       switch from "\" to "/" (search on "linux" to find the code)
c
c-----------------------------------------------------------------------
c
                     Integer*4  MaxFld, MaxBand
                     Parameter (MaxFld = 1000, MaxBand = 7)
c
      Character*5000 Line, ColNam(3), ColTyp(2), Line1, Null1, Null2,
     +               ColTyp3(2), ColTyp4(2)
      Character*500  InFNam1, InFNam2, OutFNam, FilNam, Rej1Nam,
     +               Rej2Nam, StatNam, TempDir
      Character*500  TmpNam(3), BmgNam
      Character*64   OutStr, TmpStr
      Character*25   NumStr, RAstr1, RAstr2, Decstr1, Decstr2,
     +               Field1(MaxFld), Field2(MaxFld), NamDes, NamGot
      Character*11   Vsn, FmtLine, FmtLin1, FmtLin2
      Character*9    DistStr
      Character*8    CDate, CTime
      Character*7    SrcID
      Character*5    Flag
      Character*1    Chr(5000)
      Real*8         RAmax(2), RAmin(2), Decmax(2), Decmin(2),
     +               AvgDec, DeltRA, DeltDec, Window, RTmp1, RTmp2,
     +               RTmp3, RTmp4, RTmp5, RTmp6, ran1
      Real*8         d2r, Alpha, Delta, cosD, Ang, AngMin,
     +               Cross, Dist, SumDist, SumSqDist, SumDec, SumSqDec,
     +               SumTARA, SumSqTARA, Alpha1, Delta1, DelTARA,
     +               DistMax, DelTARAmax, DeltDecMax, ADist,
     +               SumDX, SumDY, SumDZ, SumSqDX, SumSqDY, SumSqDZ,
     +               DXmax, DYmax, DZmax, dRA, dDec, dX, dY, dZ
      Integer*4      NArgs, NArg, IArgC, IA1(2), IA2(2), ID1(2), N,
     +               ID2(2), Access, LNBlnk, LRecL(2), NLines(2), I, J,
     +               NAssns, L, M, N1, N2, N2best, K, N3, Len1, IH2,
     +               N1single, N2single, NTotal, IF1a(MaxFld), LenOut,
     +               IF1b(MaxFld), IF2a(MaxFld), IF2b(MaxFld), NF1, NF2,
     +               N1multi, N2multi, Izero, IDBcoll8, IS1(2), IS2(2),
     +               IC1(2), IC2(2), IPa1, IPa2, IPb1, IPb2, NB1, NB2,
     +               IDB2(MaxBand), ISr1(MaxBand), ISr2(MaxBand), IF1,
     +               IF2, ISys, System, NBands, nColTyp(2), Nct,
     +               mdetID1, mdetID2, IM1(2), IM2(2)
c
      Real*8, allocatable :: X1(:), Y1(:), Z1(:), X2(:), Y2(:), Z2(:)
      Real*8, allocatable :: Coord1(:), Coord2(:)
      Integer*4, allocatable :: Indx1(:), Indx2(:)
      Byte, allocatable ::   Mstate(:)
      Logical        NeedHelp, GotIn1, GotIn2, GotOut, GotA1, GotA2,
     +               GotD1, GotD2, TblFil(2), GotRD, AllSrc1, AllAssn,
     +               OK, DecOrder, GotMatch, AllSrc2, Reject1, Reject2,
     +               MergOut, IndxOut, BmgColl8, Best2Come, BmgHed,
     +               GotStat, CatWISE, SkipErr, nm2, ForceDec
c
      Equivalence   (Line, Chr(1))
c
      Common / VDT / CDate, CTime, Vsn
c
      data Vsn/'2.96 B90509'/, NeedHelp/.False./, GotIn1/.False./,
     +     GotIn2/.False./, GotA1/.False./, GotA2/.False./, IH2/60/,
     +     GotOut/.False./, GotD1/.False./, GotD2/.False./, NLines/2*0/,
     +     GotRD/.False./,AllSrc1/.True./, AllAssn/.True./, N1single/0/,
     +     TmpNam/'gsa%%%tmp1','gsa%%%tmp2','gsa%%%tmp3'/, IS2/2*0/,
     +     d2r/1.74532925199433d-2/, NTotal/0/, SumDist/0.0d0/, Len1/0/,
     +     SumSqDist/0.0d0/, SumDec/0.0d0/, SumSqDec/0.0d0/, IC2/2*0/,
     +     SumTARA/0.0d0/, SumSqTARA/0.0d0/, N1multi/0/, N2multi/0/,
     +     DistMax/0.0d0/, DelTARAmax/0.0d0/, DeltDecMax/0.0d0/,
     +     SumDX/0.0d0/, SumDY/0.0d0/, SumDZ/0.0d0/, SumSqDX/0.0d0/,
     +     SumSqDY/0.0d0/, SumSqDZ/0.0d0/, DXmax/0.0d0/, DYmax/0.0d0/,
     +     DZmax/0.0d0/, dRA/0.0d0/, dDec/0.0d0/, dX/0.0d0/, dY/0.0d0/,
     +     dZ/0.0d0/, AllSrc2/.True./, Reject1/.False./, Izero/0/,
     +     Reject2/.False./, MergOut/.True./, IndxOut/.False./, Nct/0/,
     +     BmgColl8/.False./, BmgNam/'gsa%%%%bmg'/, IS1/-9,-9/,
     +     IC1/-9,-9/, NAssns/0/, N2single/0/, BmgHed/.False./,
     +     IF2/0/, IPa1/0/, IPa2/0/, IPb1/0/, IPb2/0/, NB1/1/,
     +     GotStat/.false./, CatWISE/.false./, IM1,Im2/4*-9/,
     +     SkipErr/.false./, TempDir/'./'/, nm2/.false./,
     +     ForceDec/.true./
c
c-----------------------------------------------------------------------
c
      NArgs = IArgC()
1     If ((NArgs .lt. 8) .or. NeedHelp) then
        print *,'gsa: General Source Association vsn ', Vsn
        print *
        print *,'usage:  gsa <flags> <specifications>'
        print *
        print *,'where <flags> <specifications> must be: '
        print *
        print *,'    -t intblnam     (name of input table file)'
        print *,'    -s intxtnam     (name of input simple text file)'
        print *,'    -o outfilnam    (name of output association file)'
        print *,
     +  '    -ra1 RAspec1    (RA specification for first input file)'
        print *,
     +  '    -ra2 RAspec2    (RA specification for second input file)'
        print *,
     +  '    -dec1 Decspec1  (Dec specification for first input file)'
        print *,
     +  '    -dec2 Decspec2  (Dec specification for second input file)'
        print *,
     +  '    -r AssnDist     (radial association distance, arcsec)'
        print *,
     +  '    -a1             (optional; output only closest match)'
        print *,
     +  '    -aa             (optional; output all matches; default)'
        print *,
     +'    -ns             (optional; don''t include unmatched sources'
        print *,
     +'                     in merged output)'
        print *,
     +'    -ns1            (optional; don''t include unmatched sources'
        print *,
     +'                     from primary file in merged output)'
        print *,
     +'    -ns2            (optional; don''t include unmatched sources'
        print *,
     +'                     from secondary file in merged output)'
        print *,
     +'    -as             (optional; output all sources, matched or'
        print *,
     +'                     not; default)'
        print *,
     +  '    -o2  statsnam   (optional; name of table file containing'
        print *,
     +'                     statistical summary information)'
        print *,
     +  '    -rf1 reject1nam (optional; name of file for unmatched'
        print *,
     +'                     sources from primary file)'
        print *,
     +  '    -rf2 reject2nam (optional; name of file for unmatched'
        print *,
     +'                     sources from secondary file)'
        print *,
     +  '    -dRA RAoffset   (optional; RA offset in arcsec to be'
        print *,
     +  '                     subtracted from secondary-file RA)'
        print *,
     +  '    -dDec Decoffset (optional; Dec offset in arcsec to be'
        print *,
     +  '                     subtracted from secondary-file Dec)'
        print *,
     +  '    -dX Xoffset     (optional; unit vector X offset to be'
        print *,
     +  '                     subtracted from secondary-file vector)'
        print *,
     +  '    -dY Yoffset     (optional; unit vector Y offset to be'
        print *,
     +  '                     subtracted from secondary-file vector)'
        print *,
     +  '    -dZ Zoffset     (optional; unit vector Z offset to be'
        print *,
     +  '                     subtracted from secondary-file vector)'
        print *,
     +  '    -x indexopt     (optional; indexopt = 1 --> generate an'
        print *,
     +'                     index file of match record numbers instead'
        print *,
     +'                     of merged output; indexopt = 2 --> generate'
        print *,
     +'                     both; the "-o" option is still required;'
        print *,
     +'                     the index file name will have "_index"'
        print *,
     +'                     appended)'
        print *,
     +'    -se             (optional; skip read errors on match fields;'
        print *,
     +'                     treat as unmatchable source)'
        print *,
     +'    -td  tmpdir     (optional; temporary scratch directory name)'
        print *,
     +'    -ar             (optional; allow RA search order if faster'
        print *
        print *,
     +  'For general source association, the first eight specifications'
        print *,
     +  'are required (i.e., two input files must be specified, either'
        print *,
     +  '"-t" must be used twice, or "-s" must be used twice, or each'
        print *,
     +  'must be used once), and the "-b" option must not be used. For'
        print *,
     +  'table files, RA and Dec specifications are column names; for'
        print *,
     + 'simple text files, they are column ranges (e.g., 11-20) with no'
        print *,
     +  'embedded blanks.'
        print *
        print *,
     + 'Note that because gsa uses direct-access disk I/O for rapid'
        print *,
     +'random access, data rows in both types of file must all be the'
        print *,
     +'same length within their own file; a line is considered to end'
        print *,
     +'at the last non-blank character.'
        print *
        Print *,'(Spitzer-specific)'
        print *
        print *,
     +  '    -b band#        (band number for collating two bandmerged'
        print *,
     +  '                     files via srcid in the specified band)'
        print *,
     +  '    -bh #bands      (optional; include bandmerge header'
        print *,
     +  '                     keywords NBands and Total_Merged_Number)'
        print *
        print *,
     + 'To collate two bandmerged files, the first three specifications'
        print *,
     +'are required, "-b" must specify which band the two files have in'
        print *,
     + 'common, and all other options are ignored '
     +    //'except for "-as", "-x",'
        print *,
     + '"-ns", "-ns1", "-ns2", "-rf1", and "-rf2", which are optional.'
        print *,
     + 'Collating bandmerged files is a Spitzer-specific process.'
        print *
        print *,
     + 'To place the bandmerge header keywords Total_Merged_Number and'
        print *,
     +'NBands in the output header (assuming table file format), the'
        print *,
     +'"-bh" option should be used to specify the value of NBands. This'
        print *,
     +'is intended to allow non-bandmerged files to be processed by'
        print *,
     +'certain utilities; it should not be used with the "-b" option.'
        print *
        Print *,'(CatWISE-specific)'
        print *
        print *,
     +  '    -cw             (optional; CatWISE-specific processing)'
        print *,
     +  '    -nm2            (optional; no multi-matches in file #2)'
        print *
        print *,
     +  'For CatWISE-specific processing ("-cw"), the first eight'
        print *,
     +  'specifications are required; the two mdex files must be'
        print *,
     +  'specified with the "-t" flag; in addition to the usual radial-'
        print *,
     +  'distance requirement, the "mdetID" values must match.'
        print *,
     +  'If "-nm2" is specified, sources in the second input file may'
        print *,
     +  'match only the closest source in the first input file.'
        call exit(32)
      end if
c
      call signon('gsa')
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      NArg = 0
4     NArg = NArg + 1
      call GetArg(NArg,Flag)
c
      If (Flag .eq. '-t') then
        NArg = NArg + 1
        if (NArg .le. NArgs) then
          call GetArg(NArg,FilNam)
        else
          go to 3010
        end if
        if (GotIn2) then
          print *,'ERROR: more than two input files specified'
          NeedHelp = .True.
          print *
          go to 1
        else if (GotIn1) then
          InFNam2 = FilNam
          GotIn2 = .True.
          TblFil(2) = .True.
        else
          InFNam1 = FilNam
          GotIn1 = .True.
          TblFil(1) = .True.
        end if
        if (Access(FilNam(1:LNBlnk(FilNam)),' ') .ne. 0) then
          print *,'File not found: ',FilNam(1:LNBlnk(FilNam))
          call exit(64)
        end if
      else if (Flag .eq. '-s') then
        NArg = NArg + 1
        if (NArg .le. NArgs) then
          call GetArg(NArg,FilNam)
        else
          go to 3010
        end if
        if (GotIn2) then
          print *,'ERROR: more than two input files specified'
          NeedHelp = .True.
          print *
          go to 1
        else if (GotIn1) then
          InFNam2 = FilNam
          GotIn2 = .True.
          TblFil(2) = .False.
        else
          InFNam1 = FilNam
          GotIn1 = .True.
          TblFil(1) = .False.
        end if
        if (Access(FilNam(1:LNBlnk(FilNam)),' ') .ne. 0) then
          print *,'File not found: ',FilNam(1:LNBlnk(FilNam))
          call exit(64)
        end if
c
      Else If (Flag .eq. '-o') then
        NArg = NArg + 1
        if (NArg .le. NArgs) then
          call GetArg(NArg,OutFNam)
        else
          go to 3010
        end if
        GotOut = .True.
c
      Else If (Flag .eq. '-o2') then
        NArg = NArg + 1
        if (NArg .le. NArgs) then
          call GetArg(NArg,StatNam)
        else
          go to 3010
        end if
        GotStat = .True.
c
      Else If (Flag .eq. '-rf1') then
        NArg = NArg + 1
        if (NArg .le. NArgs) then
          call GetArg(NArg,Rej1Nam)
        else
          go to 3010
        end if
        Reject1 = .True.
c
      Else If (Flag .eq. '-rf2') then
        NArg = NArg + 1
        if (NArg .le. NArgs) then
          call GetArg(NArg,Rej2Nam)
        else
          go to 3010
        end if
        Reject2 = .True.
c
      Else If (Flag .eq. '-r') then
        NArg = NArg + 1
        if (NArg .le. NArgs) then
          call GetArg(NArg,NumStr)
        else
          go to 3010
        end if
        Read (NumStr,*,Err=5) Dist
        GotRD = .True.
        Go to 80
5       print *,'ERROR: illegal -r specification: ', NumStr
        NeedHelp = .True.
        Print *
        Go to 1
c
      Else If (Flag .eq. '-b') then
        NArg = NArg + 1
        if (NArg .le. NArgs) then
          call GetArg(NArg,NumStr)
        else
          go to 3010
        end if
        Read (NumStr,*,Err=6) IDBcoll8
        if ((IDBcoll8 .lt. 1) .or. (IDBcoll8 .gt. 99)) go to 6
        BmgColl8 = .True.
        BmgHed = .False.
        Go to 80
6       print *,'ERROR: illegal -b specification: ', NumStr
        NeedHelp = .True.
        Print *
        Go to 1
c
      Else If (Flag .eq. '-ra1') then
        NArg = NArg + 1
        if (NArg .le. NArgs) then
          call GetArg(NArg,NumStr)
        else
          go to 3010
        end if
        GotA1 = .True.
        if (TblFil(1)) then
          RAstr1 = NumStr
        else
          Read (NumStr(1:Index(NumStr,'-')-1),*,Err=10) IA1(1)
          Read (NumStr(Index(NumStr,'-')+1:LNBlnk(NumStr)),*,Err=10)
     +          IA2(1)
          Go to 80
10        print *,'ERROR: illegal -ra1 specification: ', NumStr
          NeedHelp = .True.
          Print *
          Go to 1
        end if
c
      Else If (Flag .eq. '-ra2') then
        NArg = NArg + 1
        if (NArg .le. NArgs) then
          call GetArg(NArg,NumStr)
        else
          go to 3010
        end if
        GotA2 = .True.
        if (TblFil(2)) then
          RAstr2 = NumStr
        else
          Read (NumStr(1:Index(NumStr,'-')-1),*,Err=20) IA1(2)
          Read (NumStr(Index(NumStr,'-')+1:LNBlnk(NumStr)),*,Err=20)
     +          IA2(2)
          Go to 80
20        print *,'ERROR: illegal -ra2 specification: ', NumStr
          NeedHelp = .True.
          Print *
          Go to 1
        end if
c
      Else If (Flag .eq. '-dec1') then
        NArg = NArg + 1
        if (NArg .le. NArgs) then
          call GetArg(NArg,NumStr)
        else
          go to 3010
        end if
        GotD1 = .True.
        if (TblFil(1)) then
          Decstr1 = NumStr
        else
          Read (NumStr(1:Index(NumStr,'-')-1),*,Err=30) ID1(1)
          Read (NumStr(Index(NumStr,'-')+1:LNBlnk(NumStr)),*,Err=30)
     +          ID2(1)
          Go to 80
30        print *,'ERROR: illegal -dec1 specification: ', NumStr
          NeedHelp = .True.
          Print *
          Go to 1
        end if
c
      Else If (Flag .eq. '-dec2') then
        NArg = NArg + 1
        if (NArg .le. NArgs) then
          call GetArg(NArg,NumStr)
        else
          go to 3010
        end if
        GotD2 = .True.
        if (TblFil(2)) then
          Decstr2 = NumStr
        else
          Read (NumStr(1:Index(NumStr,'-')-1),*,Err=40) ID1(2)
          Read (NumStr(Index(NumStr,'-')+1:LNBlnk(NumStr)),*,Err=40)
     +          ID2(2)
          Go to 80
40        print *,'ERROR: illegal -dec2 specification: ', NumStr
          NeedHelp = .True.
          Print *
          Go to 1
        end if
c
      Else If (Flag .eq. '-bh') then
        NArg = NArg + 1
        if (NArg .le. NArgs) then
          call GetArg(NArg,NumStr)
        else
          go to 3010
        end if
        Read (NumStr,*,Err=46) NBands
        if ((NBands .lt. 2) .or. (NBands .gt. 99)) go to 46
        BmgColl8 = .False.
        BmgHed = .True.
        Go to 80
46      print *,'ERROR: illegal -bh specification: ', NumStr
        NeedHelp = .True.
        Print *
        Go to 1
c
      Else If (Flag .eq. '-a1') then
        AllAssn = .False.
c
      Else If (Flag .eq. '-aa') then
        AllAssn = .True.
c
      Else If (Flag .eq. '-ns') then
        AllSrc1 = .False.
        AllSrc2 = .False.
c
      Else If (Flag .eq. '-ns1') then
        AllSrc1 = .False.
c
      Else If (Flag .eq. '-ns2') then
        AllSrc2 = .False.
c
      Else If (Flag .eq. '-as') then
        AllSrc1 = .True.
        AllSrc2 = .True.
c
      Else If (Flag .eq. '-dRA') then
        NArg = NArg + 1
        if (NArg .le. NArgs) then
          call GetArg(NArg,NumStr)
        else
          go to 3010
        end if
        Read (NumStr,*,Err=50) dRA
        dRA = dRA/3600.0d0
        Go to 80
50      print *,'ERROR: illegal -dRA specification: ', NumStr
        NeedHelp = .True.
        Print *
        Go to 1
c
      Else If (Flag .eq. '-dDec') then
        NArg = NArg + 1
        if (NArg .le. NArgs) then
          call GetArg(NArg,NumStr)
        else
          go to 3010
        end if
        Read (NumStr,*,Err=55) dDec
        dDec = dDec/3600.0d0
        Go to 80
55      print *,'ERROR: illegal -dDec specification: ', NumStr
        NeedHelp = .True.
        Print *
        Go to 1
c
      Else If (Flag .eq. '-dX') then
        NArg = NArg + 1
        if (NArg .le. NArgs) then
          call GetArg(NArg,NumStr)
        else
          go to 3010
        end if
        Read (NumStr,*,Err=60) dX
        Go to 80
60      print *,'ERROR: illegal -dX specification: ', NumStr
        NeedHelp = .True.
        Print *
        Go to 1
c
      Else If (Flag .eq. '-dY') then
        NArg = NArg + 1
        if (NArg .le. NArgs) then
          call GetArg(NArg,NumStr)
        else
          go to 3010
        end if
        Read (NumStr,*,Err=65) dY
        Go to 80
65      print *,'ERROR: illegal -dY specification: ', NumStr
        NeedHelp = .True.
        Print *
        Go to 1
c
      Else If (Flag .eq. '-dZ') then
        NArg = NArg + 1
        if (NArg .le. NArgs) then
          call GetArg(NArg,NumStr)
        else
          go to 3010
        end if
        Read (NumStr,*,Err=70) dZ
        Go to 80
70      print *,'ERROR: illegal -dZ specification: ', NumStr
        NeedHelp = .True.
        Print *
        Go to 1
c
      Else If (Flag .eq. '-x') then
        NArg = NArg + 1
        if (NArg .le. NArgs) then
          call GetArg(NArg,NumStr)
        else
          go to 3010
        end if
        if (NumStr .eq. '1') then
          MergOut = .False.
        else if (NumStr .ne. '2') then
          print *,'ERROR: illegal -x specification: ', NumStr
          NeedHelp = .True.
          Print *
          Go to 1
        end if
        IndxOut = .True.
c
      Else If (Flag .eq. '-cw') then
        CatWISE = .true.
c
      Else If (Flag .eq. '-nm2') then
        nm2 = .true.
c
      Else If (Flag .eq. '-se') then
        SkipErr = .true.
c
      Else If (Flag .eq. '-ar') then
        ForceDec = .false.
c
      Else if (Flag .eq. '-td') then
        NArg = NArg + 1
        call GetArg(NArg,TempDir)
        k = lnblnk(TempDir)
        if (TempDir(k:k) .ne. '/') TempDir(k+1:k+1) = '/'
c
      Else
        print *,'Unrecognized flag: ',Flag
        NeedHelp = .True.
        Print *
        Go to 1
      End If
c
80    If (NArg .lt. NArgs) Go to 4
c
      If (BmgColl8) then
        If (.not.(GotIn1 .and. GotIn2 .and. GotOut))
     +  then
         If (.not.GotIn1) print *,'ERROR: no input file specified'
         If (.not.GotIn2) print *,'ERROR: only one input file specified'
         If (.not.GotOut) print *,'ERROR: no output file specified'
         NeedHelp = .True.
         Print *
         Go to 1
        end if
        go to 90
      end if
c
      If (.not.(GotIn1 .and. GotIn2 .and. GotOut .and. GotA1 .and. GotA2
     +    .and. GotD1  .and. GotD2  .and. GotRD))
     + then
        If (.not.GotIn1) print *,'ERROR: no input file specified'
        If (.not.GotIn2) print *,'ERROR: only one input file specified'
        If (.not.GotOut) print *,'ERROR: no output file specified'
        If (.not.GotRD)  print *,'ERROR: no association distance given'
        If (.not.GotA1)  print *,'ERROR: "-ra1" specification missing'
        If (.not.GotA2)  print *,'ERROR: "-ra2" specification missing'
        If (.not.GotD1)  print *,'ERROR: "-dec1" specification missing'
        If (.not.GotD2)  print *,'ERROR: "-dec2" specification missing'
        NeedHelp = .True.
        Print *
        Go to 1
      end if
c
c-----------------------------------------------------------------------
c
90    Open (10, File = InFNam1)
	  call GetTmpNam(TmpNam(1),TempDir,OutFNam,OK)
      if (.not.OK) then
	    print *,'Unable to create temporary work file no. 1'
        call exit(64)
      end if
      Open (11, File = TmpNam(1))
      if (BmgColl8) then
        SrcID = 'srcid'
        if (IDBcoll8 .le. 9) then
          write (SrcID(6:6),'(I1)') IDBcoll8
        else
          write (SrcID(6:7),'(I2)') IDBcoll8
        end if
        RAstr1  = 'X'
        RAstr2  = 'X'
        Decstr1 = 'Y'
        Decstr2 = 'Y'
      end if
c
      call Setup(10, 11, ColNam(1), ColTyp(1), TblFil(1),
     +           RAstr1, Decstr1, IA1(1), IA2(1), ID1(1), ID2(1),
     +           RAmax(1), RAmin(1), Decmax(1), Decmin(1), nColTyp(1),
     +           Field1, IF1a, IF1b, NF1, ColTyp3(1), ColTyp4(1),
     +           NLines(1), LRecL(1), InFNam1, SkipErr, OK)
      close(10)
      close(11)
      if (.not.OK) then
	    print *,'Setup failed for input file no. 1'
        call DelTmp(TmpNam(1))
        call exit(64)
      end if
c
      if (BmgHed) then                 ! Check for "srcid", convert
        N = Index(ColNam(1),'srcid')   !  to "srcid1" if found
        if (N .lt. 2) go to 94
        if (ColNam(1)(N+5:N+5) .eq. ' ') then
          ColNam(1)(N+5:N+5) = '1'
        else if (ColNam(1)(N-1:N-1) .eq. ' ') then
          ColNam(1)(N-1:N+4) = 'srcid1'
        end if
      end if
c
94    Open (12, File = InFNam2)
	  call GetTmpNam(TmpNam(2),TempDir,OutFNam,OK)
      if (.not.OK) then
	    print *,'Unable to create temporary work file no. 2'
        call DelTmp(TmpNam(1))
        call exit(64)
      end if
      Open (13, File = TmpNam(2))
      call Setup(12, 13, ColNam(2), ColTyp(2), TblFil(2),
     +           RAstr2, Decstr2, IA1(2), IA2(2), ID1(2), ID2(2),
     +           RAmax(2), RAmin(2), Decmax(2), Decmin(2), nColTyp(2),
     +           Field2, IF2a, IF2b, NF2, ColTyp3(2), ColTyp4(2),
     +           NLines(2), LRecL(2), InFNam2, SkipErr, OK)
      ColNam(3) = ColNam(2)
      close(12)
      close(13)
      if (.not.OK) go to 3009
c
c-----------------------------------------------------------------------
c
      allocate(X1(NLines(1)))          !  allocate mem for X1
      if (.not.allocated(X1)) then
        print *,'ERROR: Unable to allocate memory for X1 array'
        go to 3009
      end if
      allocate(Y1(NLines(1)))          !  allocate mem for Y1
      if (.not.allocated(Y1)) then
        print *,'ERROR: Unable to allocate memory for Y1 array'
        go to 3009
      end if
      allocate(Z1(NLines(1)))          !  allocate mem for Z1
      if (.not.allocated(Z1)) then
        print *,'ERROR: Unable to allocate memory for Z1 array'
        go to 3009
      end if
      allocate(X2(NLines(2)))          !  allocate mem for X2
      if (.not.allocated(X2)) then
        print *,'ERROR: Unable to allocate memory for X2 array'
        go to 3009
      end if
      allocate(Y2(NLines(2)))          !  allocate mem for Y2
      if (.not.allocated(Y2)) then
        print *,'ERROR: Unable to allocate memory for Y2 array'
        go to 3009
      end if
      allocate(Z2(NLines(2)))          !  allocate mem for Z2
      if (.not.allocated(Z2)) then
        print *,'ERROR: Unable to allocate memory for Z2 array'
        go to 3009
      end if
      allocate(Indx1(NLines(1)))       !  allocate mem for Indx1
      if (.not.allocated(Indx1)) then
        print *,'ERROR: Unable to allocate memory for Indx1 array'
        go to 3009
      end if
      allocate(Indx2(NLines(2)))       !  allocate mem for Indx2
      if (.not.allocated(Indx2)) then
        print *,'ERROR: Unable to allocate memory for Indx2 array'
        go to 3009
      end if
      allocate(Coord1(NLines(1)))      !  allocate mem for Coord1
      if (.not.allocated(Coord1)) then
        print *,'ERROR: Unable to allocate memory for Coord1 array'
        go to 3009
      end if
      allocate(Coord2(NLines(2)))      !  allocate mem for Coord2
      if (.not.allocated(Coord2)) then
        print *,'ERROR: Unable to allocate memory for Coord2 array'
        go to 3009
      end if
      allocate(Mstate(NLines(2)))      !  allocate mem for Mstate
      if (.not.allocated(Mstate)) then
        print *,'ERROR: Unable to allocate memory for Mstate array'
        go to 3009
      end if
c
c-----------------------------------------------------------------------
c
      if (BmgColl8) then
        DeltRA = (RAmax(1) - RAmin(1))
     +          * cos(d2r*(DecMax(1)+DecMin(1))/2.0)
        DeltDec = Decmax(1) - Decmin(1)
        DecOrder = DeltDec .gt. DeltRA
        DecOrder = DecOrder
     +       .or. (DecMin(1) .lt. -80.0) .or. (DecMax(1) .gt. 80.0)
        do 96 N = 1, NF1
          if (Field1(N) .eq. SrcID) then
            IS1(1) = IF1a(N)
            IS2(1) = IF1b(N)
          end if
          if ((Field1(N) .eq. 'CnfFlg') .or. (Field1(N) .eq. 'CF')) then
            IC1(1) = IF1a(N)
            IC2(1) = IF1b(N)
          end if
96      continue
        if (IS1(1) .le. 0) then
          print *,'ERROR: can''t find '//SrcID//' in header line:'
          print *,ColNam(1)(1:LNBlnk(ColNam(1)))
          go to 3009
        end if
        if (IC1(1) .le. 0) then
          print *,'ERROR: can''t find CnfFlg or CF in header line:'
          print *,ColNam(1)(1:LNBlnk(ColNam(1)))
          go to 3009
        end if
        NB1 = (1 + NInt(sqrt(1.0 + 4.0*float(IC2(1)-IC1(1)))))/2 ! No. bands
c
        M = 0
        do 97 N = 1, NF2
          if (Field2(N) .eq. SrcID) then
            IS1(2) = IF2a(N)
            IS2(2) = IF2b(N)
            IPa2   = IF2a(N) - 1 ! end of file-2 output piece #1
          end if
          if ((Field2(N) .eq. 'CnfFlg') .or. (Field2(N) .eq. 'CF')) then
            IC1(2) = IF2a(N)
            IC2(2) = IF2b(N)
            IPa1   = IF2a(N)     ! beginning of file-2 output piece #1
            NamDes = Field2(N)
            call MakeNam(NamDes,NamGot,ColNam(2),IF2a(N),IF2b(N),
     +                   Field1,NF1,0)
            Field2(N) = NamGot
          end if
          if (Index(Field2(N),'srcid') .gt. 0) then
            M = M + 1
            if (M .gt. MaxBand) then
              print *,'ERROR: too many srcid values in file 2;'
              print *,'       max =',MaxBand,'; found',M
              go to 3009
            end if
            ISr1(M) = IF2a(N)
            ISr2(M) = IF2b(N)
          end if
97      continue
        if (IS1(2) .le. 0) then
          print *,'ERROR: can''t find '//SrcID//' in header line:'
          print *,ColNam(2)(1:LNBlnk(ColNam(2)))
          go to 3009
        end if
        if (IC1(2) .le. 0) then
          print *,'ERROR: can''t find CnfFlg or CF in header line:'
          print *,ColNam(2)(1:LNBlnk(ColNam(2)))
          go to 3009
        end if
        NB2 = (1 + NInt(sqrt(1.0 + 4.0*float(IC2(2)-IC1(2)))))/2 ! No. bands
        call GetIDB(ColNam(2),NB2,IDB2,OK)        ! Get band numbers
        if (.not.OK) go to 3009
        if (IDBColl8 .eq. IDB2(1)) then           ! collating on first srcid
          IPa2 = ISr1(1) - 1
          IPb1 = ISr1(2)
          IPb2 = IF2b(NF2)
          IF2  = IPb1 - 1
        else if (IDBColl8 .eq. IDB2(NB2)) then    ! collating on last srcid
          IPa2 = ISr1(1) - 1
          IPb1 = ISr1(1)
          IPb2 = ISr2(NB2-1)                      ! assumes nothing follows
          IF2  = IF2b(NF2)                        !  last srcid fields
        else                                      ! collating on mid srcid
          IPb1 = IS2(2) + 1
          IPb2 = IF2b(NF2)
          IF2  = IPb1 - 1
        end if
        IF1 = IS1(1) + IF2 - IS1(2)
        go to 100
      end if ! if (BmgColl8)
c      
c                                      ! must be mdex files in same format,
      if (CatWISE) then                ! hence NF1=NF2, mdetID in same column
        do 98 N = 1, NF1
          if (nm2) then
            if (Field1(N) .eq. 'src') then
              IM1(1) = IF1a(N)
              IM2(1) = IF1b(N)
            end if
          else
            if (Field1(N) .eq. 'mdetID') then
              IM1(1) = IF1a(N)
              IM2(1) = IF1b(N)
            end if          
          end if
98      continue
        if (IM1(1) .le. 0) then
          if (nm2) then
            print *,'ERROR: can''t find "src" in file 1 header line:'
          else
            print *,'ERROR: can''t find "mdetID" in header line:'
          end if
          print *,ColNam(1)(1:LNBlnk(ColNam(1)))
          go to 3009
        end if
        do 99 N = 1, NF2
          if (nm2) then
            if (Field2(N) .eq. 'src') then
              IM1(2) = IF2a(N)
              IM2(2) = IF2b(N)
            end if
          else
            if (Field2(N) .eq. 'mdetID') then
              IM1(2) = IF2a(N)
              IM2(2) = IF2b(N)
            end if
          end if
99      continue
        if (IM1(2) .le. 0) then
          if (nm2) then
            print *,'ERROR: can''t find "src" in file 2 header line:'
          else
            print *,'ERROR: can''t find "mdetID" in header line:'
          end if
          print *,ColNam(2)(1:LNBlnk(ColNam(2)))
          go to 3009
        end if
      end if
c
      If (NLines(1) .gt. NLines(2)) then
        AvgDec = 0.5*(Decmax(1) + Decmin(1))
        if (RAmin(1) .lt. 0.0) RAmin(1) = RAmin(1) + 360.0
        if (RAmax(1) .lt. 0.0) RAmax(1) = RAmax(1) + 360.0
        DeltRA = abs(RAmax(1) - RAmin(1))
        DeltDec = Decmax(1) - Decmin(1)
c       print *,'RA min/max: ', RAmin(1), RAmax(1)
c       print *,'Dec min/max:',Decmin(1), Decmax(1)
      else
        AvgDec = 0.5*(Decmax(2) + Decmin(2))
        if (RAmin(2) .lt. 0.0) RAmin(2) = RAmin(2) + 360.0
        if (RAmax(2) .lt. 0.0) RAmax(2) = RAmax(2) + 360.0
        DeltRA = abs(RAmax(2) - RAmin(2))
        DeltDec = Decmax(2) - Decmin(2)
c       print *,'RA min/max: ', RAmin(1), RAmax(1)
c       print *,'Dec min/max:',Decmin(1), Decmax(1)
      end if
c     print *,'RA delta (true angle):', cos(d2r*AvgDec)*DeltRA
c     print *,'Dec delta:            ', DeltDec
      DecOrder = DeltDec .gt. 0.5*cos(d2r*AvgDec)*DeltRA
      DecOrder = DecOrder .or. (DeltRA .gt. 180.0) .or. CatWISE
     +      .or. ForceDec
      if (cos(d2r*AvgDec) .ne. 0.0) dRA = dRA/cos(d2r*AvgDec)
      if (DecOrder) then
        print *,'Output will be in ascending Dec order'
      else
        print *,'Output will be in ascending RA order'
      end if
c
100   Open(unit = 11, file = TmpNam(1), access = 'DIRECT',
     +     recl = LRecl(1), form = 'UNFORMATTED',
     +     status = 'OLD', err = 3000)
      Open(unit = 13, file = TmpNam(2), access = 'DIRECT',
     +     recl = LRecl(2), form = 'UNFORMATTED',
     +     status = 'OLD', err = 3001)
c
      do 105 N = 1, NLines(1)
        read(11, rec=N, err = 3002) (Chr(L), L = 1, LRecl(1))
        if (DecOrder) then
          read (Line(ID1(1):ID2(1)), *, err = 101) Coord1(N)
          go to 105
101       if (.not.SkipErr) go to 3006
          Coord1(N) = 9.9e9            ! set aside, reject later
        else
          read (Line(IA1(1):IA2(1)), *, err = 102) Coord1(N)
          go to 105
102       if (.not.SkipErr) go to 3003
          Coord1(N) = 9.9e9            ! set aside, reject later
        end if
105   continue
      Call TJISort(NLines(1),Coord1,Indx1) ! Gen. index array for
                                           !  accessing in asc. order
c
      do 110 N = 1, NLines(2)
        read(13, rec=N, err = 3004) (Chr(L), L = 1, LRecl(2))
        if (DecOrder) then
          read (Line(ID1(2):ID2(2)), *, err = 106) Coord2(N)
          go to 110
106       if (.not.SkipErr) go to 3007
          Coord2(N) = 9.9e9            ! set aside, reject later
        else
          read (Line(IA1(2):IA2(2)), *, err = 107) Coord2(N)
          go to 110
107       if (.not.SkipErr) go to 3005
          Coord2(N) = 9.9e9            ! set aside, reject later
        end if
110   continue
      Call TJISort(NLines(2),Coord2,Indx2) ! Gen. index array for
                                           !  accessing in asc. order
c
c-----------------------------------------------------------------------
c
      if (IndxOut)
     + Open (17, File = OutFNam(1:LNBlnk(OutFNam))//'_index')
      if (.not.MergOut) go to 240
c
      if (BmgColl8 .or. BmgHed) then
	    call GetTmpNam(BmgNam,TempDir,OutFNam,OK)
        if (.not.OK) then
	      print *,'Unable to create temporary work file Bandmerge'
          call DelTmp(TmpNam(1))
          call DelTmp(TmpNam(2))
          call exit(64)
        end if
        Open (14, File = BmgNam)
      else
        Open (14, File = OutFNam)
      end if
      if (.not.(TblFil(1) .or. TblFil(2))) then
        Len1   = LRecl(1)
        LenOut = LRecl(1)+LRecl(2)+9
        go to 200                          ! Neither file is table fmt
      end if
      write(14,'(A58)') '\ Generated by gsa vsn '//vsn//' on '
     +                   //CDate//' at '//CTime
      write (14, 6001)
      call MakeFmt(FmtLine,LNBlnk(InFNam1))
      write (14, FmtLine) InFNam1(1:LNBlnk(InFNam1))
      write (14, 6002)
      call MakeFmt(FmtLine,LNBlnk(InFNam2))
      write (14, FmtLine) InFNam2(1:LNBlnk(InFNam2))
      if (BmgColl8) go to 195
      write (14, 6003) Dist
c
c-----------------------------------------------------------------------
c                                            ! Both files are table fmt
      if (TblFil(1) .and. TblFil(2)) then    ! check file-2 column names
        do 120 N = 1, NF2                    ! first against file-1 names
          NamDes = Field2(N)
          call MakeNam(NamDes,NamGot,ColNam(2),IF2a(N),IF2b(N),
     +                 Field1,NF1,0)
          Field2(N) = NamGot
120     continue
        if (NF2 .lt. 2) go to 190
        do 140 N = 2, NF2                    ! then against lower-position
          NamDes = Field2(N)                 ! file-2 names
          call MakeNam(NamDes,NamGot,ColNam(2),IF2a(N),IF2b(N),
     +                 Field2,N-1,0)
          Field2(N) = NamGot
140     continue
c
        if (nColTyp(1) .eq. nColTyp(2)) then
          Nct = nColTyp(1)
        else if (nColTyp(1) .gt. nColTyp(2)) then
          Nct = nColTyp(1)
          if (NColTyp(2) .lt. 2) then
            ColTyp(2) = ColNam(2)
            call ClrCT(ColTyp(2))
          end if
          if (NColTyp(2) .lt. 3) then
            ColTyp3(2) = ColTyp(2)
            call ClrCT(ColTyp3(2))
          end if
          if (NColTyp(2) .lt. 4) then
            ColTyp4(2) = ColTyp(2)
            call ClrCT(ColTyp4(2))
          end if
        else if (nColTyp(1) .lt. nColTyp(2)) then
          Nct = nColTyp(2)
          if (NColTyp(1) .lt. 2) then
            ColTyp(1) = ColNam(1)
            call ClrCT(ColTyp(1))
          end if
          if (NColTyp(1) .lt. 3) then
            ColTyp3(1) = ColTyp(1)
            call ClrCT(ColTyp3(1))
          end if
          if (NColTyp(1) .lt. 4) then
            ColTyp4(1) = ColTyp(1)
            call ClrCT(ColTyp4(1))
          end if
        end if
c
        go to 190
      end if
c
c-----------------------------------------------------------------------
c
c                                 ! only one file is table fmt; set up
      do 150 N = 1, 5000          !  headers for the one that isn't
        Chr(N) = ' '              ! blank out Line
150   continue
c
      if (.not.TblFil(1))  then
        I = 1
        J = 2
      else
        I = 2
        J = 1
      end if
c
      ColNam(I) = Line                      ! set up table headers
      ColTyp(I) = Line
      ColNam(I)(1:1) = '|'
      ColTyp(I)(1:1) = '|'
      ColNam(I)(LRecl(I)+1:LRecl(I)+1) = '|'
      ColTyp(I)(LRecl(I)+1:LRecl(I)+1) = '|'
      ColNam(I)(IA1(I):IA1(I)) = '|'
      ColTyp(I)(IA1(I):IA1(I)) = '|'
      ColNam(I)(IA2(I)+1:IA2(I)+1) = '|'
      ColTyp(I)(IA2(I)+1:IA2(I)+1) = '|'
      ColNam(I)(ID1(I):ID1(I)) = '|'
      ColTyp(I)(ID1(I):ID1(I)) = '|'
      ColNam(I)(ID2(I)+1:ID2(I)+1) = '|'
      ColTyp(I)(ID2(I)+1:ID2(I)+1) = '|'
      if (nColTyp(J) .gt. 2) ColTyp3(I) = ColTyp(I)
      if (nColTyp(J) .gt. 3) ColTyp4(I) = ColTyp(I)
c
      if (IA2(I)-IA1(I) .lt. 4) then
        K = (IA1(I)+IA2(I)+1)/2
        ColTyp(I)(K:K) = 'r'
      else
        K = (IA2(I)-IA1(I)-4)/2 + 1
        ColTyp(I)(IA1(I)+K:IA1(I)+K+3) = 'real'
      end if
c
      if (ID2(I)-ID1(I) .lt. 4) then
        K = (ID1(I)+ID2(I)+1)/2
        ColTyp(I)(K:K) = 'r'
      else
        K = (ID2(I)-ID1(I)-4)/2 + 1
        ColTyp(I)(ID1(I)+K:ID1(I)+K+3) = 'real'
      end if
c
      NamDes = 'RA'
      if (I .eq. 1) then
        call MakeNam(NamDes,NamGot,ColNam(1),IA1(1),IA2(1),Field2,NF2,1)
        NF1 = 1
        Field1(1) = NamGot
      else
        call MakeNam(NamDes,NamGot,ColNam(2),IA1(2),IA2(2),Field1,NF1,1)
        NF2 = 1
        Field2(1) = NamGot
      end if
      NamDes = 'Dec'
      if (I .eq. 1) then
        call MakeNam(NamDes,NamGot,ColNam(1),ID1(1),ID2(1),Field2,NF2,1)
        NF1 = 2
        Field1(2) = NamGot
      else
        call MakeNam(NamDes,NamGot,ColNam(2),ID1(2),ID2(2),Field1,NF1,1)
        NF2 = 2
        Field2(2) = NamGot
      end if
c
      if (IA1(I) .lt. ID1(I)) then        ! if RA precedes Dec, then
        if (IA1(I) .gt. 1) then           ! if RA not 1st column, then
          NamDes = 'data1'
          if (IA1(I) .lt. 7) then         ! if not enough room for
            do 155 n = IA1(I)-1, 5        !  'data1', clip on right
              NamDes(n:n) = ' '
155         continue
          end if
          if (I .eq. 1) then
            call MakeNam(NamDes,NamGot,ColNam(1),2,IA1(1)-1,
     +                   Field2,NF2,1)
            NF1 = NF1 + 1
            Field1(NF1) = NamGot
          else
            call MakeNam(NamDes,NamGot,ColNam(2),2,IA1(2)-1,
     +                   Field1,NF1,1)
            NF2 = NF2 + 1
            Field2(NF2) = NamGot
          end if
          if (IA1(I) .lt. 6) then
            K = (IA1(I)+1)/2
            ColTyp(I)(K:K) = 'c'
          else
            K = (IA1(I)-4)/2 + 1
            ColTyp(I)(K:K+3) = 'char'
          end if
        end if
        if (ID1(I) .gt. IA2(I)+2) then    ! if there's space between
          NamDes = 'data2'                ! 'RA' & 'Dec', then
          if (ID1(I)-IA2(I) .lt. 7) then  ! if not enough room for
            do 160 n = ID1(I)-IA2(I)-1, 5 !  'data2', clip on right
              NamDes(n:n) = ' '
160         continue
          end if
          if (I .eq. 1) then
            call MakeNam(NamDes,NamGot,ColNam(1),IA2(1)+2,ID1(1)-1,
     +                   Field2,NF2,1)
            NF1 = NF1 + 1
            Field1(NF1) = NamGot
          else
            call MakeNam(NamDes,NamGot,ColNam(2),IA2(2)+2,ID1(2)-1,
     +                   Field1,NF1,1)
            NF2 = NF2 + 1
            Field2(NF2) = NamGot
          end if
          if (ID1(I)-IA2(I) .lt. 5) then
            K = (ID1(I)+IA2(I)+1)/2
            ColTyp(I)(K:K) = 'c'
          else
            K = (ID1(I)-IA2(I)-4)/2 + 1
            ColTyp(I)(IA2(I)+K:IA2(I)+K+3) = 'char'
          end if
          NamDes = 'data3'
        else
          NamDes = 'data2'
        end if
        if (ID2(I) .lt. LRecl(I)) then      ! if there's data beyond Dec
          if (LRecl(I)-ID2(I) .lt. 7) then  ! if not enough room for
            do 165 n = LRecl(I)-ID2(I)-1, 5 !  'data3', clip on right
              NamDes(n:n) = ' '
165         continue
          end if
          if (I .eq. 1) then
            call MakeNam(NamDes,NamGot,ColNam(1),ID2(1)+2,LRecl(1),
     +                   Field2,NF2,1)
            NF1 = NF1 + 1
            Field1(NF1) = NamGot
          else
            call MakeNam(NamDes,NamGot,ColNam(2),ID2(2)+2,LRecl(2),
     +                   Field1,NF1,1)
            NF2 = NF2 + 1
            Field2(NF2) = NamGot
          end if
          if (LRecl(I)-ID2(I) .lt. 5) then
            K = (ID2(I)+LRecl(I)+1)/2
            ColTyp(I)(K:K) = 'c'
          else
            K = (LRecl(I)-ID2(I)-4)/2 + 1
            ColTyp(I)(ID2(I)+K:ID2(I)+K+3) = 'char'
          end if
        end if
c
      else                                ! Dec precedes RA, so
        if (ID1(I) .gt. 1) then           ! if Dec not 1st column, then
          NamDes = 'data1'
          if (ID1(I) .lt. 7) then           ! if not enough room for
            do 170 n = IA1(I)-1, 5          !  'data1', clip on right
              NamDes(n:n) = ' '
170         continue
          end if
          if (I .eq. 1) then
            call MakeNam(NamDes,NamGot,ColNam(1),2,ID1(1)-1,
     +                   Field2,NF2,1)
            NF1 = NF1 + 1
            Field1(NF1) = NamGot
          else
            call MakeNam(NamDes,NamGot,ColNam(2),2,ID1(2)-1,
     +                   Field1,NF1,1)
            NF2 = NF2 + 1
            Field2(NF2) = NamGot
          end if
          if (ID1(I) .lt. 6) then
            K = (ID1(I)+1)/2
            ColTyp(I)(K:K) = 'c'
          else
            K = (ID1(I)-4)/2 + 1
            ColTyp(I)(K:K+3) = 'char'
          end if
        end if
        if (IA1(I) .gt. ID2(I)+2) then    ! if there's space between
          NamDes = 'data2'                ! 'Dec' & 'RA', then
          if (IA1(I)-ID2(I) .lt. 7) then  ! if not enough room for
            do 175 n = IA1(I)-ID2(I)-1, 5 !  'data2', clip on right
              NamDes(n:n) = ' '
175         continue
          end if
          if (I .eq. 1) then
            call MakeNam(NamDes,NamGot,ColNam(1),ID2(1)+2,IA1(1)-1,
     +                   Field2,NF2,1)
            NF1 = NF1 + 1
            Field1(NF1) = NamGot
          else
            call MakeNam(NamDes,NamGot,ColNam(2),ID2(2)+2,IA1(2)-1,
     +                   Field1,NF1,1)
            NF2 = NF2 + 1
            Field2(NF2) = NamGot
          end if
          if (IA1(I)-ID2(I) .lt. 5) then
            K = (IA1(I)+ID2(I)+1)/2
            ColTyp(I)(K:K) = 'c'
          else
            K = (IA1(I)-ID2(I)-4)/2 + 1
            ColTyp(I)(ID2(I)+K:ID2(I)+K+3) = 'char'
          end if
          NamDes = 'data3'
        else
          NamDes = 'data2'
        end if
        if (IA2(I) .lt. LRecl(I)) then      ! if there's data beyond RA
          if (LRecl(I)-IA2(I) .lt. 7) then  ! if not enough room for
            do 180 n = LRecl(I)-IA2(I)-1, 5 !  'data3', clip on right
              NamDes(n:n) = ' '
180         continue
          end if
          if (I .eq. 1) then
            call MakeNam(NamDes,NamGot,ColNam(1),IA2(1)+2,LRecl(1),
     +                   Field2,NF2,1)
            NF1 = NF1 + 1
            Field1(NF1) = NamGot
          else
            call MakeNam(NamDes,NamGot,ColNam(2),IA2(2)+2,LRecl(2),
     +                   Field1,NF1,1)
            NF2 = NF2 + 1
            Field2(NF2) = NamGot
          end if
          if (LRecl(I)-IA2(I) .lt. 5) then
            K = (IA2(I)+LRecl(I)+1)/2
            ColTyp(I)(K:K) = 'c'
          else
            K = (LRecl(I)-IA2(I)-4)/2 + 1
            ColTyp(I)(IA2(I)+K:IA2(I)+K+3) = 'char'
          end if
        end if
      end if     !  RA precedes Dec
c
c-----------------------------------------------------------------------
c                                            ! write full (concatenated)
190   NamDes = 'dist'                        ! header lines, including
      Line1 = ' '                            ! 'dist' parameter
      call MakeNam(NamDes,NamGot,Line1,1,8,Field1,NF1,1)
      NamDes = NamGot
      Line1 = ' '
      call MakeNam(NamDes,NamGot,Line1,1,8,Field2,NF2,1)
c
195   if (BmgColl8) then
        call MakeFmt(FmtLine,LNBlnk(ColNam(1))+IPA2-IPa1+IPb2-IPb1+2)
        write (14, FmtLine) ColNam(1)(1:LNBlnk(ColNam(1))-1)
     +                    //ColNam(2)(IPa1:IPa2)//ColNam(2)(IPb1:IPb2+1)
        write (14, FmtLine) ColTyp(1)(1:LNBlnk(ColTyp(1))-1)
     +                    //ColTyp(2)(IPa1:IPa2)//ColTyp(2)(IPb1:IPb2+1)
        Len1   = LNBlnk(ColNam(1)) - 1
        LenOut = Len1 + IPA2 - IPa1 + IPb2 - IPb1 + 2
        go to 200
      end if
c
      call MakeFmt(FmtLine,LNBlnk(ColNam(1))+LNBlnk(ColNam(2))+8)
      write (14, FmtLine) ColNam(1)(1:LNBlnk(ColNam(1)))//Line1(1:8)
     +                  //ColNam(2)(1:LNBlnk(ColNam(2)))
      write (14, FmtLine) ColTyp(1)(1:LNBlnk(ColTyp(1)))//'  real  '
     +                  //ColTyp(2)(1:LNBlnk(ColTyp(2)))
      if (Nct .gt. 2)
     + write (14, FmtLine) ColTyp3(1)(1:LNBlnk(ColTyp3(1)))//'  asec  '
     +                   //ColTyp3(2)(1:LNBlnk(ColTyp3(2)))
      if (Nct .gt. 3)
     + write (14, FmtLine) ColTyp4(1)(1:LNBlnk(ColTyp4(1)))//'  null  '
     +                   //ColTyp4(2)(1:LNBlnk(ColTyp4(2)))
      Len1   = LNBlnk(ColNam(1)) - 1
      LenOut = Len1 + LRecl(2) + 9
c
c-----------------------------------------------------------------------
c                                            ! set up null records
200   if (.not.AllSrc2) go to 220
      N1 = NLines(1)/2 - 25
      N2 = NLines(1)/2 + 25
      if (N1 .lt. 1) N1 = 1
      if (N2 .gt. NLines(1)) N2 = NLines(1)
      do 205 N = 1, 5000
        Chr(N) = ' '              ! blanks out Line
205   continue
      Null1 = Line
c
      do 210 N = N1, N2
        read(11, rec=N, err = 3002) (Chr(L), L = 1, LRecl(1))
        do 208 L = 1, LRecl(1)
          if (Chr(L) .ne. ' ') Null1(L:L) = '9'
208     continue
210   continue
c
      do 215 L = 1, LRecl(1)-2
        if ((Null1(L:L) .eq. ' ') .and. (Null1(L+1:L+1) .eq. '9')) then
          if (Null1(L+2:L+2) .ne. ' ') then
            Null1(L+1:L+1) = '-'
          else
            Null1(L+1:L+1) = '0'
          end if
        end if
215   continue
c
      if (BmgColl8) then
        do 217 N = 1, NF1
          if (Index(Field1(N),'srcid') .gt. 0) then
            do 216 L = IF1a(N), IF1b(N)-1
              Null1(L:L) = ' '
216         continue
            Null1(IF1b(N):IF1b(N)) = '0'
          end if
217     continue
        do 218 L = IC1(1)+1, IC2(1)
          Null1(L:L) = '0'
218     continue
      end if
c
220   if (.not.AllSrc1) go to 240
      N1 = NLines(2)/2 - 25
      N2 = NLines(2)/2 + 25
      if (N1 .lt. 1) N1 = 1
      if (N2 .gt. NLines(2)) N2 = NLines(2)
      do 225 N = 1, 5000
        Chr(N) = ' '              ! blanks out Line
225   continue
      Null2 = Line
c
      do 230 N = N1, N2
        read(13, rec=N, err = 3004) (Chr(L), L = 1, LRecl(2))
        do 228 L = 1, LRecl(2)
          if (Chr(L) .ne. ' ') Null2(L:L) = '9'
228     continue
230   continue
c
      do 235 L = 1, LRecl(2)-2
        if ((Null2(L:L) .eq. ' ') .and. (Null2(L+1:L+1) .eq. '9')) then
          if (Null2(L+2:L+2) .ne. ' ') then
            Null2(L+1:L+1) = '-'
          else
            Null2(L+1:L+1) = '0'
          end if
        end if
235   continue
c
      if (BmgColl8) then
        do 237 N = 1, NF2
          if (Index(Field2(N),'srcid') .gt. 0) then
            do 236 L = IF2a(N), IF2b(N)-1
              Null2(L:L) = ' '
236         continue
            Null2(IF2b(N):IF2b(N)) = '0'
          end if
237     continue
        do 238 L = IC1(2)+1, IC2(2)
          Null2(L:L) = '0'
238     continue
      end if
c
c-----------------------------------------------------------------------
c                                             ! get file vectors
240   do 250 N = 1, NLines(1)
        read(11, rec=N, err = 3002) (Chr(L), L = 1, LRecl(1))
        if (BmgColl8) then
          read (Line(IS1(1):IS2(1)), *, err = 3011) X1(N)
          Y1(N) = 0.0
          Z1(N) = 0.0
          go to 250
        end if
        read (Line(IA1(1):IA2(1)), *, err = 241) Alpha
        go to 242
241     if (.not.SkipErr) go to 3003
        go to 244       
242     read (Line(ID1(1):ID2(1)), *, err = 243) Delta
        go to 246
243     if (.not.SkipErr) go to 3006
244     X1(N) = 0.0d0
        Y1(N) = 0.0d0
        Z1(N) = 0.0d0
        go to 250 
246     cosD  =  dcos(d2r*Delta)
        X1(N) =  dsin(d2r*Delta)
        Y1(N) = -cosD*dsin(d2r*Alpha)
        Z1(N) =  cosD*dcos(d2r*Alpha)
250   continue
c
      do 260 N = 1, NLines(2)
        read(13, rec=N, err = 3004) (Chr(L), L = 1, LRecl(2))
        if (BmgColl8) then
          read (Line(IS1(2):IS2(2)), *, err = 3012) X2(N)
          if (X2(N) .le. 0.0) X2(N) = -9.9
          Y2(N) = 0.0
          Z2(N) = 0.0
          go to 255
        end if
        read (Line(IA1(2):IA2(2)), *, err = 251) Alpha
        go to 252
251     if (.not.SkipErr) go to 3005
        go to 254       
252     read (Line(ID1(2):ID2(2)), *, err = 253) Delta
        go to 256
253     if (.not.SkipErr) go to 3007
254     X2(N) = 1.0d0
        Y2(N) = 1.0d0
        Z2(N) = 1.0d0
        go to 255 
256     Alpha = Alpha - dRA
        Delta = Delta - dDec
        cosD =  dcos(d2r*Delta)
        X2(N) =  dsin(d2r*Delta)      - dX
        Y2(N) = -cosD*dsin(d2r*Alpha) - dY
        Z2(N) =  cosD*dcos(d2r*Alpha) - dZ
255     Mstate(N) = 0           ! Init state to "unmatchable"
260   continue
c
c-----------------------------------------------------------------------
c                                              ! set up reject files
265   if (Reject1) then
        open (15, File = Rej1Nam)
        if (TblFil(1)) then
          write(15,'(A58)') '\ Generated by gsa vsn '//vsn//' on '
     +                       //CDate//' at '//CTime
          write (15, 6004)
          call MakeFmt(FmtLine,LNBlnk(InFNam1))
          write (15, FmtLine) InFNam1(1:LNBlnk(InFNam1))
          call MakeFmt(FmtLine,LNBlnk(ColNam(1)))
          write (15, FmtLine) ColNam(1)(1:LNBlnk(ColNam(1)))
          if (nColTyp(1) .gt. 1)
     +    write (15, FmtLine) ColTyp(1)(1:LNBlnk(ColTyp(1)))
          if (nColTyp(1) .gt. 2)
     +    write (15, FmtLine) ColTyp3(1)(1:LNBlnk(ColTyp3(1)))
          if (nColTyp(1) .gt. 3)
     +    write (15, FmtLine) ColTyp4(1)(1:LNBlnk(ColTyp4(1)))
        end if
        call MakeFmt(FmtLin1,Lrecl(1))
      end if
c
      if (Reject2) then
        open (16, File = Rej2Nam)
        if (TblFil(2)) then
          write(16,'(A58)') '\ Generated by gsa vsn '//vsn//' on '
     +                       //CDate//' at '//CTime
          write (16, 6004)
          call MakeFmt(FmtLine,LNBlnk(InFNam2))
          write (16, FmtLine) InFNam2(1:LNBlnk(InFNam2))
          call MakeFmt(FmtLine,LNBlnk(ColNam(3)))
          write (16, FmtLine) ColNam(3)(1:LNBlnk(ColNam(3)))
          if (nColTyp(2) .gt. 1)
     +    write (16, FmtLine) ColTyp(2)(1:LNBlnk(ColTyp(2)))
          if (nColTyp(2) .gt. 2)
     +    write (16, FmtLine) ColTyp3(2)(1:LNBlnk(ColTyp3(2)))
          if (nColTyp(2) .gt. 3)
     +    write (16, FmtLine) ColTyp4(2)(1:LNBlnk(ColTyp4(2)))
        end if
        call MakeFmt(FmtLin2,Lrecl(2))
      end if
c
c-----------------------------------------------------------------------
c
      if (BmgColl8) then
        ADist = 0.0
        Dist  = 0.5
      else
        ADist = Dist/3600.0           ! Dist in deg
        Dist  = d2r*Dist/3600.0d0     ! Convert Dist to radians
      end if
      Window = 1.5*Dist               ! Coarse window slightly larger
c
c                                     ! identify matchable
      do 320 M = 1, NLines(2)         ! file-2 sources
        do 310 N = 1, NLines(1)
          if (abs(X1(N)-X2(M)) .gt. Window) go to 310
          if (BmgColl8) then
            Mstate(M) = 1
            go to 320
          end if
          if (abs(Y1(N)-Y2(M)) .gt. Window) go to 310
          if (abs(Z1(N)-Z2(M)) .gt. Window) go to 310
          Cross = dsqrt((Y1(N)*Z2(M) - Y2(M)*Z1(N))**2
     +                + (Z1(N)*X2(M) - Z2(M)*X1(N))**2
     +                + (X1(N)*Y2(M) - X2(M)*Y1(N))**2)
          Ang = dasin(Cross)
          if (Ang .le. Dist) then
            Mstate(M) = 1       ! Set state to "matchable"
            go to 320
          end if
310     continue
320   continue
c
c-----------------------------------------------------------------------
c                                              ! do the associating
      if (MergOut) call MakeFmt(FmtLine,LenOut)
c
      do 400 N = 1, NLines(1)
        N1 = Indx1(N)
        if (SkipErr .and. (X1(N1)+Y1(N1)+Z1(N1) .eq. 0.0d0)) go to 400
        read(11, rec=N1, err = 3002) (Chr(L), L = 1, LRecl(1))
        Line1 = Line(1:LRecl(1))//'                                    '
        read (Line1(IA1(1):IA2(1)), *, err = 3003) Alpha1
        read (Line1(ID1(1):ID2(1)), *, err = 3006) Delta1
        if (CatWISE) read (Line1(IM1(1):IM2(1)), *, err = 3013) mdetID1
        AngMin = 9.9e9
        N2best = 0
        GotMatch = .False.
c
        do 360 M = 1, NLines(2)
          N2 = Indx2(M)
          if (Mstate(N2) .eq. 0) go to 360 ! skip unmatchable recs
          if (Mstate(N2) .eq. 2) go to 360 ! skip output unmatchable recs
          if ((Mstate(N2) .eq. 3) .and.    ! skip already-matched recs
     +        .not.AllAssn)                !  if multi-matches not
     +        go to 360                    !  allowed
          if (BmgColl8) then
            if (X1(N1) .ne. X2(N2)) go to 360
            Ang = 0.0
            go to 330
          end if
          if (abs(X1(N1)-X2(N2)) .gt. Window) go to 360
          if (abs(Y1(N1)-Y2(N2)) .gt. Window) go to 360
          if (abs(Z1(N1)-Z2(N2)) .gt. Window) go to 360
          Cross = sqrt((Y1(N1)*Z2(N2) - Y2(N2)*Z1(N1))**2
     +               + (Z1(N1)*X2(N2) - Z2(N2)*X1(N1))**2
     +               + (X1(N1)*Y2(N2) - X2(N2)*Y1(N1))**2)
          Ang = dasin(Cross)
          if (Ang .gt. Dist) go to 360
          if (CatWISE) then
            read (13, rec=N2, err = 3004) (Chr(L), L = 1, LRecl(2))
            read (Line(IM1(2):IM2(2)), *, err = 3014) mdetID2
            if (mdetID1 .ne. mdetID2) go to 360
          end if
          if (.not.AllAssn) then              ! This match is acceptable on
            if (Best2Come(X1,Y1,Z1,           !  the basis of position; if
     +          X2,Y2,Z2,Indx1,N2,Window,     !  single best-match only, check
     +          NLines,N,ADist,Ang,           !  whether a downstream primary
     +          mdetID2,CatWISE,LRecL(1),     !  matches the secondary better
     +          IM1(1),IM2(1),
     +          Coord1,Coord2))
     +      go to 360
          end if
330       if (GotMatch) N1multi = N1multi + 1
          GotMatch = .True.
          if (MState(N2) .eq. 3) N2multi = N2multi + 1
          if (AllAssn) then                ! output all associations as we go
            if ((AllSrc2 .or. Reject2)     ! check for unmatchable secondaries
     +           .and. (M .gt. 1)) then    !   ahead of the matched one
              do 340 K = 1, M-1
                N3 = Indx2(K)
                if (Mstate(N3) .eq. 0) then
                  N2single = N2single + 1
                  read(13, rec=N3, err = 3004) (Chr(L), L = 1, LRecl(2))
                  if (AllSrc2) then
                    NTotal = NTotal + 1
                    if (MergOut) then
                      if (BmgColl8) then
                        TmpStr = Line(1:IH2)
                        call FillIn(Null1,TmpStr,1,IH2)
                        TmpStr = Line(IS1(2):IF2)
                        call FillIn(Null1,TmpStr,IS1(1),IF1)
                        do 335 L = 2, 7
                          if (Null1(L:L) .eq. ' ') Null1(L:L) = '9'
335                     continue
                        write (14, FmtLine)
     +                  Null1(1:Len1)//Line(IPa1:IPa2)//Line(IPb1:IPb2)
                      else
                        TmpStr = Line(IA1(2):IA2(2))
                        call FillIn(Null1,TmpStr,IA1(1),IA2(1))
                        TmpStr = Line(ID1(2):ID2(2))
                        call FillIn(Null1,TmpStr,ID1(1),ID2(1))
                        write (14, FmtLine)
     +                  Null1(1:Len1)//'     -1.0'//Line(1:LRecl(2))
                      end if
                    end if
                    if (IndxOut) write (17, 6005) Izero, N3
                  end if
                  if (Reject2) write (16, FmtLin2) Line(1:LRecl(2))
                  Mstate(N3) = 2           ! tag: unmatchable & output
                end if
340           continue
            end if
            NAssns = NAssns + 1
            NTotal = NTotal + 1
            read(13, rec=N2, err = 3004) (Chr(L), L = 1, LRecl(2))
            if (BmgColl8) go to 350
            write (DistStr,'(F9.3)') 3600.0d0*Ang/d2r
            read (Line(IA1(2):IA2(2)), *, err = 3005) Alpha
            read (Line(ID1(2):ID2(2)), *, err = 3007) Delta
            Alpha = Alpha - dRA
            Delta = Delta - dDec
            DelTARA = Alpha - Alpha1
            if (DelTARA .lt. -180.0d0) DelTARA = DelTARA + 360.0d0
            if (DelTARA .gt.  180.0d0) DelTARA = DelTARA - 360.0d0
            DelTARA = dcos(0.5*d2r*(Delta1 + Delta))*DelTARA
            if (abs(DelTARA) .gt. abs(DelTARAmax)) DelTARAmax = DelTARA
            SumTARA = SumTARA + DelTARA
            SumDist = SumDist + Ang
            if (Ang .gt. DistMax) DistMax = Ang
            SumDec = SumDec  + (Delta - Delta1)
            if (abs(Delta-Delta1) .gt. abs(DeltDecMax))
     +        DeltDecMax = Delta - Delta1
            SumSqTARA = SumSqTARA + DelTARA**2
            SumSqDist = SumSqDist + Ang**2
            SumSqDec  = SumSqDec  + (Delta - Delta1)**2
            SumDX = SumDX + X2(N2) - X1(N1)
            SumDY = SumDY + Y2(N2) - Y1(N1)
            SumDZ = SumDZ + Z2(N2) - Z1(N1)
            SumSqDX = SumSqDX + (X2(N2) - X1(N1))**2
            SumSqDY = SumSqDY + (Y2(N2) - Y1(N1))**2
            SumSqDZ = SumSqDZ + (Z2(N2) - Z1(N1))**2
            if (abs(X2(N2) - X1(N1)) .gt. abs(DXmax))
     +        DXmax = X2(N2) - X1(N1)
            if (abs(Y2(N2) - Y1(N1)) .gt. abs(DYmax))
     +        DYmax = Y2(N2) - Y1(N1)
            if (abs(Z2(N2) - Z1(N1)) .gt. abs(DZmax))
     +        DZmax = Z2(N2) - Z1(N1)
350         Mstate(N2) = 3                 ! Set state to "matched & output"
            if (MergOut) then
              if (BmgColl8) then
                write (14, FmtLine)
     +          Line1(1:Len1)//Line(IPa1:IPa2)//Line(IPb1:IPb2)
              else
                write (14, FmtLine)
     +          Line1(1:Len1)//DistStr//Line(1:LRecl(2))
              end if
            end if
            if (IndxOut) write (17, 6005) N1, N2
          else if (Ang .lt. AngMin) then   ! keep best for later
            AngMin = Ang
            N2best = N2
          end if
360     continue
c
        if (.not.GotMatch .and. (AllSrc2 .or. Reject2)) then
          do 370 K = 1, NLines(2)       ! check for unmatched secondaries
            N3 = Indx2(K)
            if (Coord2(N3) .gt. Coord1(N1)) go to 375
            if ((Mstate(N3) .eq. 0) .or.
     +          ((Mstate(N3) .eq. 1) .and. (.not.AllAssn) .and.
     +           (Coord1(N1)-Coord2(N3) .gt. ADist))) then
              N2single = N2single + 1
              read(13, rec=N3, err = 3004) (Chr(L), L = 1, LRecl(2))
              if (AllSrc2) then
                NTotal = NTotal + 1
                if (MergOut) then
                  if (BmgColl8) then
                    TmpStr = Line(1:IH2)
                    call FillIn(Null1,TmpStr,1,IH2)
                    TmpStr = Line(IS1(2):IF2)
                    call FillIn(Null1,TmpStr,IS1(1),IF1)
                    do 365 L = 2, 7
                      if (Null1(L:L) .eq. ' ') Null1(L:L) = '9'
365                 continue
                    write (14, FmtLine)
     +              Null1(1:Len1)//Line(IPa1:IPa2)//Line(IPb1:IPb2)
                  else
                    TmpStr = Line(IA1(2):IA2(2))
                    call FillIn(Null1,TmpStr,IA1(1),IA2(1))
                    TmpStr = Line(ID1(2):ID2(2))
                    call FillIn(Null1,TmpStr,ID1(1),ID2(1))
                    write (14, FmtLine)
     +              Null1(1:Len1)//'     -1.0'//Line(1:LRecl(2))
                  end if
                end if
                if (IndxOut) write (17, 6005) Izero, N3
              end if
              if (Reject2) write (16, FmtLin2) Line(1:LRecl(2))
              Mstate(N3) = 2           ! tag: unmatchable & output
            end if
370       continue
        end if
c
375     if (.not.GotMatch .and. (AllSrc1 .or. Reject1)) then
          N1single = N1single + 1
          if (AllSrc1) then
            NTotal   = NTotal + 1
            if (MergOut) then
              if (BmgColl8) then
                write (14, FmtLine)
     +          Line1(1:Len1)//Null2(IPa1:IPa2)//Null2(IPb1:IPb2)
              else
                TmpStr = Line1(IA1(1):IA2(1))
                call FillIn(Null2,TmpStr,IA1(2),IA2(2))
                TmpStr = Line1(ID1(1):ID2(1))
                call FillIn(Null2,TmpStr,ID1(2),ID2(2))
                write (14, FmtLine)
     +          Line1(1:Len1)//'     -2.0'//Null2(1:LRecl(2))
              end if
            end if
            if (IndxOut) write (17, 6005) N1, Izero
          end if
          if (Reject1) write (15, FmtLin1) Line1(1:LRecl(1))
          go to 400
        end if
        if (.not.AllAssn .and. (N2best .gt. 0))   ! we have a match and it
     +  then                                      ! has not been output yet
          if (AllSrc2 .or. Reject2) then  ! check for unmatched secondaries
            do 380 K = 1, NLines(2)
              N3 = Indx2(K)
              if (Coord2(N3) .gt. Coord1(N1)) go to 385
              if ((Mstate(N3) .eq. 0) .or.
     +            ((Mstate(N3) .eq. 1) .and. (.not.AllAssn) .and.
     +             (Coord1(N1)-Coord2(N3) .gt. ADist))) then
                N2single = N2single + 1
                read(13, rec=N3, err = 3004) (Chr(L), L = 1, LRecl(2))
                if (AllSrc2) then
                  NTotal   = NTotal + 1
                  if (MergOut) then
                    if (BmgColl8) then
                      TmpStr = Line(1:IH2)
                      call FillIn(Null1,TmpStr,1,IH2)
                      TmpStr = Line(IS1(2):IF2)
                      call FillIn(Null1,TmpStr,IS1(1),IF1)
                      do 376 L = 2, 7
                        if (Null1(L:L) .eq. ' ') Null1(L:L) = '9'
376                   continue
                      write (14, FmtLine)
     +                Null1(1:Len1)//Line(IPa1:IPa2)//Line(IPb1:IPb2)
                    else
                      TmpStr = Line(IA1(2):IA2(2))
                      call FillIn(Null1,TmpStr,IA1(1),IA2(1))
                      TmpStr = Line(ID1(2):ID2(2))
                      call FillIn(Null1,TmpStr,ID1(1),ID2(1))
                      write (14, FmtLine)
     +                Null1(1:Len1)//'     -1.0'//Line(1:LRecl(2))
                    end if
                  end if
                  if (IndxOut) write (17, 6005) Izero, N3
                end if
                if (Reject2) write (16, FmtLin2) Line(1:LRecl(2))
                Mstate(N3) = 2           ! tag: unmatchable & output
              end if
380         continue
          end if
385       write (DistStr,'(F9.3)') 3600.0d0*AngMin/d2r
          NAssns = NAssns + 1
          NTotal = NTotal + 1
          read(13, rec=N2best, err = 3004) (Chr(L), L = 1, LRecl(2))
          read (Line(IA1(2):IA2(2)), *, err = 3005) Alpha
          read (Line(ID1(2):ID2(2)), *, err = 3007) Delta
          Alpha = Alpha - dRA
          Delta = Delta - dDec
          DelTARA = Alpha - Alpha1
          if (DelTARA .lt. -180.0d0) DelTARA = DelTARA + 360.0d0
          if (DelTARA .gt.  180.0d0) DelTARA = DelTARA - 360.0d0
          DelTARA = dcos(0.5*d2r*(Delta1 + Delta))*DelTARA
          if (abs(DelTARA) .gt. abs(DelTARAmax)) DelTARAmax = DelTARA
          SumTARA = SumTARA + DelTARA
          SumDist = SumDist + AngMin
          if (AngMin .gt. DistMax) DistMax = AngMin
          SumDec = SumDec  + (Delta - Delta1)
          if (abs(Delta-Delta1) .gt. abs(DeltDecMax))
     +      DeltDecMax = Delta-Delta1
          SumSqTARA = SumSqTARA + DelTARA**2
          SumSqDist = SumSqDist + AngMin**2
          SumSqDec  = SumSqDec  + (Delta - Delta1)**2
          SumDX = SumDX + X2(N2best) - X1(N1)
          SumDY = SumDY + Y2(N2best) - Y1(N1)
          SumDZ = SumDZ + Z2(N2best) - Z1(N1)
          SumSqDX = SumSqDX + (X2(N2best) - X1(N1))**2
          SumSqDY = SumSqDY + (Y2(N2best) - Y1(N1))**2
          SumSqDZ = SumSqDZ + (Z2(N2best) - Z1(N1))**2
          if (abs(X2(N2best) - X1(N1)) .gt. abs(DXmax))
     +      DXmax = X2(N2best) - X1(N1)
          if (abs(Y2(N2best) - Y1(N1)) .gt. abs(DYmax))
     +      DYmax = Y2(N2best) - Y1(N1)
          if (abs(Z2(N2best) - Z1(N1)) .gt. abs(DZmax))
     +      DZmax = Z2(N2best) - Z1(N1)
          Mstate(N2best) = 3               ! Set state to "matched & output"
          if (MergOut) then
            if (BmgColl8) then
              write (14, FmtLine)
     +        Line1(1:Len1)//Line(IPa1:IPa2)//Line(IPb1:IPb2)
            else
              write (14, FmtLine)
     +        Line1(1:Len1)//DistStr//Line(1:LRecl(2))
            end if
          end if
          if (IndxOut) write (17, 6005) N1, N2best
        end if
400   continue
c
      if (AllSrc2 .or. Reject2) then  ! check for remaining
        do 420 K = 1, NLines(2)       !  unmatched secondaries
          N3 = Indx2(K)
          if (Mstate(N3) .lt. 2) then
            N2single = N2single + 1
            read(13, rec=N3, err = 3004) (Chr(L), L = 1, LRecl(2))
            if (AllSrc2) then
              NTotal = NTotal + 1
              if (MergOut) then
                if (BmgColl8) then
                  TmpStr = Line(1:IH2)
                  call FillIn(Null1,TmpStr,1,IH2)
                  TmpStr = Line(IS1(2):IF2)
                  call FillIn(Null1,TmpStr,IS1(1),IF1)
                  do 410 L = 2, 7
                    if (Null1(L:L) .eq. ' ') Null1(L:L) = '9'
410               continue
                  write (14, FmtLine)
     +            Null1(1:Len1)//Line(IPa1:IPa2)//Line(IPb1:IPb2)
                else
                  TmpStr = Line(IA1(2):IA2(2))
                  call FillIn(Null1,TmpStr,IA1(1),IA2(1))
                  TmpStr = Line(ID1(2):ID2(2))
                  call FillIn(Null1,TmpStr,ID1(1),ID2(1))
                  write (14, FmtLine)
     +            Null1(1:Len1)//'     -1.0'//Line(1:LRecl(2))
                end if
              end if
              if (IndxOut) write (17, 6005) Izero, N3
            end if
            if (Reject2) write (16, FmtLin2) Line(1:LRecl(2))
            Mstate(N3) = 2
          end if
420     continue
      end if
c
c-----------------------------------------------------------------------
c
      close(11)
      close(13)
      call DelTmp(TmpNam(1))
      call DelTmp(TmpNam(2))
      print *
      print *,'No. of input sources:        ',NLines
      print *,'No. of associations:         ',NAssns
      if (AllSrc1 .or. Reject1)
     +  print *,'No. unmatched in file 1:     ',N1single
      if (AllAssn .and. (NAssns .gt. 0))
     + print *,'No. multi-matched in file 1: ',N1multi
      if (AllSrc2 .or. Reject2)
     +  print *,'No. unmatched in file 2:     ',N2single
      if (AllAssn .and. (NAssns .gt. 0))
     + print *,'No. multi-matched in file 2: ',N2multi
      print *,'Total output data lines:     ',NTotal
      if (NAssns .eq. 0) go to 2000
      if (nm2 .and. (N2Multi .gt. 0)) then
1001    n = index(InFNam2,'/')
        if (n .gt. 0) then
          do 1002 i = 1, n
          InFNam2(i:i) = ' '
1002      continue
          go to 1001
        end if
        InFNam2 = AdjustL(InFNam2)
        open (33, file = 'ERROR_MESSAGE_gsa-'
     +                  //InFNam2(1:lnblnk(InFNam2))//'.txt')
        write (33,'(a)')
     +  'ERROR: "-nm2" was specified, but there were multi-matches in'
        write (33,'(a,i9)') 'file #2; no. multi-matches =',N2multi
      end if
c
      if (BmgColl8) go to 2000
      print *
      print *,
     +'Position discrepancy statistics for associated sources (arcsec):'
      print *,'                                      Standard'
      print *,
     +'Parameter                 Mean       Deviation      Max Mag'
      RTmp1 = SumDist/DFloat(NAssns)
      RTmp2 = sqrt(abs(SumSqDist/DFloat(NAssns) - RTmp1**2))
      RTmp1 = 3600.0*RTmp1/d2r
      RTmp2 = 3600.0*RTmp2/d2r
      DistMax = 3600.0*DistMax/d2r
      write(OutStr,'(''Radial Distance: '',3F14.4)')
     +  RTmp1, RTmp2, DistMax
      print *,OutStr
      RTmp3 = SumTARA/DFloat(NAssns)
      RTmp4 = sqrt(abs(SumSqTARA/DFloat(NAssns) - RTmp3**2))
      RTmp3 = 3600.0*RTmp3
      RTmp4 = 3600.0*RTmp4
      DelTARAmax = 3600.0*DelTARAmax
      write(OutStr,'(''RA (True Angle): '',3F14.4)')
     +  RTmp3, RTmp4, DelTARAmax
      print *,OutStr
      RTmp5 = SumDec/DFloat(NAssns)
      RTmp6 = sqrt(abs(SumSqDec/DFloat(NAssns) - RTmp5**2))
      RTmp5 = 3600.0*RTmp5
      RTmp6 = 3600.0*RTmp6
      DeltDecMax = 3600.0*DeltDecMax
      write(OutStr,'(''Declination:     '',3F14.4)')
     +  RTmp5, RTmp6, DeltDecMax
      print *,OutStr
c
      if (GotStat) then
        open (20, file = StatNam)
        write(20,'(A58)') '\ Generated by gsa vsn '//vsn//' on '
     +                     //CDate//' at '//CTime
        OutStr = '| mnradd | stdradd|  mnra  |  stdra |  mndec '
     +         //'| stddec |Nassns|'
        write (20,'(A62)') OutStr(1:62)
        OutStr = '|    r   |    r   |    r   |    r   |    r   '
     +         //'|    r   |   i  |'
        write (20,'(A62)') OutStr(1:62)
        OutStr = '|  asec  |  asec  |  asec  |  asec  |  asec  '
     +         //'|  asec  |      |'
        write (20,'(A62)') OutStr(1:62)
        write (20,'(6F9.4,I7)') RTmp1, RTmp2, RTmp3, RTmp4,
     +    RTmp5, RTmp6, NAssns
      end if
c
      print *
      print *,
     +'Unit-vector discrepancy statistics for associated sources:'
      print *,'                                      Standard'
      print *,
     +'Axis                      Mean       Deviation      Max Mag'
      RTmp1 = SumDX/DFloat(NAssns)
      RTmp2 = sqrt(abs(SumSqDX/DFloat(NAssns) - RTmp1**2))
      write(OutStr,'(''X2-X1:           '',1pE14.4,1pE16.4,1pE13.3)')
     +  RTmp1, RTmp2, DXmax
      print *,OutStr
      RTmp1 = SumDY/DFloat(NAssns)
      RTmp2 = sqrt(abs(SumSqDY/DFloat(NAssns) - RTmp1**2))
      write(OutStr,'(''Y2-Y1:           '',1pE14.4,1pE16.4,1pE13.3)')
     +  RTmp1, RTmp2, DYmax
      print *,OutStr
      RTmp1 = SumDZ/DFloat(NAssns)
      RTmp2 = sqrt(abs(SumSqDZ/DFloat(NAssns) - RTmp1**2))
      write(OutStr,'(''Z2-Z1:           '',1pE14.4,1pE16.4,1pE13.3)')
     +  RTmp1, RTmp2, DZmax
      print *,OutStr
c
2000  print *
      if (BmgColl8 .or. BmgHed) then
	  call GetTmpNam(TmpNam(3),TempDir,OutFNam,OK)
      if (.not.OK) then
	    print *,'Unable to create temporary work file no. 3'
        call exit(64)
      end if
        Open (18, File = TmpNam(3))
        if (BmgColl8) then
          M = NB1 + NB2 - 1
        else
          M = NBands
        end if
        write(18,'(''\INT NBands = '',I2)') M
        write(18,'(''\INT Total_Merged_Number ='',I7)') NTotal
        close(18)
        ISys = System('cat '//TmpNam(3)//' '//BmgNam//' > '//OutFNam)
        ISys = System('rm '//BmgNam)
        ISys = System('rm '//TmpNam(3))
      end if
c
      call signoff('gsa')
      stop
c
c-----------------------------------------------------------------------
c
3000  print *,'ERROR opening direct-access file #1'
      go to 3009
c
3001  print *,'ERROR opening direct-access file #2'
      go to 3009
c
3002  print *,'ERROR read direct-access file #1, record ', N
      go to 3009
c
3003  print *,'ERROR decoding RA in file #1, record ', N
      print *,Line(IA1(1):IA2(1))
      go to 3009
c
3004  print *,'ERROR read direct-access file #2, record ', N
      go to 3009
c
3005  print *,'ERROR decoding RA in file #2, record ', N
      print *,Line(IA1(2):IA2(2))
      go to 3009
c
3006  print *,'ERROR decoding Dec in file #1, record ', N
      print *,Line(ID1(1):ID2(1))
      go to 3009
c
3007  print *,'ERROR decoding Dec in file #2, record ', N
      print *,Line(ID1(2):ID2(2))
      go to 3009
c
3009  call DelTmp(TmpNam(1))
      call DelTmp(TmpNam(2))
      call exit(64)
      stop
c
3010  print *,'ERROR: expected command-line specification after ',Flag
      print *,'       but nothing found'
      NeedHelp = .True.
      go to 1
c
3011  print *,'ERROR decoding srcid in file #1, record ', N
      print *,Line(ID1(1):ID2(1))
      go to 3009
c
3012  print *,'ERROR decoding srcid in file #2, record ', N
      print *,Line(ID1(1):ID2(1))
      go to 3009
c
3013  print *,'ERROR decoding mdetID in file #1, record ', N
      print *,Line(IM1(1):IM2(1))
      go to 3009
c
3014  print *,'ERROR decoding mdetID in file #2, record ', N
      print *,Line(IM1(2):IM2(2))
      go to 3009
c
c-----------------------------------------------------------------------
c
6001  Format('\   primary input source file: '$)
6002  Format('\ secondary input source file: '$)
6003  Format('\ association radius (arcsec):',F10.3)
6004  Format('\ unmatched sources from file: '$)
6005  Format(2I10)
c
c-----------------------------------------------------------------------
c
      End
c
      include 'gsagflds.f'
      include 'makenam.f'
c
c=======================================================================
c
      subroutine DelTmp(TmpNam)
c
      Character*500 TmpNam
      Integer*4 System, ISys
c
c-----------------------------------------------------------------------
c
      ISys = System('rm '//TmpNam)
      return
      end
c
c=======================================================================
c
      Subroutine MakeFmt(FmtLine,L)
c
      Character*11 FmtLine
      Character*4  Flen
      Integer*4    L
c
c-----------------------------------------------------------------------
c
      write(Flen,'(I4)') L
      FmtLine = '(A'//Flen//')'
      return
      end
c
c=======================================================================
c
      Subroutine FillIn(Null,Str,I1,I2)
c
      Character*5000 Null
      Character*64   Str
      Integer        I1, I2, N, LNBlnk, L
c
c-----------------------------------------------------------------------
c
10    L = LNBlnk(Str)
      if (Str(1:1) .eq. ' ') then
        do 20 N = 1, L-1
          Str(N:N) = Str(N+1:N+1)
20      continue
        Str(L:L) = ' '
        go to 10
      end if
c
      If (L .ge. I2-I1) then            ! too big to fit
        Null(I1:I1) = ' '
        do 30 N = 1, I2-I1
          Null(I1+N:I1+N) = Str(N:N)
30      continue
      else
        do 40 N = 1, I2-I1-L+1
          Null(N+I1-1:N+I1-1) = ' '
40      continue
        do 50 N = 1, L
          Null(N+I2-L:N+I2-L) = Str(N:N)
50      continue
      end if
c
      return
      end
c
c=======================================================================
c
c                                  from Numerical Recipes via T. Jarrett
      SUBROUTINE TJISORT(N,RA,ND)
c
      Integer*4 N,L,IR,J,I,ND(N),NDA
      Real*8 RA(N),RRA
c
      if (n .lt. 1) return
      Do 5 I = 1, N
        ND(I) = I
5     Continue
      if (n .lt. 2) return
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(ND(L))
          NDA = ND(L)
        ELSE
          RRA=RA(ND(IR))
          NDA = ND(IR)
          ND(IR)=ND(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            ND(1)=NDA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(ND(J)).LT.RA(ND(J+1)))J=J+1
          ENDIF
          IF(RRA.LT.RA(ND(J)))THEN
            ND(I)=ND(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        ND(I) = NDA
      GO TO 10
      END
c
c=======================================================================
c
      subroutine ClrCT(Str)
      character*(*) Str
      integer*4 I, lnblnk
c
      do 10 I = 1, lnblnk(Str)
        if (Str(I:I) .ne. '|') str(I:I) = ' '
10    continue
c
      return
      end
c
c=======================================================================
c
      subroutine GetIDB(Line,NBands,IDBand,OK)
c
      Character*5000 Line, TmpLin
      Character*1    NumChr
      Integer*4      NBands, IDBand(NBands), N, K, LNBlnk
      Logical        OK
c
c-----------------------------------------------------------------------
c                              Assumes band 1 <= band number <= 99
      N = 1
      TmpLin = Line
      OK = .True.
c
10    K = Index(TmpLin,'srcid')
      if (K .lt. 1) then
        print *, 'ERROR: can''t find srcid for ordinal band #', N,
     +           ' in header line:'
        print *,Line(1:LNBlnk(Line))
        OK = .False.
        return
      end if
c
      NumChr = TmpLin(K+5:K+5)
      if (NumChr .eq. '1') then
        IDBand(N) = 1
      else if (NumChr .eq. '2') then
        IDBand(N) = 2
      else if (NumChr .eq. '3') then
        IDBand(N) = 3
      else if (NumChr .eq. '4') then
        IDBand(N) = 4
      else if (NumChr .eq. '5') then
        IDBand(N) = 5
      else if (NumChr .eq. '6') then
        IDBand(N) = 6
      else if (NumChr .eq. '7') then
        IDBand(N) = 7
      else if (NumChr .eq. '8') then
        IDBand(N) = 8
      else if (NumChr .eq. '9') then
        IDBand(N) = 9
      else
        print *,'ERROR: non-numeric srcid suffix for ordinal band #', N,
     +          ' in header line:'
        print *,Line(1:LNBlnk(Line))
        OK = .False.
        return
      end if
c
      NumChr = TmpLin(K+6:K+6)
      if (NumChr .eq. '1') then
        IDBand(N) = 10*IDBand(N) + 1
      else if (NumChr .eq. '2') then
        IDBand(N) = 10*IDBand(N) + 2
      else if (NumChr .eq. '3') then
        IDBand(N) = 10*IDBand(N) + 3
      else if (NumChr .eq. '4') then
        IDBand(N) = 10*IDBand(N) + 4
      else if (NumChr .eq. '5') then
        IDBand(N) = 10*IDBand(N) + 5
      else if (NumChr .eq. '6') then
        IDBand(N) = 10*IDBand(N) + 6
      else if (NumChr .eq. '7') then
        IDBand(N) = 10*IDBand(N) + 7
      else if (NumChr .eq. '8') then
        IDBand(N) = 10*IDBand(N) + 8
      else if (NumChr .eq. '9') then
        IDBand(N) = 10*IDBand(N) + 9
      end if
c
      if (N .eq. NBands) return
c
      N = N + 1
      TmpLin(K:K) = '%'
      go to 10
c
      end
      include 'signon.f'
c
c=======================================================================
c
      subroutine GetTmpNam(TmpNam, TempDir,OutFNam, OK)
c      
      character*500 TmpNam0, TmpNam, WrkStrg, TempDir, OutFNam
      Character*11  Vsn
      Character*8   CDate, CTime
      integer*4     nTries, idum, k, lnblnk, Access, nCall, n
      real*4        ran1
      logical*4     OK
      data          nCall/0/
c      
      Common / VDT / CDate, CTime, Vsn
c
      if (nCall .eq. 0) then
        idum = -1
        nTries  = ran1(idum)        ! initialize random number generator
        nCall = nCall + 1
      end if
      TmpNam0 = TmpNam
      TmpNam  = OutFNam
      OK      = .true.
      nTries  = 0
c
1     n = index(TmpNam,'/')        ! for linux version
c1    n = index(TmpNam,'\')        ! for DOS version
      if (n .gt. 0) then           !
        do 2 k = 1, n              ! NOTE: there is still a loophole:
          TmpNam(k:k) = ' '        !       if two runs start simultaneously
2       continue                   !       in the same working/tempdir directory
        go to 1                    !       directory with output files named the
      end if                       !       the same but with different paths
      TmpNam  = AdjustL(TmpNam)
      TmpNam0 = TmpNam0(1:lnblnk(TmpNam0))//TmpNam
c     print *,'TmpNam0: ',tmpnam0(1:lnblnk(tmpnam0))   ! debug
c
10    nTries = nTries + 1
      if (nTries .gt. 999) then
        OK = .false.
        return
      end if
      k = 1.e6*ran1(idum)
      write (WrkStrg,'(i8)') k
20    k = index(WrkStrg,' ')
      if (k .lt. lnblnk(WrkStrg)) then 
        WrkStrg(k:k) = '%'
        go to 20
      end if      
      WrkStrg = WrkStrg(1:lnblnk(WrkStrg))//CDate//CTime
      k = index(WrkStrg,':')
      WrkStrg(k:k) = '%'
      k = index(WrkStrg,':')
      WrkStrg(k:k) = '%'
30    k = index(WrkStrg,' ')
      if (k .lt. lnblnk(WrkStrg)) then 
        WrkStrg(k:k) = '%'
        go to 30
      end if      
      TmpNam = TempDir(1:lnblnk(TempDir))//TmpNam0(1:lnblnk(TmpNam0))
     +         //WrkStrg(1:lnblnk(WrkStrg))
      if (Access(TmpNam(1:lnblnk(TmpNam)),' ') .eq. 0) go to 10
c     print *,TmpNam(1:lnblnk(TmpNam))     ! debug
      return
      end
c
c=======================================================================
c
      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL*4 ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END
c
c=======================================================================
c
      subroutine Setup(IUnitIn, IUnitOut, ColNam, ColTyp, TblFil,
     +                 RAstr, Decstr, IA1, IA2, ID1, ID2,
     +                 RAmax, RAmin, Decmax, Decmin, nColTyp,
     +                 Field, IFa, IFb, NF, ColTyp3, ColTyp4,
     +                 NLines, LRecL, FilNam, SkipErr, OK)
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
      real*8         RAmax, RAmin, Decmax, Decmin, RTmp
      integer*4      IUnitIn, IUnitOut, NLines, LRecL, IA1, IA2, ID1,
     +               ID2, LNBlnk, MaxLen, NTblHdr, N, IFa(MaxFld),
     +               IFb(MaxFld), NF, nColTyp
      Logical        TblFil, OK, GotCN, GotCT, SkipErr
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
         read (Line(IA1:IA2), *, err = 420) RTmp
         go to 430
420      if (.not.SkipErr) go to 3007         
430      if (RTmp .lt. RAmin) RAmin = RTmp
         if (RTmp .gt. RAmax) RAmax = RTmp
         read (Line(ID1:ID2), *, err = 440) RTmp
         go to 450
440      if (.not.SkipErr) go to 3008         
450      if (RTmp .lt. Decmin) Decmin = RTmp
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
c
c=======================================================================
c
      Function Best2Come(X1,Y1,Z1,
     +         X2,Y2,Z2,Indx1,N2,Window,
     +         NLines,N,ADist,Ang,
     +         mdetID2,CatWISE,LRecL,IM1,IM2,
     +         Coord1,Coord2)
c=======================================================================
c
      Character*1 Chr(5000)
      Character*5000 Line
      Integer*4 NLines(2), Indx1(NLines(1)), N, N1, N2, NN, mdetID1,
     +          mdetID2, LRecL, IM1,IM2, L
      Real*8    X1(NLines(1)), Y1(NLines(1)), Z1(NLines(1)), Ang,
     +          X2(NLines(2)), Y2(NLines(2)), Z2(NLines(2)), ADist,
     +          Cross, Ang2
      Real*8    Coord1(NLines(1)), Coord2(NLines(2)), Window
      Logical   Best2Come, CatWISE
      Equivalence   (Line, Chr(1))
c
c=======================================================================
c
      if (N .ge. NLines(1)) go to 1000
      do 100 NN = N+1, NLines(1)
        N1 = Indx1(NN)
        if (abs(Coord1(N1)-Coord2(N2)) .gt. ADist) go to 100
        if (abs(X1(N1)-X2(N2)) .gt. Window) go to 100
        if (abs(Y1(N1)-Y2(N2)) .gt. Window) go to 100
        if (abs(Z1(N1)-Z2(N2)) .gt. Window) go to 100
        if (CatWISE) then
          read(11, rec=N1, err = 1000) (Chr(L), L = 1, LRecl)
          read (Line(IM1:IM2), *, err = 1000) mdetID1
          if (mdetID1 .ne. mdetID2) go to 1000
        end if
        Cross = sqrt((Y1(N1)*Z2(N2) - Y2(N2)*Z1(N1))**2
     +             + (Z1(N1)*X2(N2) - Z2(N2)*X1(N1))**2
     +             + (X1(N1)*Y2(N2) - X2(N2)*Y1(N1))**2)
        Ang2 = dasin(Cross)
        if (Ang2 .lt. Ang) then
          Best2Come = .True.
          return
        end if
100   continue
c
1000  Best2Come = .False.
c
      return
c
      end
