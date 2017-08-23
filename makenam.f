      subroutine MakeNam(NamDes,NamGot,Line,I01,I2,NamAvd,NumNA,MT)
c-----------------------------------------------------------------------
c
c   Make a name as much like NamDes as possible that fits into a space
c   from I01+1 to I2 inclusive and avoids duplicating any name in the
c   array NamAvd, which contains NumNA entries; if MT <> 0 insert the
c   name into the string Line at the location indicated, centered as
c   much as possible. If space permits, add a suffix '2' and increment
c   as needed through 'z', repeat as needed until space runs out. If no
c   space for suffix, increment characters, starting on the right and
c   moving left as needed. ONLY INCREMENT characters; this avoids
c   duplicating a previously accepted name that avoids a previously
c   processed list. Accept failure if necessary, but warn.
c
c   If MT = 0, do not insert the name if no change was necessary; if a
c   change was necessary, clear the field range before inserting.
c
c-----------------------------------------------------------------------
c
                     Integer*4  MaxFld
                     Parameter (MaxFld = 1000)
c
      character*5000 Line
      character*25   NamDes, NamGot, NamAvd(MaxFld)
      character*1    Chr1
      integer*4      I01, I1, I2, NumNA, N, LNB, LNBlnk, NTries, NCols,
     +               NTmax, NPos, Nslack, MT
      byte           IChr
      logical        Tight, OK2Add
c
      data           NTmax/1000/
c
      equivalence   (Chr1, IChr)
c
c-----------------------------------------------------------------------
c
      NamGot = NamDes
      LNB = LNBlnk(NamGot)
      NPos = LNB                             ! place to increment
      OK2Add = .True.
      I1 = I01 + 1                           ! I01 : pipe column
      NTries = 0
      NCols = I2-I01                         ! Name must fit in columns
      if (LNBlnk(NamGot) .gt. (NCols)) then  !  I1-I2 inclusive
        do 10 N = NCols+1, 25                ! Trim iff necessary
          NamGot(N:N) = ' '
10      continue
        LNB = NCols
      end if
c
20    Tight = LNB .eq. NCols                 ! Tight: no room for suffix
      do 100 n = 1, NumNA
        if (NamGot .eq. NamAvd(n)) then
          NTries = NTries + 1
          if (NTries .gt. NTmax) go to 90
          if (Tight) then                    ! no suffix possible, must
30          Chr1 = NamGot(NPos:NPos)         !  increment a character
            if (IChr .ge. 122) then          ! 'z', can't increment
              if (NPos .eq. 1) then          ! can't move left; bail
                NamGot = NamDes              !  out with warning
                LNB = LNBlnk(NamGot)
                go to 90
              else
                NPos = NPos - 1              ! move left and try again
                go to 30
              end if
            else
              IChr = IChr + 1                ! increment character at
              if ((IChr .ge. 58) .and.       !  NPos, avoid ':' - '`'
     +            (IChr .le. 96))  IChr = 97
              NamGot(NPos:NPos) = Chr1
              go to 20
            end if
          else           ! not Tight
40          if (OK2Add) then
              OK2Add = .False.               ! append a suffix,
              LNB = LNB + 1                  !  starting with '2'
              NamGot(LNB:LNB) = '2'
              NPos = LNB
              go to 20
            else
50            Chr1 = NamGot(NPos:NPos)       !  increment a character
              if (IChr .ge. 122) then        ! 'z', can't increment
                if (NPos .eq. 1) then        ! can't move left; add
                  OK2Add = .True.            !  another suffix
                  go to 40
                else
                  NPos = NPos - 1            ! move left and try again
                  go to 50
                end if
              else
                IChr = IChr + 1              ! increment character at
                if ((IChr .ge. 58) .and.     !  NPos, avoid ':' - '`'
     +              (IChr .le. 96))  IChr = 97
                NamGot(NPos:NPos) = Chr1
                go to 20
              end if
            end if
          end if
90        print *,
     +     'WARNING: unable to avoid duplicating the table-file'
          print *,
     +     '         header column name ',NamGot(1:LNB)
          go to 1000
        end if           ! NamGot .eq. NamAvd(n)
100   continue
      if ((NTries .eq. 0) .and. (MT .eq. 0))  ! return if no conflict
     +     return                             ! MT = 0 --> no insert
      go to 1000                              ! in such a case
c
1000  if (Tight) then
        do 1010 n = 1, LNB
          Line(I01+n:I01+n) = NamGot(n:n)
1010    continue
      else
        if (MT .eq. 0) then                   ! blank out the field
          do 1020 n = I1, I2                  ! before inserting
            Line(n:n) = ' '
1020      continue
        end if
        Nslack = (NCols - LNB)/2
        do 1025 n = 1, LNB
          Line(I01+Nslack+n:I01+Nslack+n) = NamGot(n:n)
1025    continue
      end if
c
      return
      end
