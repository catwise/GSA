      subroutine GetFlds(ColNam,Field,IFa,IFb,NF)
c-----------------------------------------------------------------------
c
c  Get fields in a table-file header line
c
c-----------------------------------------------------------------------
                     Integer*4  MaxFld
                     Parameter (MaxFld = 1000)
c
      character*5000 ColNam
      Character*300  Line
      character*25   Field(MaxFld)
      integer*4      IFa(MaxFld), IFb(MaxFld), NF, N, M, L, K, LNBlnk,
     +               LastErr
c
c-----------------------------------------------------------------------
c
      N = 0
      K = 0
      LastErr = 0
      do 100 M = 1, LNBlnk(ColNam)
        if (ColNam(M:M) .eq. '|') then
          N = N + 1
          NF = N - 1
          if (N .gt. 1) IFb(N-1) = M-1
          if (N .gt. MaxFld) return
          IFa(N) = M
          do 10 L = 1, 25
            Field(N)(L:L) = ' '
10        continue
          K = 0
        else
          if (ColNam(M:M) .ne. ' ') then
            K = K + 1
            if (K .le. 25) then
              Field(N)(K:K) = ColNam(M:M)
            else
              if (LastErr .ne. N) then
                write(Line,*) N
                Line = 'GetFlds - Table column name no. '
     +               //Line(1:lnblnk(Line))//' longer than 25 '
     +               //'characters: '//Field(N)//'....; excess ignored'
c               call ErrMsg(0,Line)
                print *,Line(1:lnblnk(line))
                LastErr = N
              end if
            end if
          end if
        end if
100   continue
c
      return
      end
