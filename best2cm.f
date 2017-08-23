      Function Best2Come(X1,Y1,Z1,
     +         X2,Y2,Z2,Indx1,N2,Window,
     +         NLines,N,ADist,Ang,
     +         Coord1,Coord2)
c=======================================================================
c
      Integer*4 NLines(2), Indx1(NLines(1)), N, N1, N2, NN
      Real*8    X1(NLines(1)), Y1(NLines(1)), Z1(NLines(1)), Ang,
     +          X2(NLines(2)), Y2(NLines(2)), Z2(NLines(2)), ADist,
     +          Cross, Ang2
      Real*4    Coord1(NLines(1)), Coord2(NLines(2)), Window
      Logical   Best2Come
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
