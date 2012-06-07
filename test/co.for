        program cod
        double precision  co(10,21),pi
        integer lmax,l,al, am, m

        pi=4.*datan(1.d0)

        lmax=10
        do 10 l=1,lmax
            al=l
            co(l,l+1)=dsqrt((2.*al+1.)/2./pi)
            do 20 m=1,l-1
                am=m
                co(l,l+1+m)=co(l,l+m)*dsqrt(1.d0/(al-am+1)/(al+am))
                co(l,l+1-m)=co(l,l+1+m)
20          continue
            co(l,2*l+1)=co(l,2*l)*dsqrt(1.d0/(2.*al))
            co(l,1)=co(l,2*l+1)
            co(l,l+1)=dsqrt((2.*al+1.)/4./pi)
10      continue

        do 30 l=1,lmax
            do 40 m=1,2*l+1
                print *, 'co(',l,',',m,')=',co(l,m)
40          continue
30      continue

        end program cod

