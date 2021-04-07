~~~
set tics font "Helvetica,20"
unset border
unset xtics
unset ytics
unset legend
unset key
set xrange [0:2*pi]; plot(cos(2*x-pi/2)+cos(x-pi/4)) lw 12 linecolor rgb "black"
set yrange [0:5]
plot [-4:4] 1/abs(x) lw 4, erf(x)/x lw 3 dt 2 , (1-erf(abs(x)))/abs(x) lw 3 dt 2
~~~
