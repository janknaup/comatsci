#!/bin/bash

#
# animate mueller-brown test plot, argument is frame count
#

echo "generating plot frames"
for ((i=0 ; i < $1 ; i++))
do
	index=`echo $i|awk '{printf "%04d",$1}'`
#	echo ${index}
	cp ranges.gplt temp.gplt
	echo "set term png transparent size 640,480; set out 'plot-${index}.png'; set view map; splot 'evol.dat' u 1:2:(0.) i ${index} w lp lt 2 t ''; set out " >> temp.gplt
	gnuplot < temp.gplt
	composite -compose src-over plot-${index}.png background.png gpacomposite${index}.png;
done

echo "composing gif animation from frames"

convert -delay 5 gpacomposite????.png gpamovie.gif
#ffmpeg -i gpacomposite%04d.png -r 25 -f avi -vcodec msmpeg4 -b 3M gpamovie.avi
