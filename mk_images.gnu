set term png size 1600,1200
do for [t=0:12]{
   #set output "image.".t.".png"
   if (t < 10) {
      set output "image.000".t.".png"}
   else{
      if (t < 100) {
         set output "image.00".t.".png"}
      else {
         if (t < 1000) {
            set output "image.0".t.".png"}
         else {
            if (t < 10000) {
               set output "image".t.".png"}}}}
   p [0:6.28] sin(x+t*3.14/50.0) w l ti 't='.t.' second'
   #plot [2:8] "data.txt" every :::1:1 using 1:3 with lines
   }

#set output "image.number.png"
##set yrange[]
#plot [2:8 "filename" u 1:3 w l
