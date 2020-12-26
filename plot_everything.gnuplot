set terminal postscript eps enhanced color solid 'Times' 24 linewidth 2
set output 'power.eps'
set xlabel 'Frequency (Hz)'
set ylabel 'Power'
# Power spectrum plot
plot 'power.dat' with lines linecolor 'royalblue' title ''
! convert  -density 150  'power.eps'   -background white -alpha remove  'power.png'
# Hm periodogram plot
set output 'Hm.eps'
set ylabel 'H_m'
plot 'Hm.dat' with lines linecolor 'forest-green' title ''
! convert  -density 150  'Hm.eps'   -background white -alpha remove  'Hm.png'
# Phased lightcurve plot
set output 'phase_folded_and_binned_lightcurve.eps'
set xlabel 'Phase'
set ylabel 'Photons'
plot \
"phase_folded_and_binned_lightcurve.dat" with boxes fill solid 1.0 linecolor 'red' title '', \
"phase_folded_and_binned_lightcurve.dat" with boxes fill empty linewidth 2.0 linecolor 'black' title ''
! convert  -density 150  'phase_folded_and_binned_lightcurve.eps'   -background white -alpha remove  'phase_folded_and_binned_lightcurve.png'
# Time lightcurve plot
set output 'binned_lightcurve_time.eps'
set xlabel 'Time (ks)'
set ylabel 'Count rate (cts/s)'
stats 'binned_lightcurve_time.dat' u ($1)/1000 name 'LCks' nooutput
plot \
[ LCks_min - (LCks_max-LCks_min)/20 : LCks_max + (LCks_max-LCks_min)/20] \
'binned_lightcurve_time.dat' u ($1)/1000:2:3 with errorbars linecolor 'dark-violet' linewidth 1.0 pointsize 1.0 pointtype 7 title ''
! convert  -density 150  'binned_lightcurve_time.eps'   -background white -alpha remove  'binned_lightcurve_time.png'
