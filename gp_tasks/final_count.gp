count = 990
step = 10

set term gif animate delay 10 size 2000, 1000
set yrange [0:1.5]

set output "gifs/final_count.gif"

do for [n = 0 : count : step] {
    filename = sprintf("data/final_count/total_count.dat, n)
    plot filename with lines title ""
}
