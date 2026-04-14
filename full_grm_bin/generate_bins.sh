awk 'BEGIN {                 
  for (x = -0.05; x < 0.02; x += 0.001)
    printf "%.3f %.3f\n", x, x + 0.001;
  for (x = 0.02; x < 0.2; x += 0.005)
    printf "%.3f %.3f\n", x, x + 0.005;
  for (x = 0.2; x < 1.5; x += 0.01)
    printf "%.3f %.3f\n", x, x + 0.01;
}' > bins.txt