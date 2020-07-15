###Introduction Exercises

library(UsingR)
data("father.son",package="UsingR")

#1
mean(father.son$sheight)

#2
mean(father.son$sheight[round(father.son$fheight) == 71])

#3
#Which of the following can't be written as a linear model?
  #Y = a + b^t + e

#4
#Which of the following best describes what  e  represents?
  #Between-individual variability: people of the same weight vary in their height.