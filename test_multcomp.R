library(multcomp)
hsb2 <- read.csv("http://www.ats.ucla.edu/stat/data/hsb2.csv")

hsb2$ses <- factor(hsb2$ses)
hsb2$female <- factor(hsb2$female)


m1 <- lm(read ~ socst + ses * female, data = hsb2)
summary(m1)


model.matrix(read ~ socst + ses * female, data = hsb2)

# difference between ses = 2 and ses =3 when female = 0
K <- matrix(c(0, 0, 1, -1, 0, 0, 0), 1)
t <- glht(m1, linfct = K)
summary(t)



# difference between ses = 2 and ses =3 when female = 1
K <- matrix(c(0, 0, 1, -1, 0, 1, -1), 1)
t <- glht(m1, linfct = K)
summary(t)




