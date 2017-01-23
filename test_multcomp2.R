
library("ISwR")

thuesen.lm <- lm(short.velocity ~ blood.glucose, data = thuesen)
summary(thuesen.lm)


library("multcomp")
thuesen.mc <- glht(thuesen.lm, linfct = diag(2))
summary(thuesen.mc, test = adjusted(type = "bonferroni"))

summary(thuesen.mc)

C <- diag(2)

thuesen.mc <- glht(thuesen.lm, linfct = C)
summary(thuesen.mc)

summary(thuesen.mc, test = adjusted(type = "Westfall"))


# ----------------------------------
# LM - testing one coefficient - lm and glht p-values for one coeff are the same

thuesen.lm <- lm(short.velocity ~ blood.glucose, data = thuesen)
summary(thuesen.lm)

thuesen.mc <- glht(thuesen.lm, linfct = matrix(c(1, 0), 1))
summary(thuesen.mc)

# ----------------------------------
# GLM norm - testing one coefficient - glm and glht p-values for one coeff are DIFFERENT!!!

thuesen.glm <- glm(short.velocity ~ blood.glucose, data = thuesen)
summary(thuesen.glm)

thuesen.gmc <- glht(thuesen.glm, linfct = matrix(c(1, 0), 1))
summary(thuesen.gmc)

# ----------------------------------
# GLM logit - testing one coefficient - glm and glht p-values for one coeff are the same

library("coin")

data("alzheimer", package = "coin")

y <- factor(alzheimer$disease == "Alzheimer", labels = c("other", "Alzheimer"))

alzheimer.glm <- glm(y ~ smoking * gender, data = alzheimer, family = binomial())
summary(alzheimer.glm)

alzheimer.gmc <- glht(alzheimer.glm, linfct = matrix(c(0, 0, 0, 1, 0, 0, 0, 0), 1))
summary(alzheimer.gmc)

































