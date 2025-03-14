setwd("C:/Users/franc/OneDrive/Desktop/Categorical Data Analysis")
library(openxlsx)
library(RefManageR)
library(bibtex)
library(ggplot2)
library(vcd)
library(package = binom)
library(package = PropCIs)
library(car)
library(dplyr)
library(emmeans)
library(knitr)
library(caret)
library(MASS)
library(reshape2)
library(caTools)
library(yardstick)

# Dataset e rinominazione di alcune variabili
pcos <- read.xlsx("PCOS_data_without_infertility.xlsx", sheet = 2)
pcos <- pcos %>%
  rename(
    Pimples = `Pimples(Y/N)`,
    Weight.gain = `Weight.gain(Y/N)`,
    PRG = `PRG(ng/mL)`,
    WHR = `Waist:Hip.Ratio`,
    AMH = `AMH(ng/mL)`,
    FSH_LH = `FSH/LH`,
    Cycle = `Cycle(R/I)`,
    Aborptions = `No..of.aborptions`,
    Hair.loss = `Hair.loss(Y/N)`,
    Hair.growth = `hair.growth(Y/N)`,
    FolliclesL = `Follicle.No..(L)`,
    FolliclesR = `Follicle.No..(R)`,
    Reg.Exercise = `Reg.Exercise(Y/N)`,
    Endometrium = `Endometrium.(mm)`,
    TSH = `TSH.(mIU/L)`,
    RBS = `RBS(mg/dl)`,
    Weight = `Weight.(Kg)`,
    Pulse.rate = `Pulse.rate(bpm)`,
    RR = `RR.(breaths/min)`,
    Hb = `Hb(g/dl)`,
    Cycle.length = `Cycle.length(days)`,
    Marriage.status = `Marraige.Status.(Yrs)`,
    Pregnant = `Pregnant(Y/N)`,
    I.beta.HCG = `I.beta-HCG(mIU/mL)`,
    II.beta.HCG = `II.beta-HCG(mIU/mL)`,
    FSH = `FSH(mIU/mL)`,
    LH = `LH(mIU/mL)`,
    Hip = `Hip(inch)`,
    Waist = `Waist(inch)`,
    Vit.D3 = `Vit.D3.(ng/mL)`,
    Skin.darkening = `Skin.darkening.(Y/N)`,
    Fast.food = `Fast.food.(Y/N)`,
    BP.systolic = `BP._Systolic.(mmHg)`,
    BP.diastolic = `BP._Diastolic.(mmHg)`,
    L.size = `Avg..F.size.(L).(mm)`,
    R.size = `Avg..F.size.(R).(mm)`,
    Height = `Height(Cm)`,
    PRL = `PRL(ng/mL)`
  )

# Modifiche su alcune variabili e gestione degli NA
pcos$FSH_LH <- 1/pcos$FSH_LH
pcos <- pcos %>%
  rename(LHFSH = `FSH_LH`)
pcos$Cycle[which(pcos$Cycle == 5)] <- 1
pcos$Follicles <- pmax(pcos$FolliclesL, pcos$FolliclesR)
pcos$F.size <- pmin(pcos$L.size, pcos$R.size)
pcos$AMH <- as.numeric(pcos$AMH)
pcos$AMH[is.na(pcos$AMH)] <- median(pcos$AMH[!is.infinite(pcos$AMH) & !is.na(pcos$AMH)], na.rm = TRUE)
pcos$LHFSH[is.infinite(pcos$LHFSH)] <- median(pcos$LHFSH[!is.infinite(pcos$LHFSH) & !is.na(pcos$LHFSH)], na.rm = TRUE)
pcos$II.beta.HCG <- as.numeric(pcos$II.beta.HCG)
pcos$II.beta.HCG[is.infinite(pcos$II.beta.HCG)] <- median(pcos$II.beta.HCG, na.rm = TRUE)
pcos$Pulse.rate[pcos$Pulse.rate <= 20] <- median(pcos$Pulse.rate, na.rm = TRUE)
pcos$Fast.food[is.na(pcos$Fast.food)] <- as.numeric(names(which.max(table(pcos$Fast.food))))
pcos$Marriage.status[is.na(pcos$Marriage.status)] <- median(pcos$Marriage.status, na.rm = TRUE)
pcos$II.beta.HCG[is.na(pcos$II.beta.HCG)] <- median(pcos$II.beta.HCG, na.rm = TRUE)
sum(is.na(pcos))
sum(is.infinite(as.matrix(pcos)))

# Rimozione di colonne non utili (id delle pazienti)
pcos <- pcos[, 3:46]

#####################################################################################

# RELAZIONE TRA PCOS E PESO CORPOREO

#####################################################################################


# Tabella di contingenza HighBMI e PCOS
pcos$HighBMI <- ifelse(pcos$BMI >= 25,  1, 0)
table <- table(HighBMI = pcos$HighBMI, PCOS = pcos$PCOS)
table <- table[, c(2, 1)]
pi.hat.table <- table / rowSums(table)
addmargins(table)

# Tabella di contingenza con le proporzioni 
pi.hat.table

# Intervalli di confindenza per le stime delle proporzioni 
binom.confint(x = table[1,1], n = 314, conf.level = 0.95, methods = c("asymptotic"))
binom.confint(x = table[2,1], n = 227, conf.level = 0.95, methods = c("asymptotic"))

# Differenza tra proporzioni e test di ipotesi 
wald2ci(x1 = table[1,1], n1 = sum(table[1,]), x2 = table[2,1], n2 = sum(table[2,]), conf.level = 0.95, adjust = "Wald")
prop.test(x = table[,1], n = rowSums(table), conf.level = 0.95, correct = FALSE)

# Stima del rischio relativo e intervallo di confidenza
alpha <- 0.05
n1 <- as.integer(sum(table[1,]))
n2 <- as.integer(sum(table[2,]))
pi.hat1 <- pi.hat.table[1, 1]
pi.hat2 <- pi.hat.table[2, 1]
ci <- exp(log(pi.hat1/pi.hat2) + qnorm(p = c(alpha/2, 1-alpha/2)) * sqrt((1-pi.hat1)/(n1*pi.hat1) + (1- pi.hat2)/(n2*pi.hat2)))
rev(round(1/ci, 4))

# Stima dell'odds ratio e intervallo di confidenza
OR.hat <- (table[1,1] * table[2,2]) / (table[2,1] * table[1,2])
var.log.or <- 1/table[1,1] + 1/table[1,2] + 1/table[2,1] + 1/table[2,2]
OR.CI <- exp(log(OR.hat) + qnorm(p = c(alpha/2, 1-alpha/2)) * sqrt(var.log.or))
round(1/OR.hat, 2)
rev(round(1/OR.CI, 2))

# Grafico della proporzione di PCOS per categoria di BMI
pcos$BMIlevel <- cut(pcos$BMI, breaks = c(-Inf, 18.5, 25, 30, Inf), labels = c("Underweight", "Normal Weight", "Overweight", "Obese"), right = FALSE)
prop_table <- as.data.frame(prop.table(table(pcos$BMIlevel, pcos$PCOS), margin = 1))
colnames(prop_table) <- c("BMIlevel", "PCOS", "Proportion")
ggplot(prop_table, aes(x = BMIlevel, y = Proportion, fill = factor(PCOS, labels = c("No PCOS", "PCOS")))) +
  geom_bar(stat = "identity", position = "fill", alpha = 0.8) +
  scale_fill_manual(values = c("#69b3a2", "#404080")) +
  labs(title = "Proporzione di PCOS per categoria di BMI", x = "Categoria di BMI", y = "Proporzione", fill = "Diagnosi") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Modello di regressione logistica PCOS ~ BMIlevel
pcos$BMIlevel <- relevel(x = pcos$BMIlevel, ref = "Normal Weight")
mod.fit <- glm(formula = PCOS ~ BMIlevel, family = binomial(link = logit), data = pcos)
summary(mod.fit)
calc.est <- emmeans(object = mod.fit, specs = ~ BMIlevel , type = "response")
confint(object = contrast(object = calc.est, method = "revpairwise"), adjust = "none", level = 0.95)

#####################################################################################

# VERIFICA DELL'INDIPENDENZA TRA HighLHFSH E BMIlevel

#####################################################################################


# Creazione della variabile HighLHFSH e sottoinsieme per sole pazienti con PCOS
pcos$HighLHFSH <- ifelse(pcos$LHFSH >= 2,  1, 0) 
pcos_subset <- subset(pcos, PCOS == 1)

# Tabella di contingenza HighFHLSH e BMIlevel
table1 <- table(HighFHLSH = pcos_subset$HighLHFSH, BMIlevel = pcos_subset$BMIlevel)
ind.test <- chisq.test(x = table1, correct = FALSE)
addmargins(table1)

# Tabella con le proporzioni
table1 / sum(table1)

# Test di indipendenza tra HighFHLSH e BMIlevel
assocstats(x = table1)

# Numero atteso di conteggi
ind.test$expected

#####################################################################################

# SVILUPPO DI UN CLASSIFICATORE PER IL TIPO DI PCOS

#####################################################################################


# Costruzione della variabile che indica il tipo di PCOS
pcos$PCOStype <- NA
pcos$PCOStype[pcos$PCOS == 0] <- "None"
pcos$PCOStype[pcos$PCOS == 1 & pcos$AMH >= 5] <- "Anovulatory"
pcos$PCOStype[pcos$PCOS == 1 & pcos$AMH < 5] <- "Ovulatory"
table(pcos$PCOStype)

# Rimozione di PCOS
pcos2 <- pcos[, -c(1)]
# Rimozione di alcune variabili 
pcos2 <- pcos2[, -c(2, 3, 11, 16, 17, 19, 20, 37, 38, 39, 40)]
pcos2$PCOStype <- as.factor(pcos2$PCOStype)
pcos2$PCOStype <- relevel(x = pcos2$PCOStype, ref = "None")
# Rimozione di AMH (ridondante rispetto alla variabile target)
pcos2 <- pcos2[, -c(16)]
# Rimozione di alcune variabili precedentemente aggiunte per le altre analisi
pcos2 <- pcos2[, -c(32, 33, 34)]

# Divisione in training e test set (80%-20%)
set.seed(23)
train_idx <- createDataPartition(pcos2$PCOStype, p = 0.8, list = FALSE)
train_data <- pcos2[train_idx, ]
test_data <- pcos2[-train_idx, ]

# Proporzioni delle classi di PCOStype nel training set e nel test set
cat("Proporzioni nel training set:\n")
print(prop.table(table(train_data$PCOStype)))
cat("\nProporzioni nel test set:\n")
print(prop.table(table(test_data$PCOStype)))

# Forward selection
# Costruzione del modello vuoto e del modello pieno
empty.mod <- multinom(formula = PCOStype ~ 1, data = train_data) 
full.mod <- multinom(formula = PCOStype ~ ., data = train_data)

forw.sel2 <- step(object = empty.mod, scope = list(upper = full.mod), direction = "forward", k = log(nrow(train_data)), trace = TRUE)


model <- multinom(formula(forw.sel2), data = train_data, trace = FALSE)
model

# Test globale per i parametri stimati
Anova(model)

# Accuratezza sul test set
preds <- predict(model, newdata = test_data)
accuracy <- mean(preds == test_data$PCOStype)
cat("\nAccuratezza finale sul test set:", accuracy)

# Metriche di valutazione
confusionMatrix(preds, test_data$PCOStype, mode = "everything")

# Calcolo della F1-score micro
test_results <- data.frame(
  truth = test_data$PCOStype, 
  prediction = preds
)
test_results$truth <- as.factor(test_results$truth)
test_results$prediction <- as.factor(test_results$prediction)
f_meas(test_results, truth, prediction, estimator = "micro")


# Matrice di confusione
conf_matrix <- table(Predicted = preds, Actual = test_data$PCOStype)
conf_matrix_df <- as.data.frame(as.table(conf_matrix))
ggplot(data = conf_matrix_df, aes(x = Predicted, y = Actual, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  geom_text(aes(label = Freq), color = "white", size = 5) +
  theme_minimal() +
  labs(title = "Matrice di confusione", x = "Classe predetta", y = "Classe reale") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#####################################################################################

# MODELLO DI REGRESSIONE DI POISSON PER IL NUMERO DI FOLLICOLI

#####################################################################################

# Modello di regressione di Poisson
pcos_subset <- subset(pcos, PCOS == 1)
mean(pcos_subset$Follicles)
var(pcos_subset$Follicles)

# Distribuzione empirica e teorica di Follicles
rel.freq <- table(pcos_subset$Follicles) / length(pcos_subset$Follicles)
rel.freq2 <- c(rel.freq, rep(0, times = 4))
y <- 1:25
prob <- round(dpois(x = y, lambda = mean(pcos_subset$Follicles)), 4)
data.frame(y, prob, rel.freq = rel.freq2)

plot(x = y-0.1, y = prob, type = "h", ylab = "Probabilità", xlab = "Numero di follicoli", lwd = 2, xaxt = "n")
axis(side = 1, at = 0:41)
lines(x = y+0.1, y = rel.freq2, type = "h", lwd = 2, lty = "solid", col = "red")
abline(h = 0)
legend(x = 17, y = 0.08, legend = c("Poisson", "Dati osservati"), lty = c("solid", "solid"), lwd = c(2,2), 
       col = c("black", "red"), bty = "n", cex = 0.5)


# Stima del modello di regressione di Poisson
mod.fit <- glm(formula = Follicles ~ LHFSH, data = pcos_subset, family = poisson(link = log))
summary(mod.fit)
mu.hat <- mod.fit$fitted.values

# Plot dei residui standardizzati di Pearson per il modello di regressione di Poisson
stand.resid <- rstandard(model = mod.fit, type = 
                             "pearson")  
plot(x = mu.hat, y = stand.resid, xlab = 
         expression(hat(mu)), ylab = "Standardized Pearson 
    residuals", ylim = c(min(c(-3, stand.resid)), max(c(3, 
                                                        stand.resid))))
abline(h = c(-3,-2,0,2,3), lty = "dotted", col = "red")

sum(abs(stand.resid) > 3)

# Intervallo di confidenza per i parametri 
sd(pcos_subset$LHFSH)
100 * (exp(sd(pcos_subset$LHFSH) * mod.fit$coefficients[2]) - 1)

beta.ci <- confint(object = mod.fit, parm = "LHFSH", level = 0.95)
100 * (exp(sd(pcos_subset$LHFSH) * beta.ci) - 1)

# Funzione per calcolare gli intervalli di confidenza
ci.mu <- function(newdata, mod.fit.obj, alpha) {
  lin.pred.hat <- predict(object = mod.fit.obj, newdata = newdata, type = "link", se = TRUE)
  lower <- exp(lin.pred.hat$fit - qnorm(1 - alpha / 2) * lin.pred.hat$se)
  upper <- exp(lin.pred.hat$fit + qnorm(1 - alpha / 2) * lin.pred.hat$se)
  list(lower = lower, upper = upper)
}

# Suddivisione dei valori di LHFSH in 8 intervalli
# Calcolo delle medie in ogni intervallo
cutoff1 <- quantile(pcos_subset$LHFSH, probs = 0:8 / 8, na.rm = FALSE)
groups2_1 <- cut(x = pcos_subset$LHFSH, cutoff1)
ybar1 <- aggregate(x = Follicles ~ groups2_1, data = pcos_subset, FUN = mean)
xbar1 <- aggregate(x = LHFSH ~ groups2_1, data = pcos_subset, FUN = mean)
count1 <- aggregate(x = Follicles ~ groups2_1, data = pcos_subset, FUN = length)
data.frame(ybar1, mean.LHFSH = xbar1$LHFSH, count = count1$Follicles)

# Grafico dell'andamento del numero di follicoli rispetto al rapporto LH:FSH 
# secondo il modello stimato, con intervalli di confidenza
plot(x = pcos_subset$LHFSH, y = pcos_subset$Follicles, xlab = "LHFSH", 
     ylab = "Numero di follicoli",
     main = "Modello di regressione di Poisson per il numero di follicoli", 
     panel.first = grid(), cex.main = 0.8)
curve(expr = exp(mod.fit$coefficients[1] + mod.fit$coefficients[2] * x), 
      col = "red", add = TRUE, lty = "solid")
curve(expr = ci.mu(newdata = data.frame(LHFSH = x), 
                   mod.fit.obj = mod.fit, alpha = 0.05)$lower,
                   col = "blue", add = TRUE, lty = "dotdash")
curve(expr = ci.mu(newdata = data.frame(LHFSH = x), 
                   mod.fit.obj = mod.fit, alpha = 0.05)$upper,
                   col = "blue", add = TRUE, lty = "dotdash")
# Aggiunta delle medie per ogni intervallo
points(x = xbar1$LHFSH, y = ybar1$Follicles, pch = 17, col = "darkgreen", cex = 1)

# Modello di regressione quasi-Poisson
mod.fit.quasi <- glm(formula = Follicles ~ LHFSH, data = 
                       pcos_subset, family = quasipoisson(link = log))
summary(mod.fit.quasi)

mu.hat.quasi <- mod.fit.quasi$fitted.values

# Plot dei residui standardizzati di Pearson per il modello quasi-Poisson
stand.resid.quasi <- rstandard(model = mod.fit.quasi, type = 
                           "pearson")  

plot(x = mu.hat.quasi, y = stand.resid.quasi, xlab = 
       expression(hat(mu)), ylab = "Standardized Pearson 
    residuals", ylim = c(min(c(-3, stand.resid.quasi)), max(c(3, stand.resid.quasi))))
abline(h = c(-3,-2,0,2,3), lty = "dotted", col = "red")

sum(abs(stand.resid.quasi) > 3) 

# Funzione per calcolare gli intervalli di confidenza per un modello quasi-Poisson
ci.mu.quasi <- function(newdata, mod.fit.obj, alpha, dispersion) {
  lin.pred.hat <- predict(object = mod.fit.obj, newdata = newdata, type = "link", se = TRUE)
  lower <- exp(lin.pred.hat$fit - qnorm(1 - alpha / 2) * lin.pred.hat$se * sqrt(dispersion))
  upper <- exp(lin.pred.hat$fit + qnorm(1 - alpha / 2) * lin.pred.hat$se * sqrt(dispersion))
  list(lower = lower, upper = upper)
}

# Grafico del modello quasi-Poisson con intervalli di confidenza 
plot(x = pcos_subset$LHFSH, y = pcos_subset$Follicles, xlab = "LHFSH", 
     ylab = "Numero di follicoli",
     main = "Modello di regressione quasi-Poisson per il numero di follicoli", 
     panel.first = grid(), cex.main = 0.8)
curve(expr = exp(mod.fit.quasi$coefficients[1] + mod.fit.quasi$coefficients[2] * x), 
      col = "red", add = TRUE, lty = "solid")
curve(expr = ci.mu.quasi(newdata = data.frame(LHFSH = x), 
                         mod.fit.obj = mod.fit.quasi, alpha = 0.05, dispersion = summary(mod.fit.quasi)$dispersion)$lower,
      col = "blue", add = TRUE, lty = "dotdash")
curve(expr = ci.mu.quasi(newdata = data.frame(LHFSH = x), 
                         mod.fit.obj = mod.fit.quasi, alpha = 0.05, dispersion = summary(mod.fit.quasi)$dispersion)$upper,
      col = "blue", add = TRUE, lty = "dotdash")
# Aggiunta delle medie per ogni intervallo
points(x = xbar1$LHFSH, y = ybar1$Follicles, pch = 17, col = "darkgreen", cex = 1)




