library(timeSeries)
library(fBasics)
library(tseries)
library(fracdiff)
library(FinTS)
library(fGarch)
library(TSA)
library(timeDate)
library(rugarch)
library(forecast)
library(moments)
library(lmtest)
library(Metrics)


## SETUP 

setwd(...)
load(file="dati_giugno_2024.csv")
y <- ts(Y)

## VISUALIZZAZIONE DELLA SERIE E DELLE CORRELAZIONI

ts.plot(y, main = expression(paste("Serie ", Y[t])))

acf(y, lag.max=50, main = expression(paste("ACF della serie ", 
                                           Y[t])))
lag.plot(y, lags=50, main = expression(paste("Lag plot della serie ", 
                                             Y[t])))

y <- ts(y, frequency=4)

## DECOMPOSIZIONE TRADIZIONALE

mydata=decompose(y, "multiplicative")
plot(mydata)


## ANALISI DELLA STAZIONARIETÀ

adf.test(y, alternative = "stationary")

y<-ts(y)

ydiff <- diff(y, differences=1)

ts.plot(ydiff, main = expression(paste("Serie ", 
                                nabla, " ", Y[t])))
adf.test(ydiff, alternative = "stationary")
acf(ydiff, lag.max=50, main = expression(paste("ACF della serie ", 
                                nabla, " ", Y[t])))

## ANALISI DELLA STAGIONALITÀ

y_seasonal_diff <- diff(ydiff, differences = 1, lag=4)


ts.plot(y_seasonal_diff, main = expression(paste("Serie ", 
                                                 nabla[4], nabla, " ", Y[t])))
acf(y_seasonal_diff, lag.max = 50, main = expression(paste("ACF della serie ", 
                                                           nabla[4], nabla, " ", Y[t])))



## STIMA DEI MODELLI ARMA

auto.arima(y_seasonal_diff, d=0, D=0, max.p = 3, max.q = 3, max.P=0, max.Q=0, stationary=TRUE, seasonal=FALSE,
           allowmean=FALSE)


modello1 <- Arima(y_seasonal_diff, order = c(1, 0, 0), include.mean = FALSE)
modello2 <- Arima(y_seasonal_diff, order = c(0, 0, 1), include.mean = FALSE)
modello3 <- Arima(y_seasonal_diff, order = c(1, 0, 1), include.mean = FALSE)
modello4 <- Arima(y_seasonal_diff, order = c(2, 0, 1), include.mean = FALSE)
modello5 <- Arima(y_seasonal_diff, order = c(1, 0, 2), include.mean = FALSE)
modello6 <- Arima(y_seasonal_diff, order = c(2, 0, 2), include.mean = FALSE)
modello7 <- Arima(y_seasonal_diff, order = c(3, 0, 0), include.mean = FALSE)

## RISULTATI DELLE STIME 

df <- data.frame(
  Modello = c("1", "2", "3", "4", "5", "6", "7"),
  AIC = c("345.609855844387", "327.679171975387", "329.654993753184", "327.507636713804", "331.620722474474", "322.639426602361", "318.779173458529"),
  AICc = c("345.71699870153", "327.78631483253", "329.8712099694", "327.87127307744", "331.98435883811", "323.189885317957", "319.142809822166"),
  BIC = c("351.099720101113", "333.169036232114", "337.889790138273", "338.487365227257", "342.600450987927", "336.364087244177", "329.758901971982"),
  Significativi = c("sì", "sì", "no", "no", "no", "no", "sì"),
  stringsAsFactors = FALSE
)
etichette <- c("ARMA(1,0)", "ARMA(0,1)", "ARMA(1,1)", "ARMA(2,1)", "ARMA(1,2)", "ARMA(2,2)", "ARMA(3,0)")
df$Modello <- etichette

df$AIC <- as.numeric(df$AIC)
df$AICc <- as.numeric(df$AICc)
df$BIC <- as.numeric(df$BIC)
df$AIC <- round(df$AIC, 2)
df$AICc <- round(df$AICc, 2)
df$BIC <- round(df$BIC, 2)

print(df)

par(mfrow=c(4,2))

plot(modello1$residuals, type="l", ylab="Residui", main = "Residui di ARMA(1,0)")
plot(modello2$residuals, type="l", ylab="Residui", main = "Residui di ARMA(0,1)")
plot(modello3$residuals, type="l", ylab="Residui", main = "Residui di ARMA(1,1)")
plot(modello4$residuals, type="l", ylab="Residui", main = "Residui di ARMA(2,0)")
plot(modello5$residuals, type="l", ylab="Residui", main = "Residui di ARMA(1,2)")
plot(modello6$residuals, type="l", ylab="Residui", main = "Residui di ARMA(2,2)")
plot(modello7$residuals, type="l", ylab="Residui", main = "Residui di ARMA(3,0)")


Boxtest1<-Box.test(modello1$residuals, lag=30, type = "Box-Pierce", fitdf=1)
Boxtest2<-Box.test(modello2$residuals, lag=30, type = "Box-Pierce", fitdf=1)
Boxtest3<-Box.test(modello3$residuals, lag=30, type = "Box-Pierce", fitdf=2)
Boxtest4<-Box.test(modello4$residuals, lag=30, type = "Box-Pierce", fitdf=3)
Boxtest5<-Box.test(modello5$residuals, lag=30, type = "Box-Pierce", fitdf=3)
Boxtest6<-Box.test(modello6$residuals, lag=30, type = "Box-Pierce", fitdf=4)
Boxtest7<-Box.test(modello7$residuals, lag=30, type = "Box-Pierce", fitdf=3)

Ljungtest1<-Box.test(modello1$residuals, lag=30, type = "Ljung-Box", fitdf=1)
Ljungtest2<-Box.test(modello2$residuals, lag=30, type = "Ljung-Box", fitdf=1)
Ljungtest3<-Box.test(modello3$residuals, lag=30, type = "Ljung-Box", fitdf=2)
Ljungtest4<-Box.test(modello4$residuals, lag=30, type = "Ljung-Box", fitdf=3)
Ljungtest5<-Box.test(modello5$residuals, lag=30, type = "Ljung-Box", fitdf=3)
Ljungtest6<-Box.test(modello6$residuals, lag=30, type = "Ljung-Box", fitdf=4)
Ljungtest7<-Box.test(modello7$residuals, lag=30, type = "Ljung-Box", fitdf=3)

BoxPierce <- c(Boxtest1$p.value, Boxtest2$p.value, Boxtest3$p.value, Boxtest4$p.value, Boxtest5$p.value, 
               Boxtest6$p.value, Boxtest7$p.value)
LjungBox <- c(Ljungtest1$p.value, Ljungtest2$p.value, Ljungtest3$p.value, Ljungtest4$p.value, Ljungtest5$p.value,
              Ljungtest6$p.value, Ljungtest7$p.value)
Incorrelati <- c("no", "no", "no", "no", "no", "no", "no")

results <- data.frame(Model = etichette, BoxPierce = BoxPierce, LjungBox = LjungBox, Incorrelati)

print(results)

## RESIDUI DELLE STIME

par(mfrow=c(1,1))
par(mfrow=c(4,2))
acf(modello1$residuals, main = "ACF dei residui di ARMA(1,0)")
acf(modello2$residuals, main = "ACF dei residui di ARMA(0,1)")
acf(modello3$residuals, main = "ACF dei residui di ARMA(1,1)")
acf(modello4$residuals, main = "ACF dei residui di ARMA(2,1)")
acf(modello5$residuals, main = "ACF dei residui di ARMA(1,2)")
acf(modello6$residuals, main = "ACF dei residui di ARMA(2,2)")
acf(modello7$residuals, main = "ACF dei residui di ARMA(3,0)")

par(mfrow=c(4,2))
hist(modello1$residuals, breaks = 10,  xlab="", ylab = "Frequenze relative",
     ylim = c(0,0.4), xlim = c(-12,12), freq = F, main="Residui di ARMA(1,0)", cex.main=0.5)
curve(dnorm(x, mean(modello1$residuals), sd(modello1$residuals)), col = "red", lwd = 3, add = T)

hist(modello2$residuals, breaks = 10, xlab="", ylab = "Frequenze relative",
     ylim = c(0,0.4), xlim = c(-12,12), freq = F, main="Residui di ARMA(0,1)", cex.main=0.5)
curve(dnorm(x, mean(modello2$residuals), sd(modello2$residuals)), col = "red", lwd = 3, add = T)


hist(modello3$residuals, breaks = 10, xlab="", ylab = "Frequenze relative",
     ylim = c(0,0.4), xlim = c(-12,12), freq = F, main="Residui di ARMA(1,1)", cex.main=0.5)
curve(dnorm(x, mean(modello3$residuals), sd(modello3$residuals)), col = "red", lwd = 3, add = T)


hist(modello4$residuals, breaks = 10, xlab="", ylab = "Frequenze relative",
     ylim = c(0,0.4), xlim = c(-12,12), freq = F, main="Residui di ARMA(2,1)", cex.main=0.5)
curve(dnorm(x, mean(modello4$residuals), sd(modello4$residuals)), col = "red", lwd = 3, add = T)

hist(modello5$residuals, breaks = 10,  xlab="", ylab = "Frequenze relative",
     ylim = c(0,0.4), xlim = c(-12,12), freq = F, main="Residui di ARMA(1,2)", cex.main=0.5)
curve(dnorm(x, mean(modello5$residuals), sd(modello5$residuals)), col = "red", lwd = 3, add = T)

hist(modello6$residuals, breaks = 10, xlab="", ylab = "Frequenze relative",
     ylim = c(0,0.4), xlim = c(-12,12), freq = F, main="Residui di ARMA(2,2)", cex.main=0.5)
curve(dnorm(x, mean(modello6$residuals), sd(modello6$residuals)), col = "red", lwd = 3, add = T)

hist(modello7$residuals, breaks = 10, xlab="", ylab = "Frequenze relative",
     ylim = c(0,0.4), xlim = c(-12,12), freq = F, main="Residui di ARMA(3,0)", cex.main=0.5)
curve(dnorm(x, mean(modello7$residuals), sd(modello7$residuals)), col = "red", lwd = 3, add = T)

par(mfrow=c(4,2))
qqnorm(modello1$residuals, main="Normal Q-Q Plot ARMA(1,0)", cex.main = 0.5)
qqline(modello1$residuals)

qqnorm(modello2$residuals, main="Normal Q-Q Plot ARMA(0,1)", cex.main = 0.5)
qqline(modello2$residuals)

qqnorm(modello3$residuals, main="Normal Q-Q Plot ARMA(1,1)", cex.main = 0.5)
qqline(modello3$residuals)

qqnorm(modello4$residuals, main="Normal Q-Q Plot ARMA(2,1)", cex.main = 0.5)
qqline(modello4$residuals)

qqnorm(modello5$residuals, main="Normal Q-Q Plot ARMA(1,2)", cex.main = 0.5)
qqline(modello5$residuals)

qqnorm(modello6$residuals, main="Normal Q-Q Plot ARMA(2,2)", cex.main = 0.5)
qqline(modello6$residuals)

qqnorm(modello7$residuals, main="Normal Q-Q Plot ARMA(3,0)", cex.main = 0.5)
qqline(modello7$residuals)

shapiro1 <- shapiro.test(modello1$residuals)
print(paste("Test Shapiro-Wilk per ARMA(1,0): p-value =", shapiro1$p.value))
shapiro2 <- shapiro.test(modello2$residuals)
print(paste("Test Shapiro-Wilk per ARMA(0,1): p-value =", shapiro2$p.value))
shapiro3 <- shapiro.test(modello3$residuals)
print(paste("Test Shapiro-Wilk per ARMA(1,1): p-value =", shapiro3$p.value))
shapiro4 <- shapiro.test(modello4$residuals)
print(paste("Test Shapiro-Wilk per ARMA(2,1): p-value =", shapiro4$p.value))
shapiro5 <- shapiro.test(modello5$residuals)
print(paste("Test Shapiro-Wilk per ARMA(1,2): p-value =", shapiro5$p.value))
shapiro6 <- shapiro.test(modello6$residuals)
print(paste("Test Shapiro-Wilk per ARMA(2,2): p-value =", shapiro6$p.value))
shapiro7 <- shapiro.test(modello7$residuals)
print(paste("Test Shapiro-Wilk per ARMA(3,0): p-value =", shapiro7$p.value))

par(mfrow=c(4,2))
ts.plot(y_seasonal_diff + modello1$residuals, y_seasonal_diff, lty = c(1, 2), col = c("red", "blue"),
        main="Serie stimata con ARMA(1,0)")
legend("bottomright", legend = c("ss fit", "ss oss"), lty = c(1, 2), col = c("red", "blue"), cex = 0.5)

ts.plot(y_seasonal_diff + modello2$residuals, y_seasonal_diff, lty = c(1, 2), col = c("red", "blue"),
        main="Serie stimata con ARMA(0,1)")
legend("bottomright", legend = c("ss fit", "ss oss"), lty = c(1, 2), col = c("red", "blue"), cex = 0.5)

ts.plot(y_seasonal_diff + modello3$residuals, y_seasonal_diff, lty = c(1, 2), col = c("red", "blue"),
        main="Serie stimata con ARMA(1,1)")
legend("bottomright", legend = c("ss fit", "ss oss"), lty = c(1, 2), col = c("red", "blue"), cex = 0.5)

ts.plot(y_seasonal_diff + modello4$residuals, y_seasonal_diff, lty = c(1, 2), col = c("red", "blue"),
        main="Serie stimata con ARMA(2,1)")
legend("bottomright", legend = c("ss fit", "ss oss"), lty = c(1, 2), col = c("red", "blue"), cex = 0.5)

ts.plot(y_seasonal_diff + modello5$residuals, y_seasonal_diff, lty = c(1, 2), col = c("red", "blue"),
        main="Serie stimata con ARMA(1,2)")
legend("bottomright", legend = c("ss fit", "ss oss"), lty = c(1, 2), col = c("red", "blue"), cex = 0.5)

ts.plot(y_seasonal_diff + modello6$residuals, y_seasonal_diff, lty = c(1, 2), col = c("red", "blue"),
        main="Serie stimata con ARMA(2,2)")
legend("bottomright", legend = c("ss fit", "ss oss"), lty = c(1, 2), col = c("red", "blue"), cex = 0.5)

ts.plot(y_seasonal_diff + modello7$residuals, y_seasonal_diff, lty = c(1, 2), col = c("red", "blue"),
        main="Serie stimata con ARMA(3,0)")
legend("bottomright", legend = c("ss fit", "ss oss"), lty = c(1, 2), col = c("red", "blue"), cex = 0.5)


## STIMA DEI MODELLI SARIMA

sarima1 <- Arima(y, order=c(0,1,3), seasonal = list(order=c(0,4,1), period=4))
sarima2 <- Arima(y, order=c(1,1,2), seasonal = list(order=c(2,4,2), period=4))
sarima3 <- Arima(y, order=c(1,1,3), seasonal = list(order=c(0,4,2), period=4))
sarima4 <- Arima(y, order=c(2,1,2), seasonal = list(order=c(2,4,2), period=4))
sarima5 <- Arima(y, order=c(0,1,1), seasonal = list(order=c(1,4,1), period=4))


## RISULTATI 

df <- data.frame(
  Modello = c("1", "2", "3", "4", "5"),
  AIC = c(sarima1$aic, sarima2$aic, sarima3$aic, sarima4$aic, sarima5$aic),
  AICc = c(sarima1$aicc, sarima2$aicc, sarima3$aicc, sarima4$aicc, sarima5$aicc),
  BIC = c(sarima1$bic, sarima2$bic, sarima3$bic, sarima4$bic, sarima5$bic),
  Significativi = c("sì", "no", "no", "no", "sì"),
  stringsAsFactors = FALSE
)
etichette <- c("SARIMA(0,1,3)x(0,4,1)", "SARIMA(1,1,2)x(2,4,2)", "SARIMA(1,1,3)x(0,4,2)", "SARIMA(2,1,2)x(2,4,2)", "SARIMA(0,1,1)x(0,4,1)")
df$Modello <- etichette


par(mfrow=c(1,1))
par(mfrow=c(3,2))
plot(sarima1$residuals, type="l", ylab="Residui", main = "Residui di SARIMA(0,1,3)x(0,4,1)")
plot(sarima2$residuals, type="l", ylab="Residui", main = "Residui di SARIMA(1,1,2)x(2,4,2)")
plot(sarima3$residuals, type="l", ylab="Residui", main = "Residui di SARIMA(1,1,3)x(0,4,2)")
plot(sarima4$residuals, type="l", ylab="Residui", main = "Residui di SARIMA(2,1,2)x(2,4,2)")
plot(sarima4$residuals, type="l", ylab="Residui", main = "Residui di SARIMA(0,1,1)x(0,4,1)") 

Boxtest1<-Box.test(sarima1$residuals, lag=30, type="Box-Pierce", fitdf=4)
Boxtest2<-Box.test(sarima2$residuals, lag=30, type="Box-Pierce", fitdf=7)
Boxtest3<-Box.test(sarima3$residuals, lag=30, type="Box-Pierce", fitdf=6)
Boxtest4<-Box.test(sarima4$residuals, lag=30, type="Box-Pierce", fitdf=8)
Boxtest5<-Box.test(sarima5$residuals, lag=30, type="Box-Pierce", fitdf=2)

Ljungtest1<-Box.test(sarima1$residuals, lag=30, type = "Ljung-Box", fitdf=4)
Ljungtest2<-Box.test(sarima2$residuals, lag=30, type = "Ljung-Box", fitdf=7)
Ljungtest3<-Box.test(sarima3$residuals, lag=30, type = "Ljung-Box", fitdf=6)
Ljungtest4<-Box.test(sarima4$residuals, lag=30, type = "Ljung-Box", fitdf=8)
Ljungtest5<-Box.test(sarima5$residuals, lag=30, type = "Ljung-Box", fitdf=2)

BoxPierce <- c(Boxtest1$p.value, Boxtest2$p.value, Boxtest3$p.value, Boxtest4$p.value, Boxtest5$p.value)
LjungBox <- c(Ljungtest1$p.value, Ljungtest2$p.value, Ljungtest3$p.value, Ljungtest4$p.value, Ljungtest5$p.value)
Incorrelati <- c("no", "no", "sì", "test discordanti", "sì" )

results <- data.frame(Model = etichette, BoxPierce = BoxPierce, LjungBox = LjungBox, Incorrelati)

print(results)

print(results)
par(mfrow=c(1,1))
par(mfrow=c(3,2))

acf(sarima1$residuals, main = "ACF dei residui di SARIMA(0,1,3)x(0,4,1)")
acf(sarima2$residuals, main = "ACF dei residui di SARIMA(1,1,2)x(2,4,2)")
acf(sarima3$residuals, main = "ACF dei residui di SARIMA(1,1,3)x(0,4,2)")
acf(sarima4$residuals, main = "ACF dei residui di SARIMA(2,1,2)x(2,4,2)")
acf(sarima5$residuals, main = "ACF dei residui di SARIMA(0,1,1)x(0,4,1)")

par(mfrow=c(1,1))
par(mfrow=c(3,2))
hist(sarima1$residuals, breaks = 10, xlab="", ylab = "Frequenze relative",
     ylim = c(0,0.4), xlim = c(-12,12), freq = F, main="Residui di SARIMA(0,1,3)x(0,4,1)", cex.main=0.5)
curve(dnorm(x, mean(sarima1$residuals), sd(sarima1$residuals)), col = "red", lwd = 3, add = T)

hist(sarima2$residuals, breaks = 10, xlab="", ylab = "Frequenze relative",
     ylim = c(0,0.4), xlim = c(-12,12), freq = F, main="Residui di SARIMA(1,1,2)x(2,4,2)", cex.main=0.5)
curve(dnorm(x, mean(sarima2$residuals), sd(sarima2$residuals)), col = "red", lwd = 3, add = T)

hist(sarima3$residuals, breaks = 10, xlab="", ylab = "Frequenze relative",
     ylim = c(0,0.4), xlim = c(-12,12), freq = F, main="Residui di SARIMA(1,1,3)x(0,4,2)", cex.main=0.5)
curve(dnorm(x, mean(sarima3$residuals), sd(sarima3$residuals)), col = "red", lwd = 3, add = T)

hist(sarima4$residuals, breaks = 10, xlab="", ylab = "Frequenze relative",
     ylim = c(0,0.4), xlim = c(-12,12), freq = F, main="Residui di SARIMA(2,1,2)x(2,4,2)", cex.main=0.5)
curve(dnorm(x, mean(sarima4$residuals), sd(sarima4$residuals)), col = "red", lwd = 3, add = T)

hist(sarima5$residuals, breaks = 10, xlab="", ylab = "Frequenze relative",
     ylim = c(0,0.4), xlim = c(-12,12), freq = F, main="Residui di SARIMA(0,1,1)x(0,4,1)", cex.main=0.5)
curve(dnorm(x, mean(sarima5$residuals), sd(sarima5$residuals)), col = "red", lwd = 3, add = T)

par(mfrow=c(1,1))
par(mfrow=c(3,2))

qqnorm(sarima1$residuals, main="Normal Q-Q Plot SARIMA(0,1,3)x(0,4,1)", cex.main = 0.5)
qqline(sarima1$residuals)

qqnorm(sarima2$residuals, main="Normal Q-Q Plot SARIMA(1,1,2)x(2,4,2)", cex.main = 0.5)
qqline(sarima2$residuals)

qqnorm(sarima3$residuals, main="Normal Q-Q Plot SARIMA(1,1,3)x(0,4,2)", cex.main = 0.5)
qqline(sarima3$residuals)

qqnorm(sarima4$residuals, main="Normal Q-Q Plot SARIMA(2,1,2)x(2,4,2)", cex.main = 0.5)
qqline(sarima4$residuals)

qqnorm(sarima5$residuals, main="Normal Q-Q Plot SARIMA(0,1,1)x(0,4,1)", cex.main = 0.5)
qqline(sarima5$residuals)

shapiro1 <- shapiro.test(sarima1$residuals)
print(paste("Test Shapiro-Wilk per SARIMA(0,1,3)x(0,4,1): p-value =", shapiro1$p.value))
shapiro2 <- shapiro.test(sarima2$residuals)
print(paste("Test Shapiro-Wilk per SARIMA(1,1,2)x(2,4,2): p-value =", shapiro2$p.value))
shapiro3 <- shapiro.test(sarima3$residuals)
print(paste("Test Shapiro-Wilk per SARIMA(1,1,3)x(0,4,2): p-value =", shapiro3$p.value))
shapiro4 <- shapiro.test(sarima4$residuals)
print(paste("Test Shapiro-Wilk per SARIMA(2,1,2)x(2,4,2): p-value =", shapiro4$p.value))
shapiro5 <- shapiro.test(sarima5$residuals)
print(paste("Test Shapiro-Wilk per SARIMA(0,1,1)x(0,4,1): p-value =", shapiro5$p.value))

par(mfrow=c(1,1))
par(mfrow=c(3,2))
ts.plot(y + sarima1$residuals, y, lty = c(1, 2), col = c("red", "blue"),
        main="Serie stimata con SARIMA(0,1,3)x(0,4,1)")
legend("bottomright", legend = c("ss fit", "ss oss"), lty = c(1, 2), col = c("red", "blue"), cex = 0.5)

ts.plot(y + sarima2$residuals, y, lty = c(1, 2), col = c("red", "blue"),
        main="Serie stimata con SARIMA(1,1,2)x(2,4,2)")
legend("bottomright", legend = c("ss fit", "ss oss"), lty = c(1, 2), col = c("red", "blue"), cex = 0.5)

ts.plot(y + sarima3$residuals, y, lty = c(1, 2), col = c("red", "blue"),
        main="Serie stimata con SARIMA(1,1,3)x(0,4,2)")
legend("bottomright", legend = c("ss fit", "ss oss"), lty = c(1, 2), col = c("red", "blue"), cex = 0.5)

ts.plot(y + sarima4$residuals, y, lty = c(1, 2), col = c("red", "blue"),
        main="Serie stimata con SARIMA(2,1,2)x(2,4,2)")
legend("bottomright", legend = c("ss fit", "ss oss"), lty = c(1, 2), col = c("red", "blue"), cex = 0.5)

ts.plot(y + sarima5$residuals, y, lty = c(1, 2), col = c("red", "blue"),
        main="Serie stimata con SARIMA(0,1,1)x(0,4,1)")
legend("bottomright", legend = c("ss fit", "ss oss"), lty = c(1, 2), col = c("red", "blue"), cex = 0.5)

coeftest(sarima5, vcov(sarima5))


## PREVISIONI

previsione_sarima = forecast(sarima5, h = 8, level = c(0.95, 0.99))
previsione_sarima

plot(previsione_sarima, include = 50, PI = TRUE, shadebars = TRUE, shaded = TRUE, main="Previsione con SARIMA(0,1,1)x(0,4,1)")

exp_window<- function(ts_data, N, n, model) {
  results <- c()  
  
  x <- length(ts_data)
  
  indice_inizio_validation <- N + 1
  
  while (x > n) {
    train_set <- window(ts_data, start = 1, end = indice_inizio_validation-1)
    test_set <- window(ts_data, start = indice_inizio_validation, end = indice_inizio_validation + n - 1)
    
    model_fit <- Arima(train_set, model = model)
    
    forecasted <- forecast(model_fit, h = n)
    
    rmse_value <- rmse(test_set, forecasted$mean)
    
    results <- c(results, rmse_value)
    
    x <- x - indice_inizio_validation 
    
    indice_inizio_validation <- indice_inizio_validation + n
  }
  
  train_set <- window(ts_data, start = 1, end = length(ts_data)-n)
  test_set <- window(ts_data, start = length(ts_data)-n+1, end = length(ts_data))
  
  model_fit <- Arima(train_set, model = model)
  
  forecasted <- forecast(model_fit, h = n)
  
  rmse_value <- rmse(test_set, forecasted$mean)
  
  results <- c(results, rmse_value)
  
  mean_rmse <- mean(results)
  
  return(mean_rmse)
}

slid_windows<- function(ts_data, N, n, model) {
  results <- c()
  
  indice_inizio_train <- 1
  indice_inizio_validation <- N+1
  
  while (length(ts_data) - indice_inizio_train >= N+n) {
    train_set <- window(ts_data, start = indice_inizio_train, end = indice_inizio_train+N-1)
    test_set <- window(ts_data, start = indice_inizio_validation, end = indice_inizio_validation + n-1)
    
    model_fit <- Arima(train_set, model = model)
    
    forecasted <- forecast(model_fit, h = n)
    
    rmse_value <- rmse(test_set, forecasted$mean)
    
    results <- c(results, rmse_value)
    
    indice_inizio_train <- indice_inizio_train + n
    indice_inizio_validation <- indice_inizio_validation + n
  }
  
  train_set <- window(ts_data, start = length(ts_data)-N-n, end = length(ts_data)-n)
  test_set <- window(ts_data, start = length(ts_data)-n+1, end = length(ts_data))
  
  model_fit <- Arima(train_set, model = model)
  
  forecasted <- forecast(model_fit, h = n)
  
  rmse_value <- rmse(test_set, forecasted$mean)
  
  results <- c(results, rmse_value)
  
  mean_rmse <- mean(results)
  
  return(mean_rmse)
}

exp_window(y, 20, 8, sarima5)
slid_windows(y, 20, 8, sarima5)


## PREVISIONI CON RETI NEURALI 


fit_nnar <- nnetar(y)
fit_nnar$model
prev <- forecast(fit_nnar, h=8)
autoplot(prev, ylab="")

exp_window_rn<- function(ts_data, N, n) {
  results <- c()  
  
  x <- length(ts_data)
  
  indice_inizio_validation <- N + 1
  
  while (x > n) {
    train_set <- window(ts_data, start = 1, end = indice_inizio_validation-1)
    test_set <- window(ts_data, start = indice_inizio_validation, end = indice_inizio_validation + n - 1)
    
    fit_nnar <- nnetar(train_set)
    fit_nnar$model
    forecasted <- forecast(fit_nnar, h=n)
    
    rmse_value <- rmse(test_set, forecasted$mean)
    
    results <- c(results, rmse_value)
    
    x <- x - indice_inizio_validation 
    
    indice_inizio_validation <- indice_inizio_validation + n
  }
  
  train_set <- window(ts_data, start = 1, end = length(ts_data)-n)
  test_set <- window(ts_data, start = length(ts_data)-n+1, end = length(ts_data))
  
  fit_nnar <- nnetar(train_set)
  fit_nnar$model
  forecasted <- forecast(fit_nnar, h=n)
  
  rmse_value <- rmse(test_set, forecasted$mean)
  
  results <- c(results, rmse_value)
  
  mean_rmse <- mean(results)
  
  return(mean_rmse)
}

slid_windows_rn<- function(ts_data, N, n) {
  results <- c()
  
  indice_inizio_train <- 1
  indice_inizio_validation <- N+1
  
  while (length(ts_data) - indice_inizio_train >= N+n) {
    train_set <- window(ts_data, start = indice_inizio_train, end = indice_inizio_train+N-1)
    test_set <- window(ts_data, start = indice_inizio_validation, end = indice_inizio_validation + n-1)
    
    fit_nnar <- nnetar(train_set)
    fit_nnar$model
    forecasted <- forecast(fit_nnar, h=n)
    
    rmse_value <- rmse(test_set, forecasted$mean)
    
    results <- c(results, rmse_value)
    
    indice_inizio_train <- indice_inizio_train + n
    indice_inizio_validation <- indice_inizio_validation + n
  }
  
  train_set <- window(ts_data, start = length(ts_data)-N-n, end = length(ts_data)-n)
  test_set <- window(ts_data, start = length(ts_data)-n+1, end = length(ts_data))
  
  fit_nnar <- nnetar(train_set)
  fit_nnar$model
  forecasted <- forecast(fit_nnar, h=n)
  
  rmse_value <- rmse(test_set, forecasted$mean)
  
  results <- c(results, rmse_value)
  
  mean_rmse <- mean(results)
  
  return(mean_rmse)
}

exp_window_rn(y, 20, 8)
slid_windows_rn(y, 20, 8)

