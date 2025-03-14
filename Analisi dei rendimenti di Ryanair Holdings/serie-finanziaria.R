library(timeSeries)
library(fBasics)
library(tseries)
library(fracdiff)
library(FinTS)
library(fGarch)
library(TSA)
library(timeDate)
library(rugarch)
library(moments)
library(FinTS)
library(ggplot2)
library(lubridate)
library(Metrics)

## SETUP 

setwd(...)
dati = read.csv("RY4C.F.csv", header=TRUE)
y<-dati[,5]

## VISUALIZZAZIONE DEI DATI 

dati$Data <- as.Date(dati$Date, format="%Y-%m-%d")
dati$Anno <- year(dati$Data)
dati$Trimestre <- quarter(dati$Data)
grafico <- ggplot(dati, aes(x = Data, y = Close)) +
  geom_line() +  
  scale_x_date(date_labels = "%Y", date_breaks = "2 years") + 
  labs(x = "Anno", y = "Valore (EUR)", title = "Valore delle azioni di Ryanair Holdings dal 2010 al 2024") +
  theme_minimal()

print(grafico)

## ESTRAZIONE DEI RENDIMENTI

r = returns(y, method="continuous")
dati$rendimenti <- r
r <- r[-1]
r <- ts(r)
y<-ts(y)

## VISUALIZZAZIONE DEI RENDIMENTI 

grafico1 <- ggplot(dati[-1, ], aes(x = Data, y = rendimenti)) +
  geom_line() +  
  scale_x_date(date_labels = "%Y", date_breaks = "2 years") + 
  labs(x = "Anno", y = "Valore", title = "Rendimenti di Ryanair Holdings dal 2010 al 2024") +
  theme_minimal()

print(grafico1)
basicStats(r)

## VERIFICA DI NORMALITÀ SUI RENDIMENTI 

T<-3660
WNG<-rnorm(T, mean(r), sd(r))
WNG=ts(WNG)

ts.plot(r, WNG, main="Grafico dei rendimenti e\ndi un processo White Noise\na confronto"
        ,xlab = ".", ylab = ".", lty = c(1, 2), col = c("red", "blue"))
legend("bottomright",                   
       legend = c("Rendimenti", "White Noise"),  
       col = c("red", "blue"),  
       lty = c(1, 2),
       cex = 0.5)     
par(mfrow=c(1,2))
hist(r, breaks = 20, col = "lightblue", xlab = "Rendimenti", 
     ylab = "Frequenza", main = "Istogramma dei rendimenti \ne distribuzione normale", freq = FALSE,
     cex.main = 0.6)
curve(dnorm(x, mean = mean(r), sd = sd(r)), add = TRUE, col = "red", lwd = 2)
qqnormPlot(r, main="Normal QQplot", cex.main = 0.8, title=FALSE)
par(mfrow=c(1,1))

ksnormTest(r, title="Test di Kolmogorov-Smirnov", 
           description="Ipotesi nulla: la distribuzione dei rendimenti è normale")
shapiroTest(r, title = "Tests di Shapiro-Wilk",
            description="Ipotesi nulla: la distribuzione dei rendimenti è normale")


## CORRELAZIONE

acf(r, lag.max=60, main="ACF dei rendimenti", cex.main=0.8)
Box.test(r, lag=6, type="Ljung-Box")
Box.test(r, lag=6, type="Box-Pierce")
lag.plot(r, lags=6, do.lines = F, main="Lag plot dei rendimenti per sfasamenti da 1 a 6", cex.main=0.8)

acf(r^2, lag.max=60, main="ACF dei rendimenti \nal quadrato", cex.main=0.8)
Box.test(r^2, lag=6, type="Ljung-Box")
Box.test(r^2, lag=6, type="Box-Pierce")
lag.plot(r^2, lags=6, do.lines = F, main="Lag plot dei rendimenti al quadrato", cex.main=0.8)

ArchTest(r)

acf(abs(r), lag.max=60, main="ACF dei rendimenti\nin valore assoluto", cex.main=0.8)
Box.test(abs(r), lag=6, type="Ljung-Box")
Box.test(abs(r), lag=6, type="Box-Pierce")

## STIMA DEI MODELLI 

spec_garch <- ugarchspec(variance.model = list(model = "sGARCH"), mean.model = list(armaOrder = c(0, 0)))
fit_garch <- ugarchfit(spec = spec_garch, data = r)
fit_garch

spec_garch2 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2,2)), mean.model = list(armaOrder = c(0, 0)))
fit_garch2 <- ugarchfit(spec = spec_garch2, data = r)
fit_garch2

spec_tgarch <- ugarchspec(variance.model = list(model = "fGARCH", submodel = "TGARCH"), mean.model = list(armaOrder = c(0, 0)))
fit_tgarch <- ugarchfit(spec = spec_tgarch, data = r)
fit_garch

spec_egarch <- ugarchspec(variance.model = list(model = "eGARCH"), mean.model = list(armaOrder = c(0, 0)))
fit_egarch <- ugarchfit(spec = spec_egarch, data = r)
fit_egarch

spec_egarch2 <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(2,2)), mean.model = list(armaOrder = c(0, 0)))
fit_egarch2 <- ugarchfit(spec = spec_egarch2, data = r)
fit_egarch2

ic_garch <- infocriteria(fit_garch)
ic_garch2 <- infocriteria(fit_garch2)
ic_tgarch <- infocriteria(fit_tgarch)
ic_egarch <- infocriteria(fit_egarch)
ic_egarch2 <- infocriteria(fit_egarch2)

aic_garch <- ic_garch["Akaike",1]
bic_garch <- ic_garch["Bayes",1]

aic_garch2 <- ic_garch2["Akaike",1]
bic_garch2 <- ic_garch2["Bayes",1]

aic_tgarch <- ic_tgarch["Akaike",1]
bic_tgarch <- ic_tgarch["Bayes",1]

aic_egarch <- ic_egarch["Akaike",1]
bic_egarch <- ic_egarch["Bayes",1]

aic_egarch2 <- ic_egarch2["Akaike",1]
bic_egarch2 <- ic_egarch2["Bayes",1]

aic_values <- c(GARCH = aic_garch, GARCH2 = aic_garch2, TGARCH = aic_tgarch, EGARCH = aic_egarch, EGARCH2 = aic_egarch)
bic_values <- c(GARCH = bic_garch, GARCH2 = bic_garch2, TGARCH = bic_tgarch, EGARCH = bic_egarch, EGARCH2 = aic_egarch2)

aic_values
bic_values

nic <- newsimpact(fit_egarch2)
plot(nic$zx, nic$zy, type = "l", xlab = "", ylab = "")

## PREVISIONI

previsioni<-ugarchforecast(fit_egarch2, n.ahead = 10)
plot(previsioni, which = 3, main="Previsioni tramite eGARCH(2,2)")

exp_window<- function(ts_data, N, n, spec) {
  results <- c()  # Vettore per salvare gli RMSE
  
  x <- length(ts_data)
  
  indice_inizio_validation <- N + 1
  
  while (x > n) {
    train_set <- window(ts_data, start = 1, end = indice_inizio_validation-1)
    test_set <- window(ts_data, start = indice_inizio_validation, end = indice_inizio_validation + n - 1)
    
    model_fit <- ugarchfit(spec = spec, data = train_set)
    
    forecasted <- ugarchforecast(model_fit, n.ahead = n)
    
    rmse_value <- rmse(test_set, as.numeric(forecasted@forecast$sigmaFor))
    
    results <- c(results, rmse_value)
    
    x <- x - indice_inizio_validation 
    
    indice_inizio_validation <- indice_inizio_validation + n
  }
  
  train_set <- window(ts_data, start = 1, end = length(ts_data)-n)
  test_set <- window(ts_data, start = length(ts_data)-n+1, end = length(ts_data))
  
  model_fit <- ugarchfit(spec = spec, data = ts_data)
  
  forecasted <- ugarchforecast(model_fit, n.ahead = n)
  
  rmse_value <- rmse(test_set, as.numeric(forecasted@forecast$sigmaFor))
  
  results <- c(results, rmse_value)
  
  mean_rmse <- mean(results)
  
  return(mean_rmse)
}

exp_window(r, 500, 100, spec_egarch2)

