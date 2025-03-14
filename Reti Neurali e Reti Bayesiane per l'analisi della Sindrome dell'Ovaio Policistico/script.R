setwd("C:/Users/franc/OneDrive/Desktop/Neural Networks")
library(openxlsx)
library(RefManageR)
library(bibtex)
library(ggplot2)
library(vcd)
library(binom)
library(PropCIs)
library(car)
library(dplyr)
library(emmeans)
library(knitr)
library(caret)
library(MASS)
library(reshape2)
library(nnet)
library(neuralnet)
library(gridExtra)
library(NeuralNetTools)
library(RSNNS)
library(lattice)
library(BiocManager)
library(gRain)
library(Rgraphviz)
library(bnlearn)

# Caricamento del dataset
pcos <- read.xlsx("PCOS_data_without_infertility.xlsx", sheet = 2)

# Rinominazione delle colonne
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

# Trasformazioni delle variabili 
pcos$FSH_LH <- 1 / pcos$FSH_LH
pcos <- pcos %>% rename(LHFSH = `FSH_LH`)
pcos$Cycle[pcos$Cycle == 5] <- 1

# Gestione dei valori mancanti
pcos$AMH <- as.numeric(pcos$AMH)
pcos$AMH[is.na(pcos$AMH)] <- median(pcos$AMH, na.rm = TRUE)
pcos$LHFSH[is.infinite(pcos$LHFSH)] <- median(pcos$LHFSH, na.rm = TRUE)
pcos$II.beta.HCG <- as.numeric(pcos$II.beta.HCG)
pcos$II.beta.HCG[is.infinite(pcos$II.beta.HCG)] <- median(pcos$II.beta.HCG, na.rm = TRUE)
pcos$Pulse.rate[pcos$Pulse.rate <= 20] <- median(pcos$Pulse.rate, na.rm = TRUE)
pcos$Fast.food[is.na(pcos$Fast.food)] <- 
  as.numeric(names(sort(table(pcos$Fast.food), decreasing = TRUE))[1])

cols_with_na <- c("Marriage.status", "II.beta.HCG")
for (col in cols_with_na) {
  pcos[[col]][is.na(pcos[[col]])] <- median(pcos[[col]], na.rm = TRUE)
}

sum(is.na(pcos))
sum(is.infinite(as.matrix(pcos)))

# Creazione di nuove variabili e gestione dei valori mancanti
pcos$Follicles <- pmax(pcos$FolliclesL, pcos$FolliclesR)
pcos$F.size <- pmin(pcos$L.size, pcos$R.size)

# Rimozione di colonne non utili (id delle pazienti)
pcos <- pcos[, 3:46]

# Costruzione della variabile target PCOStype
set.seed(77)
pcos$PCOStype <- NA
pcos$PCOStype[pcos$PCOS == 0] <- "None"
pcos$PCOStype[pcos$PCOS == 1 & pcos$AMH >= 5] <- "Anovulatory"
pcos$PCOStype[pcos$PCOS == 1 & pcos$AMH < 5] <- "Ovulatory"

# Normalizzazione dei dati 
max_vals <- as.numeric(apply(pcos[, 1:44], 2, max))
min_vals <- as.numeric(apply(pcos[, 1:44], 2, min))
scaled_data <- round(scale(pcos[, 1:44], center = min_vals, scale = max_vals - min_vals), 4)
pcos_norm <- data.frame(scaled_data, PCOStype = pcos$PCOStype)

# Rimozione colonna ridondante (PCOS da cui deriva PCOStype)
pcos_norm <- pcos_norm[, -1]
pcos_norm$PCOStype <- as.factor(pcos_norm$PCOStype)

# Matrice di correlazione 
cor.matrix <- cor(pcos_norm[, -c(44)])
correlation_matrix_melted <- melt(cor.matrix)
ggplot(correlation_matrix_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() + 
  scale_fill_gradient2(midpoint = 0, low = "blue", high = "red", mid = "white") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Ruota le etichette dell'asse x
  labs(title = "Matrice di Correlazione", x = "", y = "")

# Eliminazione delle variabili troppo correlate con altre
pcos_norm <- pcos_norm[, -c(2, 3, 11, 16, 17, 19, 20, 37, 38, 39, 40)]

# Suddivisione in training set e test set
set.seed(123)
train_idx <- createDataPartition(pcos_norm$PCOStype, p = 0.85, list = FALSE)
train_data <- pcos_norm[train_idx, ]
test_data <- pcos_norm[-train_idx, ]
train_data$PCOStype <- factor(train_data$PCOStype, levels = c("None", "Anovulatory", "Ovulatory"))
test_data$PCOStype <- factor(test_data$PCOStype, levels = c("None", "Anovulatory", "Ovulatory"))

# Distribuzione delle classi di PCOStype nel training set
train_dist <- table(train_data$PCOStype)
train_percentage <- prop.table(train_dist) * 100

# Distribuzione delle classi di PCOStype nel test set
test_dist <- table(test_data$PCOStype)
test_percentage <- prop.table(test_dist) * 100

# Creazione di una nuova colonna per indicare se il dato è nel training o nel test set
# (necessario per il plot)
train_data$Set <- "Training"
test_data$Set <- "Test"
combined_data <- rbind(train_data, test_data)

# Proporzioni di ciascuna classe di PCOStype nei due set
combined_data_summary <- combined_data %>%
  group_by(PCOStype, Set) %>%
  summarise(Count = n()) %>%
  group_by(Set) %>%
  mutate(Proportion = Count / sum(Count)) # Calcoliamo la proporzione

# Grafico a barre con le proporzioni
ggplot(combined_data_summary, aes(x = PCOStype, y = Proportion, fill = Set)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  facet_wrap(~ PCOStype, ncol = 3) +  # Definiamo 3 colonne per i grafici
  ggtitle("Proporzione delle Classi per Training e Test Set") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Eliminazione di colonne superflue dai set
train_data <- train_data[, -c(34)]
test_data <- test_data[, -c(34)]

# Definizione dell'insieme di parametri che si vuole testare
thresholds <- c(0.5, 0.1, 0.001)
neurons_in_hidden_layer <- c(8, 16, 5)
decay_val <- c(1, 0.5, 0.1, 0.001)

# Inizializzazione di alcuni parametri
best_f1_score <- 0
best_model <- NULL
best_params <- list()

# Funzione per calcolare la metrica F1-score pesata 
f1_score <- function(predicted, expected, positive.class="1") {
  predicted <- factor(as.character(predicted), levels=unique(as.character(expected)))
  expected  <- as.factor(expected)
  cm = as.matrix(table(expected, predicted))
  
  
  precision <- diag(cm) / colSums(cm)
  recall <- diag(cm) / rowSums(cm)
  
  
  f1 <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
  
  
  f1[is.na(f1)] <- 0
  
  
  class_counts <- table(expected)
  
  
  weighted_f1 <- sum(f1 * class_counts) / sum(class_counts)
  
  
  return(weighted_f1)
}

# Ciclo per testare tutte le combinazioni di parametri
  for (thresh in thresholds) {
    for (neurons in neurons_in_hidden_layer) {
      for (decay_value in decay_val){
      
      set.seed(15)
      nn_model <- nnet(PCOStype ~ ., data = train_data, 
                       size = neurons,
                       maxit = 2000,
                       decay = decay_value,
                       type = "class")
      
      
      predictions <- predict(nn_model, test_data, type = "class")
      predictions <- factor(predictions, levels = c("None", "Anovulatory", "Ovulatory"))
      
      current_f1_score <- f1_score(predictions, test_data$PCOStype)
      
      # Aggiornamento del miglior modello se l'F1 score è migliore
      if (current_f1_score > best_f1_score) {
        best_f1_score <- current_f1_score
        best_model <- nn_model
        best_params <- list(threshold = thresh, neurons = neurons, decay = decay_value)
      }
    }
  }
}

# Previsione finale del miglior modello sui dati di test
test_predictions <- predict(best_model, test_data, type = "class")
test_predictions <- factor(test_predictions, levels = c("None", "Anovulatory", "Ovulatory"))
test_data$PCOStype <- factor(test_data$PCOStype, levels = c("None", "Anovulatory", "Ovulatory"))
# Calcolo dell'F1 score sul test set
final_f1_score <- f1_score(test_predictions, test_data$PCOStype)
# Calcolo dell'accuratezza sul test set
accuracy <- sum(test_predictions == test_data$PCOStype) / length(test_data$PCOStype)

plotnet(best_model)
conf_matrix <- table(Predicted = test_predictions, Actual = test_data$PCOStype)
conf_matrix_df <- as.data.frame(as.table(conf_matrix))
ggplot(data = conf_matrix_df, aes(x = Predicted, y = Actual, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  geom_text(aes(label = Freq), color = "white", size = 5) +
  theme_minimal() +
  labs(title = "Matrice di confusione", x = "Classe predetta", y = "Classe reale") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Preparazione dei dati per la libreria RSNNS
X_train <- train_data[, -c(33)]
X_test <- test_data[, -c(33)]

Y_train <- data.frame(
  None = ifelse(train_data$PCOStype == "None", 1, 0),
  Anovulatory = ifelse(train_data$PCOStype == "Anovulatory", 1, 0),
  Ovulatory = ifelse(train_data$PCOStype == "Ovulatory", 1, 0)
)

# Sviluppo della rete neurale con due strati nascosti
set.seed(123)
net2 <- mlp(X_train, Y_train, 
             size = c(10, 5),  
             maxit = 10000, 
             learnFuncParams = c(0.1, 0.01),
             linOut = FALSE)  

# Previsioni sul test set in termini di probabilità
pred_probs <- predict(net2, X_test)
colnames(pred_probs) <- c("None", "Anovulatory", "Ovulatory") 
# Conversioni delle probabilità in classi
pred_labels <- colnames(pred_probs)[max.col(pred_probs, ties.method = "first")]
# Conversione in fattore con gli stessi livelli di y_test
pred_labels <- factor(pred_labels, levels = c("None", "Anovulatory", "Ovulatory"))
true_labels <- test_data$PCOStype
# Matrice di confusione
conf_matrix <- table(pred_labels, true_labels)
# Calcolo della F1score
f1_score2 <- f1_score(pred_labels, true_labels)
# Calcolo dell'accuratezza
accuracy2 <- sum(pred_labels == test_data$PCOStype) / length(test_data$PCOStype)
conf_matrix2 <- table(Predicted = pred_labels, Actual = true_labels)
conf_matrix2_df <- as.data.frame(as.table(conf_matrix2))
ggplot(data = conf_matrix2_df, aes(x = Predicted, y = Actual, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  geom_text(aes(label = Freq), color = "white", size = 5) +
  theme_minimal() +
  labs(title = "Matrice di confusione per la rete neurale con 2 strati nascosti", x = "Classe predetta", y = "Classe reale") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Sviluppo della rete neurale con tre strati nascosti
set.seed(72)
net3 <- mlp(X_train, Y_train, 
             size = c(20, 10, 5),  
             maxit = 10000, 
              learnFuncParams = c(0.1, 0.01),
             linOut = FALSE)  # Softmax implicita per classificazione

# Previsioni sul test set in termini di probabilità
pred_probs3 <- predict(net3, X_test)
colnames(pred_probs3) <- c("None", "Anovulatory", "Ovulatory")  # Imposta i nomi delle classi
# Conversioni delle probabilità in classi
pred_labels3 <- colnames(pred_probs3)[max.col(pred_probs3, ties.method = "first")]
pred_labels3 <- factor(pred_labels3, levels = c("None", "Anovulatory", "Ovulatory"))
true_labels3 <- test_data$PCOStype
# Matrice di confusione
conf_matrix3 <- table(pred_labels3, true_labels3)
# Calcolo della F1-score
f1_score3 <- f1_score(pred_labels3, true_labels3)
# Accuratezza
accuracy3 <- sum(pred_labels3 == test_data$PCOStype) / length(test_data$PCOStype)
conf_matrix3 <- table(Predicted = pred_labels3, Actual = true_labels3)
conf_matrix3_df <- as.data.frame(as.table(conf_matrix3))
ggplot(data = conf_matrix3_df, aes(x = Predicted, y = Actual, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  geom_text(aes(label = Freq), color = "white", size = 5) +
  theme_minimal() +
  labs(title = "Matrice di confusione per la rete neurale con tre strati nascosti", x = "Classe predetta", y = "Classe reale") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Modello di regressione multinomiale
train_data$PCOStype <- as.factor(train_data$PCOStype)
set.seed(70)
model <- multinom(PCOStype ~., data = train_data, trace = FALSE)
model
# Previsioni con il modello di regressione multinomiale
preds <- predict(model, newdata = test_data)
# Accuratezza
accuracy_multinom <- mean(preds == test_data$PCOStype)
# F1-score
f1_score_multinom <- f1_score(preds, test_data$PCOStype)
# Matrice di confusioni
conf_matrix_multinom <- table(Predicted = preds, Actual = test_data$PCOStype)
conf_matrix_multinom_df <- as.data.frame(as.table(conf_matrix_multinom))
ggplot(data = conf_matrix_multinom_df, aes(x = Predicted, y = Actual, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  geom_text(aes(label = Freq), color = "white", size = 5) +
  theme_minimal() +
  labs(title = "Matrice di confusione per il modello di regressione multinomiale", x = "Classe predetta", y = "Classe reale") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

results_final <- data.frame(
  Model = c("Rete neurale 1 strato", "Rete neurale 2 strati", "Rete neurale 3 strati", "Regressione multinomiale"),
  F1score = c(round(best_f1_score, 4), round(f1_score2, 4), round(f1_score3, 4), round(f1_score_multinom, 4)),
  Accuracy = c(round(accuracy, 4), round(accuracy2, 4), round(accuracy3, 4), round(accuracy_multinom, 4))
)

print(results_final)

par(mfrow = c(1,2))
# Definizione dei colori 
blue_colors <- c("#B0C4DE", "#4682B4", "#1E90FF", "#0D1B2A")  

# Bar plot delle F1score
bp <- barplot(
  height = results_final$F1score,            
  #names.arg = results_final$Model,           
  col = blue_colors[1:nrow(results_final)],  
  main = "F1-score per modello", 
  ylab = "F1-score", 
  ylim = c(0,1),
  cex.names = 0.8,
  xaxt = "n" 
)
abline(h=0)
# Aggiunta dei nomi
text(x = bp, y = par("usr")[3] - 0.02, labels = results_final$Model, 
     srt = 45, adj = 1, xpd = TRUE, cex = 0.8)

bp1 <- barplot(
  height = results_final$Accuracy,                    
  col = blue_colors[1:nrow(results_final)],  
  main = "Accuracy per modello", 
  ylab = "Accuracy", 
  ylim = c(0,1),
  cex.names = 0.8,
  xaxt = "n" 
)
abline(h=0)
# Aggiunta dei nomi
text(x = bp1, y = par("usr")[3] - 0.02, labels = results_final$Model, 
     srt = 45, adj = 1, xpd = TRUE, cex = 0.8)

###############################################################################

# Reti bayesiane 

# Preparazione del dataset 
pcos_bayesian <- pcos[, -c(1, 3, 4, 12, 17, 18, 20, 21, 38, 39, 40, 41)]
pcos_bayesian <- pcos_bayesian[, -c(3, 4, 5, 6, 8, 9, 10, 11, 12, 18, 21, 23, 26, 27, 28, 29)]
pcos_bayesian$BMIlevel <- cut(pcos$BMI, breaks = c(-Inf, 18.5, 25, 30, Inf), 
                              labels = c("Underweight", "Normal Weight", "Overweight", "Obese"), right = FALSE)
pcos_bayesian <- pcos_bayesian[, -c(2)]
pcos_bayesian$HighLHFSH <- ifelse(pcos_bayesian$LHFSH >= 2,  1, 0) 
pcos_bayesian <- pcos_bayesian[, -c(3)]
pcos_bayesian <- pcos_bayesian[, -c(3)]
pcos_bayesian$TSH_Category <- cut(pcos_bayesian$TSH, 
                         breaks = c(-Inf, 0.39, 4.5, Inf), 
                         labels = c("Ipertiroidismo", "Normale", "Ipotiroidismo"))
pcos_bayesian <- pcos_bayesian[, -c(3)]
pcos_bayesian$AMH_Category <- cut(pcos_bayesian$AMH, 
                         breaks = c(-Inf, 1.0, 3.5, Inf), 
                         labels = c("Bassa riserva ovarica", "Normale", "Alta riserva ovarica"))
pcos_bayesian <- pcos_bayesian[, -c(3)]
pcos_bayesian <- pcos_bayesian[, -c(3, 4)]
pcos_bayesian$Hyperglycemia <- ifelse(pcos_bayesian$RBS >= 140, 1, 0)
pcos_bayesian <- pcos_bayesian[, -c(3)]
pcos_bayesian <- pcos_bayesian[, -c(4, 5)]
pcos_bayesian$Abnormal_Endometrium <- ifelse(pcos_bayesian$Endometrium <= 5 
                                           | pcos_bayesian$Endometrium >= 14, 1, 0)
pcos_bayesian <- pcos_bayesian[, -c(4)]
pcos_bayesian$High.Follicles <- ifelse(pcos_bayesian$Follicles >= 15, 1, 0)
pcos_bayesian <- pcos_bayesian[, -c(4)]
pcos_bayesian$Immature.follicles <- ifelse(pcos_bayesian$F.size < 10, 1, 0)
pcos_bayesian <- pcos_bayesian[, -c(4)]
pcos_bayesian$Age <- ifelse(pcos_bayesian$Age <= 29, "Giovane", ifelse(pcos_bayesian$Age >= 30 
                                                                        & pcos_bayesian$Age <= 39, 
                                            "Adulta", "Matura"))

pcos_bayesian <- as.data.frame(lapply(pcos_bayesian, as.factor))


# Definizione della prima blacklist
blacklist <- data.frame(
  from = rep("PCOStype", length(names(pcos_bayesian))),
  to = names(pcos_bayesian)
)
blacklist <- rbind(
  blacklist,
  data.frame(from = names(pcos_bayesian), to = "Age"),
  data.frame(from = names(pcos_bayesian), to = "BMIlevel"),
  data.frame(from = "Hair.growth", to = names(pcos_bayesian[,-c(4)]))
)

# Definizione della whitelist 
whitelist <- data.frame(
  from = c("Age", "BMIlevel"),
  to = c("PCOStype", "PCOStype")
)
# Sviluppo del primo dag
dag <- hc(pcos_bayesian, blacklist = blacklist, whitelist = whitelist, score = "bic")
score1 <- score(dag, pcos_bayesian)
graphviz.plot(dag, shape = "ellipse", fontsize=50, main="Grafo iniziale")

# Calcolo delle misure di strength nel primo dag
strengths <- arc.strength(dag, data = pcos_bayesian, criterion = "mi")
weak_arcs <- strengths[strengths$strength < 0.2, ]
print(weak_arcs)

# Test di indipendenza condizionale degli archi deboli
test_results <- list()
for (i in 1:nrow(weak_arcs)) {
  from_var <- weak_arcs$from[i]
  to_var <- weak_arcs$to[i]
  cond_vars <- setdiff(names(pcos_bayesian), c(from_var, to_var))
  test_result <- ci.test(from_var, to_var, cond_vars, data = pcos_bayesian)
  test_results[[paste(from_var, "→", to_var)]] <- test_result$p.value
}
print(test_results)

# Costruzione della blacklist aggiornata
blacklist <- data.frame(
  from = rep("PCOStype", length(names(pcos_bayesian))),
  to = names(pcos_bayesian)
)
blacklist <- rbind(
  blacklist,
  data.frame(from = names(pcos_bayesian), to = "Age"),
  data.frame(from = names(pcos_bayesian), to = "BMIlevel"),
  data.frame(from = "Hair.growth", to = names(pcos_bayesian[,-c(4)])),
  data.frame(
    from = c("TSH_Category", "Cycle", "Cycle", "Cycle", "Cycle", "High.Follicles"),
    to = c("Hair.growth", "HighLHFSH", "AMH_Category", "Hair.growth", "High.Follicles", "Abnormal_Endometrium"))
)

# Sviluppo del dag aggiornato
dag <- hc(pcos_bayesian, blacklist = blacklist, whitelist = whitelist,  score = "bic")
score2 <- score(dag, pcos_bayesian)

# Ultimo aggiornamento della whitelist
whitelist <- rbind(whitelist, data.frame(from = c("Hyperglycemia", "AMH_Category"),
                                         to = c("HighLHFSH", "Immature.follicles")))

# Sviluppo del dag definitivo
dag <- hc(pcos_bayesian, blacklist = blacklist, whitelist = whitelist,  score = "bic")
graphviz.plot(dag, shape = "ellipse", fontsize=50)
score3 <- score(dag, pcos_bayesian)

# Formula del modello
model_formula <- modelstring(dag)
print(model_formula)

# Stima delle distribuzioni di probabilità tramite il metodo di Bayes
bn <- bn.fit(dag, pcos_bayesian,  method = "bayes", iss = 10)
junction_tree = compile(as.grain(bn))

# BMI e PCOStype
querygrain(junction_tree, nodes = c("PCOStype", "BMIlevel"), type = "conditional") 
querygrain(junction_tree, nodes = c("PCOStype"), type = "marginal") 

# BMI e PCOStype anovulatoria
querygrain(setEvidence(junction_tree, nodes = "PCOStype", states = "Anovulatory"),
           nodes = "BMIlevel", type = "conditional")
querygrain(junction_tree, nodes = c("PCOStype"), type = "marginal") 

# BMI e PCOStype ovulatoria
querygrain(setEvidence(junction_tree, nodes = "PCOStype", states = "Ovulatory"),
           nodes = "BMIlevel", type = "conditional")
querygrain(junction_tree, nodes = c("PCOStype"), type = "marginal") 

# Età e PCOStype
querygrain(junction_tree, nodes = c("PCOStype", "Age"), type = "conditional") 
querygrain(junction_tree, nodes = c("PCOStype"), type = "marginal") 

# Iperglicemia e LHFSH elevato 
querygrain(setEvidence(junction_tree, nodes = "Hyperglycemia", states = "1"), 
           nodes = "HighLHFSH", type = "conditional")
querygrain(junction_tree, nodes = c("HighLHFSH"), type = "marginal")