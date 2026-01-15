############################################################
# TP3 - Régression logistique bayésienne + calibration
# Version "rendu TP" : figures ggplot2 + tableau final (gt)
############################################################

# --- Packages
suppressPackageStartupMessages({
  library(mlbench)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(brms)
  library(pROC)
  library(gt)
})

set.seed(1)

# ==========================================================
# 1) Données (PimaIndiansDiabetes2) + split train/test
# ==========================================================
data(PimaIndiansDiabetes2)
df <- na.omit(PimaIndiansDiabetes2)

# Variable cible en 0/1 (diabetes: "pos"/"neg")
df <- df %>%
  mutate(y = ifelse(diabetes == "pos", 1, 0))

# Split 70/30
n <- nrow(df)
idx <- sample.int(n, size = floor(0.7 * n))
train <- df[idx, ]
test  <- df[-idx, ]

# Modèle du TP3 : glucose + mass + age
form <- y ~ glucose + mass + age

# ==========================================================
# 2) Fit fréquentiste (glm)
# ==========================================================
m_glm <- glm(form, data = train, family = binomial(link = "logit"))
p_glm <- predict(m_glm, newdata = test, type = "response")

# ==========================================================
# 3) Fit bayésien (brms) : Bernoulli logit + prior N(0,2.5)
# ==========================================================
priors <- c(
  prior(normal(0, 2.5), class = "b"),
  prior(normal(0, 5),   class = "Intercept")
)

m_bayes <- stan_glm(
  y ~ glucose + mass + age,
  data = train,
  family = binomial(link = "logit"),
  prior = normal(0, 2.5),
  prior_intercept = normal(0, 5),
  chains = 2, iter = 1200, warmup = 600, refresh = 0
)

# Probas bayésiennes : E[p_i | data] via posterior_epred
p_draws <- posterior_linpred(m_bayes, newdata = test, transform = TRUE)
#p_draws <- posterior_epred(m_bayes, newdata = test, transform = TRUE)
p_bayes <- colMeans(p_draws)

y_test <- test$y

# ==========================================================
#  Fonctions métriques (NLL, Brier, ECE) + bins calibration
# ==========================================================
eps <- 1e-15

nll <- function(p, y){
  p <- pmin(pmax(p, eps), 1 - eps)
  -mean(y * log(p) + (1 - y) * log(1 - p))
}

brier <- function(p, y) mean((p - y)^2)

ece <- function(p, y, nb = 10){
  rel <- reliability_bins(p, y, nb = nb) %>% filter(n > 0)
  sum((rel$n / sum(rel$n)) * abs(rel$acc - rel$conf))
}

# ==========================================================
# (Bayes vs glm)
# ==========================================================

# Figure 2 : Distribution des probabilités prédites
library(tibble)
library(dplyr)
library(ggplot2)
pred_df <- tibble(
  p = c(p_bayes, p_glm),
  model = rep(c("Bayes (brms)", "Fréquentiste (glm)"), each = length(p_glm))
)

fig_hist <- ggplot(pred_df, aes(x = p)) +
  geom_histogram(bins = 25) +
  facet_wrap(~model, ncol = 1) +
  scale_x_continuous(limits = c(0,1)) +
  labs(
    title = "Distribution des probabilités prédites (test)",
    x = "Probabilité prédite P(y=1|x)",
    y = "Effectif"
  ) +
  theme_minimal(base_size = 12)

print(fig_hist)

# Optionnel : ROC
library(pROC)
roc_bayes <- roc(y_test, p_bayes, quiet = TRUE)
roc_glm   <- roc(y_test, p_glm, quiet = TRUE)

roc_df <- bind_rows(
  tibble(
    fpr = 1 - roc_bayes$specificities,
    tpr = roc_bayes$sensitivities,
    model = "Bayes (brms)"
  ),
  tibble(
    fpr = 1 - roc_glm$specificities,
    tpr = roc_glm$sensitivities,
    model = "Fréquentiste (glm)"
  )
)

fig_roc <- ggplot(roc_df, aes(x = fpr, y = tpr, linetype = model)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_path(linewidth = 0.9) +
  scale_x_continuous(limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  labs(
    title = "Courbes ROC (test)",
    x = "False Positive Rate (1 - spécificité)",
    y = "True Positive Rate (sensibilité)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# Affichage des figures

print(fig_hist)
print(fig_roc)

# ==========================================================
#  Métriques finales + tableau formaté
# ==========================================================
results <- tibble(
  model = c("Bayes (brms)", "Fréquentiste (glm)"),
  NLL   = c(nll(p_bayes, y_test), nll(p_glm, y_test)),
  Brier = c(brier(p_bayes, y_test), brier(p_glm, y_test)),
  AUC   = c(as.numeric(auc(roc_bayes)), as.numeric(auc(roc_glm))),
  ECE10 = c(ece(p_bayes, y_test, nb = 10), ece(p_glm, y_test, nb = 10))
)


# Tableau gt : arrondis + mise en forme
# et calcul des intervalles 

#  incertitude bayésienne sur NLL/Brier (IC 95% via tirages)
NLL_draws   <- apply(p_draws, 1, nll, y = y_test)
Brier_draws <- apply(p_draws, 1, brier, y = y_test)

bayes_ci <- tibble(
  model = "Bayes (brms)",
  NLL_q025 = as.numeric(quantile(NLL_draws, 0.025)),
  NLL_q975 = as.numeric(quantile(NLL_draws, 0.975)),
  Brier_q025 = as.numeric(quantile(Brier_draws, 0.025)),
  Brier_q975 = as.numeric(quantile(Brier_draws, 0.975))
)


library(gt)
tab <- results %>%
  left_join(bayes_ci, by = "model") %>%
  mutate(
    NLL_CI = ifelse(model == "Bayes (brms)",
                    sprintf("[%.3f, %.3f]", NLL_q025, NLL_q975), ""),
    Brier_CI = ifelse(model == "Bayes (brms)",
                      sprintf("[%.3f, %.3f]", Brier_q025, Brier_q975), "")
  ) %>%
  select(model, NLL, NLL_CI, Brier, Brier_CI, AUC, ECE10) %>%
  gt() %>%
  fmt_number(columns = c(NLL, Brier, AUC, ECE10), decimals = 3) %>%
  cols_label(
    model = "Modèle",
    NLL = "NLL (log-loss)",
    NLL_CI = "IC 95% NLL (Bayes)",
    Brier = "Brier",
    Brier_CI = "IC 95% Brier (Bayes)",
    AUC = "AUC",
    ECE10 = "ECE (10 bins)"
  ) %>%
  tab_header(
    title = "Comparaison Bayes vs fréquentiste (jeu de test)",
    subtitle = "NLL/Brier : plus petit = mieux ; AUC : plus grand = mieux ; ECE : plus petit = mieux."
  )

tab
