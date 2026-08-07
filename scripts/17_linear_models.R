##### Background #####
# Here I store scripts for running linear models on data

##### Initialise #####
library(mgcv)

# Extract coefficients (lm) or the parametric term table (gam) from a named list
# of models and format into a single table with term/model/significance columns
format_model_coefs <- function(modelList){
  modelList %>%
    imap(function(model, name){
      coefTab <- if (inherits(model, "gam")) summary(model)$p.table else summary(model)$coefficients
      as.data.frame(coefTab) %>%
        rownames_to_column("term") %>%
        mutate(model = name)
    }) %>%
    bind_rows() %>%
    mutate(termType = case_when(
      term == "(Intercept)" ~ "intercept",
      grepl(":", term) ~ "interaction",
      TRUE ~ "main effect"
    )) %>%
    mutate(stars = case_when(
      `Pr(>|t|)` < 0.001 ~ "***",
      `Pr(>|t|)` < 0.01 ~ "**",
      `Pr(>|t|)` < 0.05 ~ "*",
      TRUE ~ ""
    )) %>%
    relocate(model, term, termType)
}

