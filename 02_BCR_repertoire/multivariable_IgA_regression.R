# =============================================================================
# Multivariable analysis of clonal IgA proportion
# =============================================================================

library(dplyr)
library(broom)
library(ggplot2)

# ---- 1. Multiple linear regression ------------------------------------------
model <- lm(
  iga_fraction ~ age + recurrent_tonsillitis + sdb + sex +
                 total_ige + eosinophil_pct + tonsil_size,
  data = dat
)
print(summary(model))

# ---- 2. Coefficient table for the forest plot -------------------------------
coefs <- tidy(model, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(significant = p.value < 0.05)

# ---- 3. Forest plot (beta +/- 95% CI) ---------------------------------------
ggplot(coefs, aes(x = estimate, y = reorder(term, estimate))) +
  geom_vline(xintercept = 0, linetype = "longdash", color = "grey55") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = significant),
                 height = 0.25) +
  geom_point(aes(color = significant), size = 3) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#1B4F9E")) +
  labs(x = "Effect size (beta, 95% CI)", y = NULL,
       title = "Clinical predictors of clonal IgA proportion") +
  theme_classic(base_size = 13)
ggsave("Fig_IgA_Clonal_forest.pdf", width = 7, height = 5, device = cairo_pdf)
