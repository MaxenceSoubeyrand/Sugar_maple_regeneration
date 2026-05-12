diagnostic_plots <- function(model, obs){
  #data
  diag_df <- data.frame(
    fitted = fitted(model),                         
    pearson = resid(model, type = "pearson"),
    deviance = resid(model, type = "deviance"),
    site = obs$site,
    block = obs$block
  )
  
  #Pearson vs fitted
  pd1 <- ggplot(diag_df, aes(fitted, pearson)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      x = "Adjusted probability",
      y = "Pearson residuals", 
      title="Pearson residuals vs fitted values"
    ) +
    theme_bw()
  
  #Pearson residuals QQ-plot
  pd2 <- ggplot(diag_df, aes(sample = pearson)) +
    stat_qq(alpha = 0.6) +
    stat_qq_line() +
    labs(x = "Theoritical quantiles",
         y = "Pearson residuals",
         title = "Pearson residuals QQ-plot") +
    theme_bw()
  
  #Residuals site 
  pd3 <- ggplot(diag_df, aes(site, pearson)) +
    geom_boxplot(outlier.alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      x = "Site",
      y = "Pearson residuals",
      title = "Site residuals"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90))
  
  #Radom effect qq plots
  ran_site <- ranef(model)$site[,1]
  
  pd4 <- ggplot(data.frame(effet = ran_site),
                aes(sample = effet)) +
    stat_qq() +
    stat_qq_line() +
    labs(x = "Theoritical quantiles",
         y = "Pearson residuals",
         title = "Random effect QQ-plots") +
    theme_bw()
  
  #Dharma
  library(DHARMa)
  
  sim_res <- simulateResiduals(model)
  
  df <- data.frame(
    res = sim_res$scaledResiduals,
    pred = sim_res$fittedPredictedResponse
  )
  
  qq_df <- data.frame(
    expected = sort(ppoints(length(df$res))),
    observed = sort(df$res)
  )
  
  #qqplot
  pd5 <- ggplot(qq_df, aes(expected, observed)) +
    geom_point(size = 1.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(
      title = "DHARMa QQ plot residuals",
      x = "Uniform quantiles",
      y = "DHARMa simulates residuals"
    ) +
    theme_classic()
  
  ks <- ks.test(df$res, "punif")$p.value
  
  pd5 <- pd5 + annotate(
    "text",
    x = 0.1, y = 0.8,
    label = paste0("KS test: p = ", round(ks, 3)),
    hjust = 0
  )
  
  #Dharma res vs pred
  pd6 <- ggplot(df, aes(pred, res)) +
    geom_point(shape = 1) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_hline(yintercept = 0.25, linetype = 2) +
    geom_hline(yintercept = 0.5, linetype = 2) +
    geom_hline(yintercept = 0.75, linetype = 2) +
    labs(
      title = "DHARMa residual vs. predicted",
      x = "Model predictions (rank transformed)",
      y = "DHARMa residual"
    ) +
    theme_classic()
  
  ggarrange(
    plotlist = list(pd1, pd2, pd3, pd4, pd5, pd6), ncol=2, nrow=3,
    labels = c("a", "b", "c", "d", "e")
  )
}

