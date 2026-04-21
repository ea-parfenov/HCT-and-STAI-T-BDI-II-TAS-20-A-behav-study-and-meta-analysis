# HCT and STAI-T ====
# Description: script to conduct meta-analysis investigating relationships between the HCT and STAI-T

# Create a vector of packages to use 
used_pachages <- c("dmetar", "ggplot2", "gridExtra", "metafor", "metameta", "metaviz", "readxl")
# Load packages
lapply(used_pachages, require, character.only = TRUE)

# Set working directory
setwd('/home/eparfenov/projects/HCT-questionnaires/')

# Load data and set up factors
mydata = read.table("Articles - STAI-T.tsv", sep = '\t', quote = "", header=TRUE)
mydata$Article <- as.factor(mydata$Article)
mydata$Formula_type <- as.factor(mydata$Formula_type)
mydata$Instructions_comb <- as.factor(mydata$Instructions_comb)
mydata$r_nature <- as.factor(mydata$r_nature)
# Set levels in the order they appear in the data set
mydata$Participants <- factor(mydata$Participants, levels = unique(mydata$Participants)) 

## Mixed effects model ====
res <- rma.mv(Corrected_r, Variance, random = ~ 1 | Article, data = mydata)
summary(res)

### Calculate heterogeneity ====
W <- diag(1/res$vi)
X <- model.matrix(res)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I = 100 * sum(res$sigma2) / (sum(res$sigma2) + (res$k-res$p)/sum(diag(P)))
I

## Influence analyses =====
### DFFITS ====
mydata$df_vals <- dfbetas(res)$intrcpt
which(abs(mydata$df_vals) > 3 * sqrt(1 / (length(unique(res$dat$ID)) - 1)))
### Cook's D ====
mydata$cook_vals <- cooks.distance(res)
which(mydata$cook_vals > 0.45)
### Hats ====
mydata$hat_vals <- hatvalues(res)
which(mydata$hat_vals > 3 * (1 / length(unique(res$dat$ID))))

### Plot influence analysis ====
p1 <- qplot(cook_vals, ID, data=mydata) + 
  geom_point(data = mydata[25, ], aes(color = "red")) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=8),
        axis.title.x=element_text(size=10),
        legend.position = "none") + 
  xlab("Cook's D")

p1 <- p1 + theme(plot.margin = unit(c(0.2, 0, 0, 1), "lines"))

p2 <- qplot(hat_vals, ID, data=mydata) + 
      geom_point(data = mydata[25, ], aes(color = "red")) +
      theme(axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      axis.text.x=element_text(size=8),
      axis.title.x=element_text(size=10),
      legend.position = "none") + 
  xlab("Hat value")

p2 <- p2 + theme(plot.margin = unit(c(0.2, 0, 0, 1), "lines"))

p3 <- qplot(df_vals, ID, data=mydata) + 
      geom_point(data = mydata[25, ], aes(color = "red")) +
      theme(axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      axis.text.x=element_text(size=8),
      axis.title.x=element_text(size=10),
      legend.position = "none") +
  xlab("DFFITS") + 
  scale_x_reverse()

p3 <- p3 + theme(plot.margin = unit(c(0.2, 0.2, 0, 1), "lines"))

p4 <- grid.arrange(p1, p2, p3, nrow = 1, widths=c(1.8,1.25,1.25))

ggsave('influence-analysis.png', path='plots/STAI', plot=p4, width=12, height=6, dpi=600)

# Remove influential cases
data_inf = mydata[-25, ]

## Recalculate mixed effects model ====
res <- rma.mv(Corrected_r, Variance, random = ~ 1 | Article, data = data_inf)
summary(res)
### Re-calculate heterogeneity ====
W <- diag(1/res$vi)
X <- model.matrix(res)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I_inf = 100 * sum(res$sigma2) / (sum(res$sigma2) + (res$k-res$p)/sum(diag(P)))
I_inf

## Forest plot ====
colour_palette <- c(
  "patients" = '#fc8d62',
  "healthy" = '#66c2a5'
)
# Assign population type to the specific colours
colout_colours <- colour_palette[res$data$Participants]

png(filename='plots/STAI/forest-plot_all.png', res=600, width=6, height=8, units='in')
forest_plot = forest(res, xlim = c(-5, 3.5),
                     showweights = TRUE,
                     at = seq(-1, 1), shade = "zebra",
                     ilab = (N_sub), ilab.xpos = (-2), cex = 0.5,
                     header = c("Study", "COR      [95% CI]"),
                     colout = colout_colours,
                     # col = col_colours,
                     mlab = "",
                     xlab = "",
                     slab = paste(data_inf$ID),
                     addcred = TRUE)
legend("bottomright", legend = names(colour_palette), fill = colour_palette, title = "Population type", xpd = TRUE, 
       inset = c(0.01, -0.05), 
       x.intersp = 0.5, # adjust horizontal interspacing
       y.intersp = 0.7, # adjust vertical interspacing
       cex = 0.55) # adjust font size
par(xpd = TRUE)
text(-2, res$k+2, "Total", cex = 0.5, font = 2)
text(-2, -1, paste(sum(data_inf$N_sub)), cex = 0.5, font = 2)
text(0, res$k+2, "Correlation", cex = 0.5, font = 2)
text(1.8, res$k+2, "Weight", cex = 0.5, font = 2)
text(forest_plot$xlim[1], -1.23, pos = 4,
     bquote(paste("Overall effect: Z = ", .(formatC(res$zval, digits=2, format="f")), 
                  ", p", .(ifelse(res$pval<.001, " < 0.001",
                                  paste0(" = ",formatC(res$pval, digits=3, format="f")))))),
     xpd = TRUE, cex = 0.5)
dev.off()

## Moderation analysis ====
### Population type ====
mod_part <- rma.mv(Corrected_r, Variance, random = ~ 1 | Article, mods = ~ Participants, data = data_inf,
                   tdist = TRUE # base test statistics and confidence intervals on the t-distribution
)
summary(mod_part)

# Plot and save
mod_partic <- ggplot(data_inf, aes(x = Participants, y = Corrected_r, size = 1/Variance)) + 
  geom_point(data = data_inf, shape = 1) + 
  xlab("") + ylab("Effect size") +
  geom_hline(yintercept = 0, linetype="dashed") +
  theme(axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 10, colour = 'black'),
        axis.text.x = element_text(size = 10, colour = 'black'),
        legend.position = "none",
        axis.line = element_line(color='black'),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_y_continuous(breaks = round(seq(min(data_inf$Corrected_r), max(data_inf$Corrected_r), by = 0.2), 1))
ggsave('participants-moderator.png', path='plots/STAI', plot=mod_partic, width=4, height=6, dpi=600)

### Instructions used ====
# Remove one study which didn't report type of instructions used
data_inf_instr <- data_inf[!is.na(data_inf$Instructions_comb),]
mod_instr <- rma.mv(Corrected_r, Variance, random = ~ 1 | Article, mods = ~ Instructions_comb, data = data_inf_instr, 
                    tdist = TRUE # base test statistics and confidence intervals on the t-distribution
)
summary(mod_instr)

# Check the effect each instructions type (model with no intercept)
# instr_type <- rma.mv(Corrected_r, Variance, random = ~ 1 | Article, mods = ~ Instructions_comb-1, data = data_inf_instr, 
#                     tdist = TRUE # base test statistics and confidence intervals on the t-distribution
# )
# summary(instr_type)

# Extract predicted values
preds <- predict(mod_instr, newmods = c(0,1))
# Plot and save
mod_instruct <- ggplot(data_inf_instr, aes(x = Instructions_comb, y = Corrected_r, size = 1/Variance))+
  geom_point(data = data_inf_instr, shape = 1) +
  xlab("") + ylab("") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_abline(intercept = preds$pred[1], slope = preds$pred[2], color = 'red', linewidth = 1) +
  theme(axis.text.y = element_text(size = 10, colour = 'black'),
        axis.text.x = element_text(size = 10, colour = 'black'),
        legend.position = "none",
        axis.line = element_line(color='black'),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_y_continuous(breaks = round(seq(min(data_inf$Corrected_r), max(data_inf$Corrected_r), by = 0.2), 1)) +
  scale_x_discrete(labels = c("1" = "original", "2" = "modified"))
ggsave('instructions-moderator.png', path='plots/STAI', plot=mod_instruct, width=4, height=6, dpi=600)

### Formula used ====
mod_formula <- rma.mv(Corrected_r, Variance, random = ~ 1 | Article, mods = ~ Formula_type, data = data_inf,
                      tdist = TRUE # base test statistics and confidence intervals on the t-distribution
)
summary(mod_formula)

# Plot and save
mod_form <- ggplot(data_inf, aes(x = Formula_type, y = Corrected_r, size = 1/Variance)) + 
  geom_point(data = data_inf, shape = 1) + 
  xlab("") + ylab("") +
  geom_hline(yintercept = 0, linetype="dashed") +
  theme(axis.text.y = element_text(size = 10, colour = 'black'),
        axis.text.x = element_text(size = 10, colour = 'black'),
        legend.position = "none",
        axis.line = element_line(color='black'),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  scale_y_continuous(breaks = round(seq(min(data_inf$Corrected_r), max(data_inf$Corrected_r), by = 0.2), 1)) +
  scale_x_discrete(labels = c("0" = "original", "1" = "adapted"))
ggsave('formula-moderator.png', path='plots/STAI', plot=mod_form, width=4, height=6, dpi=600)

### Correction for multiple comparison ====
p_values <- c(mod_part$QMp, mod_instr$QMp, mod_formula$QMp)
p.adjust(p_values, method = "bonferroni")

## Publication bias ====
### Egger's test (Funnel plot asymmetry) ====
reg <- regtest(Corrected_r, Variance, data=data_inf)
reg

### Funnel plot ====
png(filename='plots/STAI/funnel_all.png', res=600, width=6, height=6, units='in')
funnel_plot = funnel(res,
                     level=c(90, 95, 99),
                     xlab = "Observed correlation",
                     shade=c("white", "gray55", "gray75"),
                     legend=TRUE)
se <- seq(0,1.8,length=100) # add regression line to funnel plot
lines(coef(reg$fit)[1] + coef(reg$fit)[2]*se, se, lwd=3)
legend("bottomright", inset = 0.02, 
       lty = "solid", lwd = 3, 
       cex = 0.5, 
       bg = "white", 
       legend = "Standard Errors as Predictor")
dev.off()

## Study-level power analysis ====
data_inf$sei <- sqrt(data_inf$Variance)
data_inf$yi <- data_inf$Corrected_r
power_parf <- mapower_se(dat=data_inf, 
                         observed_es = 0.06, # observed effect size from the model
                         name = 'Parfenov & Duncan 2023')
power_pd_dat <- power_parf$dat
power_parf$power_median_dat
min(power_pd_dat$power_es_observed)
max(power_pd_dat$power_es_observed)
median(power_pd_dat$power_es_observed)
### Sunset plot ====
sunset_plot <- viz_sunset(res, true_effect = 0.06, sig_level = 0.05, power_contours = "continuous", xlab = "Observed correlation")
ggsave('plots/STAI/sunset.png', plot = sunset_plot, dpi = 600, width = 6, height = 6, units = 'in')
