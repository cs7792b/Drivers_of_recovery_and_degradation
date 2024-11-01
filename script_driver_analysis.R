#### analysis  paper 1  ####

#loading data
setwd("C:/Users/cschu/Documents/Paper1/")
data <- read.csv("table_driver_analysis.csv", sep = ",", header = T, stringsAsFactors = FALSE)

#loading libraries
library(dplyr)
library(ggplot2)
library(vegan)
library(missForest)
library(lme4)
library(MASS)
library(lmerTest)
library(randomForest)
library(glmmTMB)
library(cowplot)
library(MuMIn)
library(DHARMa)


#### start analysis ####


#### pca ####
colnames(data)

pca_result <- prcomp(data[, c(8:15)], center = TRUE, scale. = TRUE)

# Create a data frame for PCA results
pca_data <- as.data.frame(pca_result$x)
pca_data$Ecological_Status <- data$OEZK 

loadings <- as.data.frame(pca_result$rotation)
loadings$variable <- rownames(loadings)

# Plot the PCA results using ggplot2
p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Ecological_Status)) +
  geom_point(size = 5, alpha=0.8, pch=19) +
  scale_color_gradient(low = "red", high = "green") +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1*5, yend = PC2*5), 
               arrow = arrow(length = unit(0.3, "cm")), color = "black", size = 1) +
  geom_text(data = loadings, aes(x = PC1*5, y = PC2*5, label = variable), 
            vjust = 1.5, color = "black", size = 5) +
  xlim(-2, 5) +
  ylim(-5, 5) +
  labs(title = "PCA of Abiotic Stressors", x = "PC1", y = "PC2") +
  theme_minimal()

# Print the plot
print(p)


#### save 1st axis ####
pca_result$rotation
data$PCA1 <- pca_result$x[, 1]


#### random forest rel importance ####


# Define relevant columns 
relevant_columns <- c("river_type", "delta_OEZK", "delta_gower_similarity", "delta_PCA", 
                      "System", "Status", "Date", "ID","Cat_Agri_without_grassl", "Cat_Urban","CatCropWald")

# Remove rows with missing values
data_naomit <- na.omit(data[relevant_columns])

set.seed(96)  # For reproducibility
rf_model <- randomForest(delta_OEZK ~ ., data = data_naomit, importance = TRUE)

rel_model <- rf_model$importance 
rel_model[,2] <- rel_model[,2]/sum(rel_model[,2])*100

importance_values <- rel_model[, 2]
variable_names <- c("River type", "Change in competition", "Change in abiotic stress", "River system", 
                    "River category", "Date", "Site ID", "Cropland", "Urban", "Forest")

# Combine into a data frame
importance_data <- data.frame(Variable = variable_names, Importance = importance_values)




importance_data



# Print the model summary to get % variance explained and MSE
print(rf_model)

# Access % variance explained (which approximates R²)
r_squared <- rf_model$rsq[length(rf_model$rsq)] * 100  

# Display model skill metrics
cat("Model Skill Metrics:\n")
cat("R-squared (% variance explained):", round(r_squared, 2), "%\n")



# Plot the importance values with a vertical line at 10
p <- ggplot(importance_data, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 5, color = "red", linetype = "dashed", size = 1.5) +  # Increase the thickness of the line
  coord_flip() +  # Flip coordinates for horizontal bars
  labs(x = "", y = "Relative importance (%)") +  # Removed title argument
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 15),  
    axis.text.y = element_text(size = 15),  
    axis.title = element_text(size = 16),   
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA)  
  )

ggsave("importance_plot.tiff", plot = p, dpi = 300, width = 8, height = 4, units = "in")



#### Relationship between change in Ecological status and Abiotic Stress and Competition ####


#### Fit the model using the entire dataset ####
mod <- glmmTMB(delta_OEZK ~ delta_PCA + (1 | ID), 
               family = gaussian(link = 'identity'), 
               data = data)

mf <- lm(resid(mod) ~ mod$frame$delta_OEZK)$fitted.values^2

mod_mf <- glmmTMB(delta_OEZK ~ delta_PCA + (1 | ID), 
                  family = gaussian(link = 'identity'), 
                  weights = mf,
                  data = data)

r.squaredGLMM(mod_mf)
r.squaredGLMM(mod)
qqnorm(resid(mod_mf))
qqline(resid(mod_mf), col='red')

# Summary of the model
summary(mod_mf)
summary(mod)



# Generate model predictions for plotting
xl <- seq(-5, 5, 0.01)
modlines <- fixef(mod_mf)[1]$cond[1] + fixef(mod_mf)[1]$cond[2] * xl
dfline <- data.frame(x = xl, y = modlines)

modlines <- fixef(mod)[1]$cond[1] + fixef(mod)[1]$cond[2] * xl
dfline2 <- data.frame(x = xl, y = modlines)

# Plot the data and the model predictions
p1 <- ggplot(data, aes(x = delta_PCA, y = delta_OEZK)) +
  geom_point(color = "gray28") +
  geom_line(data = dfline, aes(x = x, y = y), color = "red") +
  geom_line(data = dfline2, aes(x = x, y = y), color = "blue") +
  labs(title = "a) Relationship between Change in Ecological Status and Abiotic Stress",
       x = "Change in Abiotic Stress",
       y = "Change in Ecological Status") +
  xlim(-6, 6) +
  ylim(-0.6, 0.6)  +
  annotate("text", x = 4.9, y = 0.5, label = "weighted model: R²=0.39", hjust = 1, size = 5, color="red") +
  annotate("text", x = 4.9, y = 0.4, label = "unweighted model: R²=0.14", hjust = 1, size = 5, color="blue") +
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 16))

p1



# Fit the model using the entire dataset
mod <- glmmTMB(delta_OEZK ~ delta_gower_similarity + (1 | river_type / Date), 
               family = gaussian(link = 'identity'), 
               data = data)

mf <- lm(resid(mod) ~ mod$frame$delta_OEZK)$fitted.values^2

mod_mf <- glmmTMB(delta_OEZK ~ delta_gower_similarity + (1 | river_type / Date), 
                  family = gaussian(link = 'identity'), 
                  weights = mf,
                  data = data)

r.squaredGLMM(mod_mf)

qqnorm(resid(mod_mf))
qqline(resid(mod_mf), col='red')

# Summary of the model
summary(mod_mf)
summary(mod)


# Generate model predictions for plotting
xl <- seq(-5, 5, 0.01)
modlines <- fixef(mod_mf)[1]$cond[1] + fixef(mod_mf)[1]$cond[2] * xl
dfline <- data.frame(x = xl, y = modlines)

modlines <- fixef(mod)[1]$cond[1] + fixef(mod)[1]$cond[2] * xl
dfline2 <- data.frame(x = xl, y = modlines)

# Plot the data and the model predictions
p2 <- ggplot(data, aes(x = delta_gower_similarity, y = delta_OEZK)) +
  geom_point(color = "gray28") +
  geom_line(data = dfline, aes(x = x, y = y), color = "red") +
  geom_line(data = dfline2, aes(x = x, y = y), color = "blue") +
  labs(title = "b) Relationship between Change in Ecological Status and Competition",
       x = "Change in Competition",
       y = "") +
  xlim(-0.25, 0.25) +
  ylim(-0.6, 0.6)  +
  annotate("text", x = 0.25, y = -0.4, label = "weighted model: R²=0.12", hjust = 1, size = 5, color="red")+
  annotate("text", x = 0.25, y = -0.5, label = "unweighted model: R²=0.04", hjust = 1, size = 5, color="blue") +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 16))

p2


combined_plot <- cowplot::plot_grid(p1, p2, rel_widths = c(0.5, 0.5))

# Save the plot as a TIFF file with 300 DPI
tiff("combined_plot_full123.tiff", width = 18, height = 6, units = "in", res = 300)
print(combined_plot)
dev.off()




#### split the dataset in recovery and degradation ####

data_no_NA <- data %>%
  mutate(condition = factor(ifelse(delta_OEZK > 0, "Recovery", "Degradation"), levels = c("Recovery", "Degradation")))

df <- split(data_no_NA, data_no_NA$condition)
df$Degradation$delta_OEZK <- abs(df$Degradation$delta_OEZK)
mod1 <- glmmTMB(delta_OEZK~delta_PCA + (1|ID) , zi=~1,
                family = ziGamma(link = "log"), data = df$Degradation)

mod2 <- glmmTMB(delta_OEZK~delta_PCA+(1|ID) , zi=~1,
                family = ziGamma(link = "log"), data = df$Recovery)

summary(mod1)
summary(mod2)


lab1 <- paste0("Coefficient recovery=",round(fixef(mod2)[1]$cond[2],2),"\nCoefficient degradation=",round(-1*fixef(mod1)[1]$cond[2],2))

xl        <- seq(-5, 5, 0.01)
modlines1 <- -1*exp(fixef(mod1)[1]$cond[1]+fixef(mod1)[1]$cond[2]*xl)
modlines2 <- exp(fixef(mod2)[1]$cond[1]+fixef(mod2)[1]$cond[2]*xl)

dfline <- data.frame(x=c(xl, xl), y=c(modlines2, modlines1), condition=c(rep("Recovery", length(xl)), 
                                                                         rep("Degradation", length(xl))))

p1 <- ggplot(data_no_NA, aes(x = delta_PCA, y = delta_OEZK, color = condition)) +
  geom_point() +
  annotate("text", x=2.5, y=0.5, label=lab1, size=5)+
  geom_line(data=dfline, aes(x=x, y=y, group=condition), inherit.aes = F)+
  labs(title = "a) Relationship between change in Ecological Status and Abiotic Stress",
       x = "Change in Abiotic Stress",
       y = "Change in Ecological Status",
       color = "Condition") +
  xlim(-5.5,5.5)+
  ylim(-0.6, 0.6)+
  theme_minimal()+
  theme(legend.position = "none",
        text = element_text(size=14),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        plot.title = element_text(size=16))

mod1 <- glmmTMB(delta_OEZK~delta_gower_similarity+ (1|ID), zi=~1,
                family = ziGamma(link = "log"), data = df$Degradation)

mod2 <- glmmTMB(delta_OEZK~delta_gower_similarity+ (1|ID), zi=~1,
                family = ziGamma(link = "log"), data = df$Recovery)

summary(mod1)
summary(mod2)
lab2 <- paste0("Coefficient recovery=",round(fixef(mod2)[1]$cond[2],2),"\nCoefficient degradation=",round(-1*fixef(mod1)[1]$cond[2],2))

xl        <- seq(-1, 1, 0.01)
modlines1 <- -1*exp(fixef(mod1)[1]$cond[1]+fixef(mod1)[1]$cond[2]*xl)
modlines2 <- exp(fixef(mod2)[1]$cond[1]+fixef(mod2)[1]$cond[2]*xl)

dfline <- data.frame(x=c(xl, xl), y=c(modlines2, modlines1), condition=c(rep("Recovery", length(xl)), 
                                                                         rep("Degradation", length(xl))))

p2 <- ggplot(data_no_NA, aes(x = delta_gower_similarity, y = delta_OEZK, color = condition)) +
  geom_point() +
  geom_line(data=dfline, aes(x=x, y=y, group=condition), inherit.aes = F)+
  annotate("text", x=0.18, y=0.5, label=lab2, size=5)+
  labs(title = "b) Relationship between change in Ecological Status and Competition",
       x = "Change in Competition",
       y = "Change in Ecological Status",
       color = "Condition") +
  xlim(-0.25,.25)+
  ylim(-0.6, 0.6)+
  theme_minimal()+
  theme(axis.title.y = element_blank(),
        text = element_text(size=14),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        plot.title = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))

combined_plot <- cowplot::plot_grid(p1, p2, rel_widths = c(0.46, 0.54))

# Save the plot as a TIFF file with 300 DPI
tiff("combined_plot_abiotic_competition1.tiff", width = 18, height = 6, units = "in", res = 300)
print(combined_plot)
dev.off()



#### analyses land use ####
recovery_data <- data %>% filter(delta_OEZK > 0)
degradation_data <- data %>% filter(delta_OEZK < 0)
# For recovery data
cor_test_recovery_forest <- cor.test(recovery_data$CatCropWald, recovery_data$delta_OEZK, use="pairwise.complete.obs", method="spearman")
cor_test_recovery_urban <- cor.test(recovery_data$Cat_Urban, recovery_data$delta_OEZK, use="pairwise.complete.obs", method="spearman")
cor_test_recovery_agri <- cor.test(recovery_data$Cat_Agri_without_grassl, recovery_data$delta_OEZK, use="pairwise.complete.obs", method="spearman")

# For degradation data
cor_test_degradation_forest <- cor.test(degradation_data$CatCropWald, degradation_data$delta_OEZK, use="pairwise.complete.obs", method="spearman")
cor_test_degradation_urban <- cor.test(degradation_data$Cat_Urban, degradation_data$delta_OEZK, use="pairwise.complete.obs", method="spearman")
cor_test_degradation_agri <- cor.test(degradation_data$Cat_Agri_without_grassl, degradation_data$delta_OEZK, use="pairwise.complete.obs", method="spearman")

# For overall data
cor_test_overall_forest <- cor.test(data$CatCropWald, data$delta_OEZK, use="pairwise.complete.obs", method="spearman")
cor_test_overall_urban <- cor.test(data$Cat_Urban, data$delta_OEZK, use="pairwise.complete.obs", method="spearman")
cor_test_overall_agri <- cor.test(data$Cat_Agri_without_grassl, data$delta_OEZK, use="pairwise.complete.obs", method="spearman")




