

library("L0Learn")
library("tidyr")
library("Seurat")
source("as_matrix_cpp.R")

#
load("/data/mengxu/data/L0/lung_L0_data_harmony.Rdata")
seu_list <- SplitObject(scRNA_harmony, split.by = "stage")

# normal -----------------------------------------------------------------------
mat_normal <- seu_list[[1]]@assays$SCT@counts %>%
  as.matrix() %>%
  as.data.frame()
genes_normal <- rownames(mat_normal)
samples_normal <- colnames(mat_normal) %>% as.data.frame()
samples_normal$label <- "0"

# Stage-1 -----------------------------------------------------------------------
mat_stage1 <- seu_list[[2]]@assays$SCT@counts %>%
  as.matrix() %>%
  as.data.frame()
genes_stage1 <- rownames(mat_stage1)
samples_stage1 <- colnames(mat_stage1) %>% as.data.frame()
samples_stage1$label <- "1"

# Stage-2 -----------------------------------------------------------------------
mat_stage2 <- seu_list[[3]]@assays$SCT@counts %>%
  as.matrix() %>%
  as.data.frame()
genes_stage2 <- rownames(mat_stage2)
samples_stage2 <- colnames(mat_stage2) %>% as.data.frame()
samples_stage2$label <- "2"

# Stage-3 -----------------------------------------------------------------------
mat_stage3 <- seu_list[[4]]@assays$SCT@counts %>%
  as.matrix() %>%
  as.data.frame()
genes_stage3 <- rownames(mat_stage3)
samples_stage3 <- colnames(mat_stage3) %>% as.data.frame()
samples_stage3$label <- "3"

# Stage-4 -----------------------------------------------------------------------
mat_stage4 <- seu_list[[5]]@assays$SCT@counts %>%
  as.matrix() %>%
  as.data.frame()
genes_stage4 <- rownames(mat_stage4)
samples_stage4 <- colnames(mat_stage4) %>% as.data.frame()
samples_stage4$label <- "4"


# Combine -----------------------------------------------------------------
samples_combine <- rbind(
  samples_normal,
  samples_stage1,
  samples_stage2,
  samples_stage3,
  samples_stage4
)

mat_combine <- rbind(
  t(mat_normal),
  t(mat_stage1),
  t(mat_stage2),
  t(mat_stage3),
  t(mat_stage4)
)


# -------------------------------------------------------------------------

# data_other_gene<-data_other_gene[samples_com$.,]#只保留具有标签的样本
# data_other_gene<-data.frame(data_other_gene,check.names=TRUE)#check.names是因为回归变量名不能有‘-’，比如找不到对象'HLA-B'

# X_SNV=as.matrix(mat_com)
X_SNV <- mat_combine
row.names(X_SNV) <- row.names(mat_combine)

Y_label <- samples_combine$label %>%
  as.numeric() %>%
  as.vector()
names(Y_label) <- samples_combine$.


# clear workspace ---------------------------------------------------------

rm(scRNA_harmony)
rm(seu_list)

rm(mat_normal)
rm(mat_stage1)
rm(mat_stage2)
rm(mat_stage3)
rm(mat_stage4)

rm(samples_normal)
rm(samples_stage1)
rm(samples_stage2)
rm(samples_stage3)
rm(samples_stage4)

rm(mat_combine)
rm(samples_combine)

rm(genes_normal)
rm(genes_stage1)
rm(genes_stage2)
rm(genes_stage3)
rm(genes_stage4)
gc()

# Save --------------------------------------------------------------------
save(X_SNV, Y_label, file = paste0("/data/mengxu/data/L0/lung_L0_input_data.Rdata"))
# Load
load(paste0("/data/mengxu/data/L0/lung_L0_input_data.Rdata"))

# Model -------------------------------------------------------------------
maxSNVSize <- 50

# L0L2 --------------------------------------------------------------------

print("L0 learn")
cvfit <- L0Learn.cvfit(X_SNV, Y_label,
  penalty = "L0L2", nGamma = 5, gammaMin = 0.0001,
  gammaMax = 10, maxSuppSize = maxSNVSize
)

print("L0 done")
j <- which.min(lapply(cvfit$cvMeans, min))
optimalGammaIndex <- j # index of the optimal gamma identified previously
optimalLambdaIndex <- which.min(cvfit$cvMeans[[optimalGammaIndex]])
optimalLambda <- cvfit$fit$lambda[[optimalGammaIndex]][optimalLambdaIndex]
print("predicting...")
y_cat <- predict(cvfit,
  newx = X_SNV, # predict函数要求行为样本。所以需转置
  lambda = optimalLambda, gamma = cvfit$fit$gamma[j]
)
print("predict done")
y_hat <- as.vector(y_cat)

library(Metrics)
res_rmse <- rmse(Y_label, y_hat)
res_rse <- rse(Y_label, y_hat)
r_square <- 1 - res_rse
print("rmse done")
############
temp <- coef(cvfit, lambda = optimalLambda, gamma = cvfit$fit$gamma[j])
temp <- as.vector(temp)
temp <- temp[-1] # 排除第一个位置上的intercept
temp <- which(temp != 0) # 排除系数为0的冗余特征
temp <- colnames(X_SNV)[temp]
X_Y <- cbind(X_SNV[, temp], Y_label)
X_Y_frame <- as.data.frame(X_Y)
lmfit <- lm(Y_label ~ ., data = X_Y_frame)
fit_temp <- summary(lmfit)
L0_Rsquare_L0L2 <- 1 - mse(Y_label, y_hat) / var(Y_label)
########
write.csv(fit_temp$coefficients, file = "/data/mengxu/data/L0/feature_selected_byL0L2.csv")

# L0L1 --------------------------------------------------------------------


print("L0 learn")
cvfit <- L0Learn.cvfit(X_SNV, Y_label,
  penalty = "L0L1", nGamma = 5, gammaMin = 0.0001,
  gammaMax = 10, maxSuppSize = maxSNVSize
)

# Extract coefficient at middle lambda
fit_L0_information <- as.data.frame(print(fit_L0))
fit_L0_information <- fit_L0_information[order(fit_L0_information$suppSize, decreasing = TRUE), ]
lambda_L0 <- fit_L0_information$lambda[1]
gamma_L0 <- fit_L0_information$gamma[1]

print("L0 done")
j <- which.min(lapply(cvfit$cvMeans, min))
optimalGammaIndex <- j # index of the optimal gamma identified previously
optimalLambdaIndex <- which.min(cvfit$cvMeans[[optimalGammaIndex]])
optimalLambda <- cvfit$fit$lambda[[optimalGammaIndex]][optimalLambdaIndex]
print("predicting...")
y_cat <- predict(cvfit,
  newx = X_SNV, # predict函数要求行为样本。所以需转置
  lambda = optimalLambda, gamma = cvfit$fit$gamma[j]
)
print("predict done")
y_hat <- as.vector(y_cat)

library(Metrics)
res_rmse <- rmse(Y_label, y_hat)
res_rse <- rse(Y_label, y_hat)
r_square <- 1 - res_rse
print("rmse done")
############
temp <- coef(cvfit, lambda = optimalLambda, gamma = cvfit$fit$gamma[j])
temp <- as.vector(temp)
temp <- temp[-1] # 排除第一个位置上的intercept
temp <- which(temp != 0) # 排除系数为0的冗余特征
temp <- colnames(X_SNV)[temp]
X_Y <- cbind(X_SNV[, temp], Y_label)
X_Y_frame <- as.data.frame(X_Y)
lmfit <- lm(Y_label ~ ., data = X_Y_frame)
fit_temp <- summary(lmfit)

L0_Rsquare_L0L1 <- 1 - mse(Y_label, y_hat) / var(Y_label)
########
write.csv(fit_temp$coefficients, file = "/data/mengxu/data/L0/feature_selected_byL0L1.csv")

# L0 -----------------------------------------------------------------------

cvfit <- L0Learn.fit(X_SNV, Y_label,
  penalty = "L0",
  maxSuppSize = maxSNVSize
)

print(cvfit)
# Extract the coefficients at middle lambda
lambda_value <- as.data.frame(cvfit$lambda)
names(lambda_value) <- "lambda"
lambda_value[nrow(lambda_value) / 2, ]

y_cat <- predict(cvfit,
  newx = X_SNV,
  lambda = lambda_value[nrow(lambda_value) / 2, ],
  gamma = 0
)

y_hat <- as.vector(y_cat)
plot(cvfit, gamma = cvfit$gamma, showLines = TRUE)

temp <- coef(cvfit,
  lambda = lambda_value[nrow(lambda_value) / 2, ],
  gamma = cvfit$gamma
)
temp <- as.vector(temp)
temp <- temp[-1] # 排除第一个位置上的intercept
temp <- which(temp != 0) # 排除系数为0的冗余特征
temp <- colnames(X_SNV)[temp]
temp
# 这里的ifelse是为了处理只有一个样本TFs被选择，列名为“V1”的问题
if (length(temp) == 1) {
  X_Y <- cbind(X_SNV[, temp], Y_label)
  colnames(X_Y)[1] <- temp
  X_Y_frame <- as.data.frame(X_Y)
} else {
  X_Y <- cbind(X_SNV[, temp], Y_label)
  X_Y_frame <- as.data.frame(X_Y)
}
lmfit <- lm(Y_label ~ ., data = X_Y_frame)
fit_temp <- summary(lmfit)
res_data <- as.matrix(fit_temp$coefficients)
res_data

#### 打印模型信息####
print(tidy(lmfit))


L0_Rsquare_L0 <- 1 - mse(Y_label, y_hat) / var(Y_label)
# ####评价指标####
# evaluate_L0 <- data.frame(L0_Rsquare = 1-mse(Y_label,y_hat)/var(Y_label), L0_RMSE = RMSE(Y_label,y_hat))
#
# ####只保存P小于0.05的结果####
# res_data_f <- as.data.frame(res_data)
# res_data_f <- res_data_f[which(res_data_f$`Pr(>|t|)` <= 0.05),]
# write.csv(res_data_f ,file='/data/mengxu/data/L0/feature_selected_byL0.csv')


fit_temp <- summary(lmfit)
########
write.csv(fit_temp$coefficients, file = "/data/mengxu/data/L0/feature_selected_byL0.csv")
gc()
