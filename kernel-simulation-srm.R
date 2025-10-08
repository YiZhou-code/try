library(MASS)
library(ks)
# 调用函数
source('C:/Users/Administrator/Desktop/JEF-new method/核估计方法进行估计/核估计系统性风险测度模拟/Kernel-in-sample/kernel_estimator.R')

# 定义均值向量和协方差矩阵
mean_vector <- c(0, 0)  # 二元正态分布的均值，两个变量均为0
cov_matrix <- matrix(c(0.1603, 0.0583, 0.0583, 0.0904), nrow = 2)  # 协方差矩阵

# true var covar mes
samp_true <- mvrnorm(100000, mu = mean_vector, Sigma = cov_matrix)
# 置信水平
alpha <- 0.05
beta <- 0.05
var_true <- quantile(samp_true[,1], probs = alpha)
sel_data_true <- samp_true[samp_true[, 1] < var_true, ]
covar_true <- quantile(sel_data_true[, 2], probs = beta)
mes_true <- mean(sel_data_true[, 2])

simu_num <- 100
result_matrix <- matrix(1:(simu_num*6), nrow = simu_num, ncol = 6)
colnames(result_matrix) <- c("var_emp", "covar_emp", "mes_emp", "var_kernel", "covar_kernel", "mes_kernel")

for (i in 1:simu_num) {# 生成100个二元正态分布的随机数
  i<-1
  n <- 1439  # 随机数的个数
  samp_emp <- mvrnorm(n, mu = mean_vector, Sigma = cov_matrix)
  # emprical var covar mes
  var_emp <- quantile(samp_emp[,1], probs = alpha)
  sel_data_emp <- samp_emp[samp_emp[, 1] <= var_true, ]
  covar_emp <- quantile(sel_data_emp[, 2], probs = beta)
  mes_emp <- mean(sel_data_emp[, 2])
  # Normal scale bandwidth
  mean_samp_emp_1 <- mean(samp_emp[,1])
  mean_samp_emp_2 <- mean(samp_emp[,2])
  std_samp_emp_1 <- var(samp_emp[,1])^(0.5)
  std_samp_emp_2 <- var(samp_emp[,2])^(0.5)
  samp_emp_1=(samp_emp[,1]-mean_samp_emp_1)/std_samp_emp_1
  samp_emp_2=(samp_emp[,2]-mean_samp_emp_2)/std_samp_emp_2
  samp_emp_stand <-  data.frame(Column1 = samp_emp_1, Column2 = samp_emp_2)
  ###
  kernel_est <- kernel_estimator(samp_emp_stand,alpha,beta,1440)
  
  var_kernel <- kernel_est[1]*std_samp_emp_1+mean_samp_emp_1
  covar_kernel <- kernel_est[2]*std_samp_emp_2+mean_samp_emp_2
  mes_kernel <- kernel_est[3]*std_samp_emp_2+mean_samp_emp_2
  # (Hns <- ks::Hns(x = samp_emp))
  # 
  # kde_Hns <- ks::kde(x = samp_emp, H = Hns)
  # 
  # # 自定义更密集的网格
  # x_grid <- seq(min(samp_emp[, 1]), max(samp_emp[, 1]), length.out = 1500)  # 100 个点
  # y_grid <- seq(min(samp_emp[, 2]), max(samp_emp[, 2]), length.out = 1500)  # 100 个点
  # 
  # # 计算网格间的步长
  # dx <- diff(x_grid)[1]
  # dy <- diff(y_grid)[1]
  # 
  # # 创建 eval.points 数据框
  # eval_points <- expand.grid(x = x_grid, y = y_grid)
  # 
  # # 使用更密集的网格计算核密度估计
  # kde_Hns <- ks::kde(x = samp_emp, H = Hns, eval.points = eval_points)
  # 
  # # 获取联合密度估计值 f(x, y)
  # density_values <- kde_Hns$estimate
  # grid <- kde_Hns$eval.points
  # 
  # # 将 density_values 转换为矩阵
  # density_matrix <- matrix(density_values, nrow = length(x_grid), ncol = length(y_grid), byrow = TRUE)
  # 
  # # 将 distribution_values 转换为矩阵
  # distribution_matrix <- matrix(density_values, nrow = length(x_grid), ncol = length(y_grid), byrow = TRUE)*dx*dy
  # 
  # distribution_cal <- function(x_up,y_up){
  #   # 找到小于 cv 和 v 的 x 和 y 的索引
  #   x_idx <- which(x_grid <= x_up)
  #   y_idx <- which(y_grid <= y_up)
  #   # 提取所需区域的密度值
  #   subset_density_values <- distribution_matrix[x_idx, y_idx]
  #   # 计算积分
  #   integral_result <- sum(subset_density_values) 
  #   integral_result
  # }
  # 
  # var_cal <- function(x_up){
  #   # 找到小于 cv 和 v 的 x 和 y 的索引
  #   x_idx <- which(x_grid <= x_up)
  #   y_idx <- which(y_grid <= Inf)
  #   # 提取所需区域的密度值
  #   subset_density_values <- distribution_matrix[x_idx, y_idx]
  #   # 计算积分
  #   integral_result <- sum(subset_density_values) 
  #   diff_result <- (integral_result-alpha)^2
  #   diff_result
  # }
  # 
  # # 使用 optimize 来寻找最小化 var_cal 的 x_up 值
  # var_kernel <- optimize(var_cal, interval = c(min(x_grid), max(x_grid)))$minimum
  # var_kernel <- var_kernel/10000
  # 
  # covar_cal <- function(y_up){
  #   # 找到小于 cv 和 v 的 x 和 y 的索引
  #   x_idx <- which(x_grid <= var_kernel)
  #   y_idx <- which(y_grid <= y_up)
  #   # 提取所需区域的密度值
  #   subset_density_values <- distribution_matrix[x_idx, y_idx]
  #   # 计算积分
  #   integral_result <- sum(subset_density_values) 
  #   diff_result <- (integral_result-alpha*beta)^2
  #   diff_result
  # }
  # 
  # # 使用 optimize 来寻找最小化 covar_cal 的 x_up 值
  # covar_kernel <- optimize(covar_cal, interval = c(min(y_grid), max(y_grid)))$minimum
  # covar_kernel <- covar_kernel/10000
  # 
  # # 生成一个全是 1 的向量
  # ones_vector <- rep(1, length(y_grid))
  # y_grid_matrix <- t(ones_vector %*% t(y_grid))
  # 
  # # 将 distribution_values 转换为矩阵
  # distribution_matrix_y <- y_grid_matrix*distribution_matrix
  # 
  # x_idx <- which(x_grid <= var_kernel)
  # y_idx <- which(y_grid <= Inf)
  # # 提取所需区域的密度值
  # subset_density_values <- distribution_matrix_y[x_idx, y_idx]
  # # 计算积分
  # mes_kernel <- sum(subset_density_values) / alpha
  # mes_kernel <- mes_kernel/10000 
  ##输出结果
  result_matrix[i, 1] <- var_emp
  result_matrix[i, 2] <- covar_emp
  result_matrix[i, 3] <- mes_emp
  result_matrix[i, 4] <- var_kernel
  result_matrix[i, 5] <- covar_kernel
  result_matrix[i, 6] <- mes_kernel
  
  print(paste("i=", i))
}

rmse_var_emp <-  sqrt(mean((result_matrix[,1]-var_true)^2))
rmse_covar_emp <- sqrt(mean(result_matrix[,2]-covar_true)^2)
rmse_mes_emp <-  sqrt(mean(result_matrix[,3]-mes_true)^2)
rmse_var_kernel <-  sqrt(mean((result_matrix[,4]-var_true)^2))
rmse_covar_kernel <- sqrt(mean(result_matrix[,5]-covar_true)^2)
rmse_mes_kernel <-  sqrt(mean(result_matrix[,6]-mes_true)^2)

