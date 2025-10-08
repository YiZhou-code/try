kernel_estimator <- function(samp_emp,alpha,beta,len){
  # Normal scale bandwidth
  (Hns <- ks::Hns(x = samp_emp))
  # 
  # kde_Hns <- ks::kde(x = samp_emp, H = Hns)
  
  # 自定义更密集的网格
  # x_grid <- seq(min(samp_emp[, 1]), max(samp_emp[, 1]), length.out = len)  # 100 个点
  # y_grid <- seq(min(samp_emp[, 2]), max(samp_emp[, 2]), length.out = len)  # 100 个点
  x_grid <- seq(min(samp_emp[, 1])-0.05, max(samp_emp[, 1])+0.05, length.out = len)  # 100 个点
  y_grid <- seq(min(samp_emp[, 1])-0.05, max(samp_emp[, 1])+0.05, length.out = len)  # 100 个点
  
  # 计算网格间的步长
  dx <- diff(x_grid)[1]
  dy <- diff(y_grid)[1]
  
  # 创建 eval.points 数据框
  eval_points <- expand.grid(x = x_grid, y = y_grid)
  
  # 使用更密集的网格计算核密度估计
  kde_Hns <- ks::kde(x = samp_emp, H = Hns, eval.points = eval_points)
  
  # 获取联合密度估计值 f(x, y)
  density_values <- kde_Hns$estimate
  grid <- kde_Hns$eval.points
  
  # 将 density_values 转换为矩阵
  density_matrix <- matrix(density_values, nrow = length(x_grid), ncol = length(y_grid), byrow = TRUE)
  
  # 将 distribution_values 转换为矩阵
  distribution_matrix <- matrix(density_values, nrow = length(x_grid), ncol = length(y_grid), byrow = TRUE)*dx*dy
  
  distribution_cal <- function(x_up,y_up){
    # 找到小于 cv 和 v 的 x 和 y 的索引
    x_idx <- which(x_grid <= x_up)
    y_idx <- which(y_grid <= y_up)
    # 提取所需区域的密度值
    subset_density_values <- distribution_matrix[x_idx, y_idx]
    # 计算积分
    integral_result <- sum(subset_density_values) 
    integral_result
  }
  
  var_cal <- function(x_up){
    # 找到小于 cv 和 v 的 x 和 y 的索引
    x_idx <- which(x_grid <= x_up)
    y_idx <- which(y_grid <= Inf)
    # 提取所需区域的密度值
    subset_density_values <- distribution_matrix[x_idx, y_idx]
    # 计算积分
    integral_result <- sum(subset_density_values) 
    diff_result <- (integral_result-alpha)^2
    diff_result
  }
  
  # var_cal(-0.00798)
  
  # 使用 optimize 来寻找最小化 var_cal 的 x_up 值
  var_kernel <- optimize(var_cal, interval = c(min(x_grid), max(x_grid)))$minimum
  
  covar_cal <- function(y_up){
    # 找到小于 cv 和 v 的 x 和 y 的索引
    x_idx <- which(x_grid <= var_kernel)
    y_idx <- which(y_grid <= y_up)
    # 提取所需区域的密度值
    subset_density_values <- distribution_matrix[x_idx, y_idx]
    # 计算积分
    integral_result <- sum(subset_density_values) 
    diff_result <- (integral_result-alpha*beta)^2
    diff_result
  }
  
  # 使用 optimize 来寻找最小化 covar_cal 的 x_up 值
  covar_kernel <- optimize(covar_cal, interval = c(min(y_grid), max(y_grid)))$minimum
  
  # 生成一个全是 1 的向量
  ones_vector <- rep(1, length(y_grid))
  y_grid_matrix <- t(ones_vector %*% t(y_grid))
  
  # 将 distribution_values 转换为矩阵
  distribution_matrix_y <- y_grid_matrix*distribution_matrix
  
  x_idx <- which(x_grid <= var_kernel)
  y_idx <- which(y_grid <= Inf)
  # 提取所需区域的密度值
  subset_density_values <- distribution_matrix_y[x_idx, y_idx]
  # 计算积分
  mes_kernel <- sum(subset_density_values) / alpha
  ##
  return(c(var_kernel, covar_kernel, mes_kernel))
}